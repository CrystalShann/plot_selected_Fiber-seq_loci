# selectFiberseqLociWorkflow
Analyze selected Fiber-seq loci, plot m6A and per-read footprint tracks, and
optionally cluster reads by their m6A patterns using Manhattan-distance Leiden
clustering.

## Install

Use an R version compatible with your Bioconductor release (R >= 4.4). Install
Bioconductor dependencies through BiocManager before installing from GitHub or
a source tarball:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("BiocGenerics", "GenomeInfoDb", "GenomicRanges",
    "IRanges", "S4Vectors", "SummarizedExperiment", "rtracklayer", "Rsamtools",
    "ComplexHeatmap"))
remotes::install_github("CrystalShann/plot_selected_Fiber-seq_loci")
```


```r
install.packages(c("Matrix", "data.table", "dplyr", "igraph", "withr", "ggplot2",
    "patchwork", "cowplot"))
install.packages("selectFiberseqLociWorkflow_0.99.0.tar.gz",
    repos = NULL, type = "source")
```


## Step-by-step analysis

Follow steps 1–10 in order. 
Replace the example paths and coordinates with your own data. The workflow inputs processed after Fibertools/FiberHMM, for all species and genome build.

### 1. Load the package

```r
library(selectFiberseqLociWorkflow)
library(SummarizedExperiment)
```

### 2. Prepare input files, paths and sample IDs

Use processed tracks from fibertools and/or FiberHMM, all aligned to the same
reference assembly. This package does not run these upstream pipelines. Every
supplied track must be coordinate-sorted, bgzip-compressed and tabix-indexed,
with an adjacent `.tbi` file (for example, `m6a.bed.gz.tbi`). Ordinary gzip
compression alone is insufficient.

| Track | Explicit format | Record interpretation |
| --- | --- | --- |
| m6A (required) | `bed12_fibertools` | Fibertools blocks; terminal sentinel blocks removed |
| Nucleosome (optional) | `bed12_fibertools` or `bed13_fiberhmm` | Fibertools sentinels removed; FiberHMM blocks retained |
| TF footprint (optional) | `bed15_fiberhmm_tf` or `bed6_per_sample` | FiberHMM per-block score threshold; BED6 already filtered |
| FIRE (optional) | `bed11_fire` | Per-read FIRE intervals |

Formats must match the actual file contents; they are not inferred from file
names or species. FiberHMM-only modification output cannot replace the required
fibertools m6A BED12 track. Footprint read IDs must match the original read IDs
in the corresponding sample's m6A file. No allele-specific RDS or phasing input
is required.

Create one row per sample. Here, each sample directory contains whole-genome
tracks:

```r
input_dir <- "/path/to/processed_tracks"  # change this
sample_ids <- c("sample_1", "sample_2")

sample_table <- data.frame(
    sample_name = sample_ids,
    sample_label = c("Control", "Treatment"),
    m6a = file.path(input_dir, sample_ids, "m6a.bed.gz"),
    m6a_format = "bed12_fibertools",
    nucleosome = file.path(input_dir, sample_ids, "nucleosome.bed.gz"),
    nucleosome_format = "bed13_fiberhmm",
    tf = file.path(input_dir, sample_ids, "tf.bed.gz"),
    tf_format = "bed15_fiberhmm_tf"
)
sample_table
```

`sample_name` is a unique ID using letters, numbers, underscores, hyphens or
periods (use `sample_1`, not `sample 1`). `sample_label` is a display/group
label; several samples can share a label. Omit optional track and format
columns when unavailable, or leave that sample's optional path blank.
To include FIRE, add `fire` paths and `fire_format = "bed11_fire"`.

Supported path styles include:

| Path in the sample sheet | Interpretation |
| --- | --- |
| `/data/sample_1/m6a.bed.gz` | Explicit absolute path to a whole-genome file |
| `tracks/sample_1/m6a.bed.gz` | Relative path; resolved against the TSV directory when loading a TSV, or `base_dir` for an in-memory table |
| `/data/{sample_name}/m6a.bed.gz` | One whole-genome file per sample |
| `/data/{sample_name}/m6a.{chr}.bed.gz` | One file per sample and contig |

Only `{sample_name}` and `{chr}` placeholders are supported. For example,
`/data/{sample_name}/m6a.{chr}.bed.gz` becomes
`/data/sample_1/m6a.Chr5.bed.gz` for sample `sample_1` and contig `Chr5`.
To use chromosome-split files, change the path columns before step 3:

```r
# Alternative paths; uncomment if your files follow this layout:
# sample_table$m6a <- file.path(input_dir, "{sample_name}", "m6a.{chr}.bed.gz")
# sample_table$nucleosome <- file.path(input_dir, "{sample_name}", "nucleosome.{chr}.bed.gz")
# sample_table$tf <- file.path(input_dir, "{sample_name}", "tf.{chr}.bed.gz")
```

**Test-data alternative:** instead of defining your own paths above, run:

```r
example <- example_data()  # copies bundled files into a writable directory
sample_table <- as.data.frame(example$samples)
sample_table[, c("sample_name", "sample_label")]
```

These are small subsets of the AL6/AL8 LCL tracks around ALAS1 and PRDX1,
not genome-wide data. Use the test coordinates in step 4 with these files.

### 3. Make the sample sheet

```r
samples <- make_sample_sheet(sample_table, output = "samples.tsv")
samples

# In a later R session, reload the saved sheet:
# samples <- make_sample_sheet("samples.tsv")
```

The returned object is an `S4Vectors::DataFrame`. The TSV is written in the
current working directory (`getwd()`). For relative paths in an in-memory
sample table, supply `base_dir = "/path/to/base"`; the default is `getwd()`.
File and contig checks are performed in step 5.

### 4. Prepare custom regions

Coordinate strings are **one-based, inclusive**. For example:

```r
regions <- prepare_regions(c("Chr5:1001-3000", "scaffold_42:4001-6000"))
names(regions) <- c("plant_locus_1", "plant_locus_2")
regions
```

Replace these illustrative plant coordinates with your own. Contig names must
match the tracks exactly, including case: `Chr5`, `chr5` and `5` are different.
All coordinates and annotations must use the same reference assembly as the
tracks. The returned object is a named `GRanges`.

You can also supply a **BED file**. For example, `regions.bed` can contain the
following tab-separated rows, with no header:

```text
Chr5	1000	3000	plant_locus_1
scaffold_42	4000	6000	plant_locus_2
```

BED uses **zero-based, half-open** coordinates: the first row above represents
one-based `Chr5:1001-3000`. Column 4 optionally supplies unique locus names.
Unlike input tracks, the region BED does not need bgzip compression or an index.

```r
# Alternative to coordinate strings:
# regions <- prepare_regions("regions.bed")
```

Coordinate tables require an explicit convention:

```r
# Equivalent table alternative:
# regions <- prepare_regions(
#     data.frame(region_id = "plant_locus_1", chr = "Chr5",
#                start = 1001, end = 3000),
#     coordinate_system = "one_based"
# )
# For a table containing BED coordinates, use coordinate_system = "bed".
```

To select promoters from a GTF/GFF3 annotation instead:

```r
# regions <- prepare_regions(
#     annotation = "/path/to/plant.gff3",
#     genes = c("YOUR_GENE_ID"), flank_bp = 1000
# )
```

Annotations must contain transcript/mRNA rows with gene relationships.
Canonical tags are not required. Distinct TSSs are retained; shared TSSs retain
all associated transcript IDs. The default promoter is one-based
`[TSS - 1000, TSS + 999]`; promoters extending before position 1 are rejected.
`extract_canonical_tss()` is available for explicitly selecting canonical
transcripts with a chosen annotation attribute/value.

**Test-data alternative:** use a custom window inside the bundled ALAS1 locus:

```r
regions <- prepare_regions("chr3:52197660-52198659")
names(regions) <- "ALAS1_custom"
# Or use both original example loci: regions <- example$regions
```

Choose one region-input option before continuing.

### 5. Validate the inputs

```r
validate_inputs(samples, regions)
```

This checks supplied paths, adjacent `.tbi` indexes, declared column counts and
m6A contig availability. Successful validation returns `TRUE` invisibly.
Missing optional files are errors if their paths were supplied; optional
tracks may legitimately have no calls on a selected contig.

### 6. Specify analysis parameters

Define settings once so both runs use the same footprint filters and profile
smoothing. These example values are the package defaults:

```r
window_size <- 0
k_neighbors <- 10
resolution <- 1
seed <- 1
nuc_min <- 130
nuc_max <- 160
min_tf_score <- 50
tf_length_ranges <- list(c(10, 29), c(40, 59), c(60, 79))
smoothing_bp <- 10
```

| Parameter | Meaning |
| --- | --- |
| `clustering` | Required `FALSE` or `TRUE`, chosen separately in steps 7 and 8. Clustering uses m6A patterns, not footprint tracks. |
| `window_size` | Clustering bin width in bp. `0` uses individual observed m6A sites; a positive integer, such as `50`, averages each read's calls over observed sites in each occupied bin. Bins start at the locus's first base. |
| `k_neighbors` | Positive integer number of Manhattan nearest neighbors per read; capped at retained reads minus one. |
| `resolution` | Nonnegative Leiden resolution; larger values generally favor finer partitions. |
| `seed` | Nonnegative integer random seed for reproducible Leiden clustering. |
| `nuc_min`, `nuc_max` | Inclusive nucleosome length limits in bp. Adjust these to your footprint definitions; defaults retain 130–160 bp. |
| `min_tf_score` | Minimum FiberHMM BED15 TF block score, inclusive. Custom BED6 TF inputs are assumed to be prefiltered and do not use this threshold. |
| `tf_length_ranges` | Nonempty list of nonoverlapping inclusive length ranges in bp, e.g. `list(c(30, 39), c(80, 120))`. Use `list(c(35, 35))` for exactly 35 bp or `NULL` for all positive lengths. Applies to both TF formats; defaults preserve 10–29, 40–59 and 60–79 bp. |
| `smoothing_bp` | Positive integer smoothing width in bp for accessibility profiles; independent of clustering bins. |

Binning changes clustering features; the stored m6A assay remains single-base.
The clustering settings do not change the pooled m6A profile. For an illustrative
binned run, set `window_size <- 50`, `k_neighbors <- 5`, and
`resolution <- 0.5` before running step 8. These values are not optimized for a
particular species or experiment.

### 7. Get m6A and footprint results without clustering

```r
unclustered_results <- analyze_loci(
    samples, regions,
    clustering = FALSE,
    nuc_min = nuc_min, nuc_max = nuc_max,
    min_tf_score = min_tf_score, smoothing_bp = smoothing_bp,
    tf_length_ranges = tf_length_ranges
)

locus_id <- names(regions)[1]  # choose a locus to inspect
x <- unclustered_results[[locus_id]]
x

m6a_matrix <- assay(x, "m6a")       # genomic positions x reads
read_metadata <- colData(x)
footprints <- colData(x)$footprints # per-read GRangesList
aggregate_m6a <- locus_profile(x)   # pooled profile across retained reads
head(aggregate_m6a)
S4Vectors::metadata(x)$qc
```

Each selected region is analyzed independently. Only the longest overlapping
alignment per read is considered, and it must span the entire locus. Zero-call
reads are retained and included in the pooled denominator. Sparse zeros mean
no observed m6A call, not a measured unmethylated adenine at every base.

The objects above supply the m6A heatmap, aggregate m6A profile and footprint
panels plotted in step 9. `frac` in the profile is the observed-call fraction;
`smooth_frac` is its smoothed value.

### 8. Get clustering results

```r
clustered_results <- analyze_loci(
    samples, regions,
    clustering = TRUE,
    window_size = window_size, k_neighbors = k_neighbors,
    resolution = resolution, seed = seed,
    nuc_min = nuc_min, nuc_max = nuc_max,
    min_tf_score = min_tf_score, smoothing_bp = smoothing_bp,
    tf_length_ranges = tf_length_ranges
)

xc <- clustered_results[[locus_id]]
clustering_status <- S4Vectors::metadata(xc)$clustering_status
clustering_status

if (identical(clustering_status, "complete")) {
    assignments <- colData(xc)[, c("RID", "sample_name", "cluster")]
    print(assignments)
    print(table(sample = colData(xc)$sample_label,
                cluster = colData(xc)$cluster))
    cluster_m6a <- locus_profile(xc, by_cluster = TRUE)
}
```

Reads from the supplied samples are clustered together within each locus using
Manhattan distances, a nearest-neighbor graph and Leiden communities. At least
three reads and two occupied features are required; otherwise, metadata record
why clustering was skipped. Exact pairwise distances use quadratic memory in
the number of reads.

Clustering adds assignments without changing the assay or its column order.
Settings and clustering provenance are saved in `S4Vectors::metadata(xc)`.
To change a fit, rerun `analyze_loci()`; plot options do not rerun clustering.

### 9. Plot all figures

For the selected locus, draw all three unclustered views:

```r
ComplexHeatmap::draw(plot_locus(x, "heatmap", clustered = FALSE))
print(plot_locus(x, "profile", clustered = FALSE))
print(plot_locus(x, "footprints", clustered = FALSE))
```

Then draw the corresponding clustered views when clustering completed:

```r
if (identical(clustering_status, "complete")) {
    ComplexHeatmap::draw(plot_locus(xc, "heatmap", clustered = TRUE))
    print(plot_locus(xc, "profile", clustered = TRUE))
    print(plot_locus(xc, "footprints", clustered = TRUE))
}
```

The heatmap displays reads vertically and genomic positions horizontally.
The footprint figure stacks the aggregate m6A profile, m6A read tracks and
available FIRE/nucleosome/TF tracks. Sample annotations appear farthest left,
followed by cluster in clustered views. TSS markers appear when TSS metadata
are available. `cluster_colors(colData(xc)$cluster)` provides the same palette
for custom figures after successful clustering.

**TF footprint lengths are configurable.** The default ranges retain 10–29,
40–59 and 60–79 bp for compatibility with earlier analyses. To use other lengths,
change the setting in step 6 and rerun steps 7–10:

```r
# Choose one option:
tf_length_ranges <- list(c(30, 39), c(80, 120)) # inclusive custom ranges
# tf_length_ranges <- list(c(35, 35))          # exactly 35 bp
# tf_length_ranges <- NULL                    # all positive footprint lengths
```

Bounds must be positive whole numbers, with minimum <= maximum. Ranges cannot
overlap; `list(c(10, 30), c(30, 50))` is invalid because both include 30 bp.
Custom ranges get plot labels and colors automatically; `NULL` displays all TF
lengths as one class. FiberHMM score filtering still applies independently.
Length filtering uses each original footprint's width, before plotting clips
it to the selected locus. These settings do not change m6A clustering.

For custom TF BED6 files, columns must be chromosome, zero-based start, end,
original read ID, score and strand. Supply one prefiltered file per sample and
set `tf_format = "bed6_per_sample"`. Nucleosome lengths are separately
configurable with `nuc_min` and `nuc_max`.

If a supplied nucleosome or TF track has no retained calls, the warning
`no footprints called in the selected region` is expected; analysis still
completes. Check the selected region, matching read IDs and length/score filters.
Omitting optional tracks does not trigger this warning.

To inspect another locus, select its object from each results list. Step 10
exports the figures for **all** analyzed loci automatically.

### 10. Save results

Choose new or empty output directories:

```r
output_dir <- "fiberseq_results"

write_locus_results(
    unclustered_results, file.path(output_dir, "without_clustering")
)
write_locus_results(
    clustered_results, file.path(output_dir, "with_clustering")
)
```

Each locus gets its own directory with complete RDS objects, TSV summaries and
PDF figures. The clustered run exports both unclustered views and, when
clustering succeeds, clustered views plus assignments and graph edges.
`locus.rds` preserves the assay, read metadata, per-read footprints and analysis
provenance. See the output layout below.

```r
list.files(output_dir, recursive = TRUE)

# Reload a locus from the clustered run. Use its exported directory name
# (unsafe filename characters in locus names are replaced by underscores).
exported_locus <- list.files(
    file.path(output_dir, "with_clustering"),
    pattern = "^locus[.]rds$", recursive = TRUE, full.names = TRUE
)[1]
restored <- readRDS(exported_locus)
```

## Data model and invariants

`analyze_loci()` returns a named `SimpleList` of `RangedSummarizedExperiment`
objects. Every object's sparse `m6a` assay is **genomic positions × reads**.
`rowRanges` is a single-base `GRanges` whose names/order match the assay rows.
`colData` rownames and `RID` match the assay columns exactly. A read ID is
namespaced as `sample_name::original_read_id`.

`colData$read_range` is a `GRanges` column for full alignment coordinates;
`colData$footprints` is a `GRangesList` column containing each read's footprint
and FIRE ranges. `x[positions, reads]` keeps these read annotations aligned.
All in-memory ranges and alignment coordinates are one-based inclusive;
region TSV start/end use BED conventions. Footprint metadata specify RID,
track and class. Row subsetting does not clip stored per-read intervals.

`metadata(x)` contains original region, samples, QC, parameters and clustering
provenance. Generic subsetting does not recompute this run provenance.
`locus_profile(x)` recomputes fractions from the current experiment.

Only the longest overlapping alignment per read is considered, and it must
span the entire locus. Zero-call reads are retained in both matrices and
pooled denominators. Unclustered ordering is sample-sheet order, read start,
then read ID. Sparse zeros mean no observed call, not a measured unmethylated
adenine at every genomic base. Tutorial profile smoothing uses observed-call
positions; all-zero regions return zero profiles and empty regions NA profiles.

Clustering adds assignments to `colData` without changing the assay or its
order. Unclustered plots and pooled profiles are independent of clustering
parameters. The cluster plot orders by cluster, sample, start and read ID.
Clustering requires at least three reads and two occupied features. Exact
Manhattan distances use quadratic memory in read count. This is Leiden
clustering, not topic modelling.

A supplied nucleosome or TF track with no retained calls emits exactly:

> no footprints called in the selected region

Analysis still completes, with empty per-read ranges and the status recorded.
An omitted optional track does not trigger this warning. A supplied missing
file is an input error.

## Output names

```text
results/ALAS1/
  locus.rds                    # complete RangedSummarizedExperiment
  m6a_matrix.rds               # sparse positions x reads assay
  reads.tsv                    # read metadata; range columns stay in locus.rds
  qc.tsv
  region.tsv                   # BED coordinates and one-based TSS
  unclustered/
    m6a_heatmap.pdf
    aggregate_m6a_profile.pdf
    footprint_tracks.pdf       # when optional tracks were supplied
    accessibility_profile.tsv
    read_order.tsv
    manifest.tsv
  clustered/                   # only after successful clustering
    m6a_heatmap.pdf
    aggregate_m6a_profile.pdf   # profiles separated by cluster
    footprint_tracks.pdf
    accessibility_profile.tsv
    read_order.tsv
    manifest.tsv
    assignments.tsv
    edges.tsv
```

## Further tutorials

Three executable vignettes cover all ten public functions:

- [Function reference and tutorials](vignettes/function-reference.Rmd):
  a package overview, functions grouped by task, and runnable examples for every
  public function.
- [ALAS1 with and without clustering](vignettes/ALAS1-clustered-and-unclustered.Rmd):
  all three plot types in both modes, adjustable clustering and result export.
- [Overview of the functions](vignettes/LCL-input-structures.Rmd):
  actual AL6/AL8 sample sheets and records, coordinates, annotations, optional
  tracks, and the mapping to `GRanges` and `SummarizedExperiment`.

The additional [LCL fiber-seq data example](https://crystalshann.github.io/plot_selected_Fiber-seq_loci/LCL-fiber-seq-data-example.html)
uses all 31 samples from the shared project at ALAS1, with unclustered and
optional clustered panels. Its [R Markdown source](https://crystalshann.github.io/plot_selected_Fiber-seq_loci/LCL-fiber-seq-data-example.Rmd)
uses the package functions with configurable paths and clustering settings.
Reading the HTML requires no access to the shared files; re-running this
particular example requires those files. The bundled-data vignettes above
remain portable.

After installing a built release, open them with:

```r
vignette("function-reference", package = "selectFiberseqLociWorkflow")
vignette("ALAS1-clustered-and-unclustered", package = "selectFiberseqLociWorkflow")
vignette("LCL-input-structures", package = "selectFiberseqLociWorkflow")
```

To include vignettes when installing from GitHub, set `build_vignettes = TRUE`
in `remotes::install_github()`, with `knitr`, `rmarkdown` and Pandoc available.
To knit local HTML previews from the package root after installing the package:

```r
rmarkdown::render("vignettes/function-reference.Rmd")
rmarkdown::render("vignettes/ALAS1-clustered-and-unclustered.Rmd")
rmarkdown::render("vignettes/LCL-input-structures.Rmd")
```

Local HTML previews are written beside their R Markdown sources and ignored
by Git and source-package builds. `R CMD build` builds the installed vignettes
from their R Markdown sources. You can share the self-contained HTML files
directly or attach them to a GitHub Release.
