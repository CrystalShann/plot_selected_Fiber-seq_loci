# selectFiberseqLociWorkflow
This is to plot m6a heatmap for a selected locus and have the option to do clustering on methylation signal

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
install.packages(c("Matrix", "data.table", "dplyr", "igraph", "withr", "ggplot2"))
install.packages("selectFiberseqLociWorkflow_0.99.0.tar.gz",
    repos = NULL, type = "source")
```


## Quick start

```r
library(selectFiberseqLociWorkflow)
library(SummarizedExperiment)

example <- example_data()  # copies portable inputs into a writable directory
result <- analyze_loci(example$samples, example$regions, clustering = FALSE)
x <- result[["ALAS1"]]
x
rowRanges(x)                 # one genomic base per assay row
colData(x)                  # one read per assay column
colData(x)$read_range        # GRanges of alignments, one-based inclusive
colData(x)$footprints        # GRangesList with one element per read
plot_locus(x, "heatmap")
plot_locus(x, "profile")
plot_locus(x, "footprints")
write_locus_results(result, "my-fiberseq-results")
```

Enable clustering explicitly:

```r
clustered <- analyze_loci(example$samples, example$regions,
    clustering = TRUE, k_neighbors = 10, resolution = 1, seed = 1)
plot_locus(clustered[["ALAS1"]], "heatmap", clustered = TRUE)
plot_locus(clustered[["ALAS1"]], "profile", clustered = TRUE)
write_locus_results(clustered, "my-clustered-results")
```

The full executable tutorial is
[vignettes/ALAS1-clustered-and-unclustered.Rmd](vignettes/ALAS1-clustered-and-unclustered.Rmd).
After installing a built release, open it with:

```r
vignette("ALAS1-clustered-and-unclustered", package = "selectFiberseqLociWorkflow")
```

## Your own samples and regions

```r
samples <- make_sample_sheet(data.frame(
    sample_name = c("plant_A", "plant_B"),
    sample_label = c("Control", "Treatment"),
    m6a = "/data/{sample_name}/m6a.{chr}.bed.gz",
    m6a_format = "bed12_fibertools",
    nucleosome = "/data/{sample_name}/nuc.{chr}.bed.gz",
    nucleosome_format = "bed13_fiberhmm",
    tf = "/data/{sample_name}/tf.{chr}.bed.gz",
    tf_format = "bed15_fiberhmm_tf"), output = "samples.tsv")

regions <- prepare_regions(annotation = "plant.gff3",
    genes = c("AT5G00010", "PlantGene2"), output = "regions.tsv")
# Alternatives:
# regions <- prepare_regions("regions.bed")
# regions <- prepare_regions(c("Chr5:1001-3000", "scaffold_42:4001-6000"))
# regions <- prepare_regions(data.frame(chr="Chr5", start=1000, end=3000),
#     coordinate_system="bed")
# regions <- prepare_regions(my_GRanges)

validate_inputs(samples, regions)
result <- analyze_loci(samples, regions, clustering = FALSE)
```

Sample sheets are `S4Vectors::DataFrame` objects and can be written to TSV.
Relative paths in a TSV are resolved against the TSV directory. In-memory
relative paths use `base_dir`. Only `{sample_name}` and `{chr}` placeholders
are supported. Omit optional track columns or leave their paths blank.
Supplied paths must exist with adjacent `.tbi` indexes.

| Track | Explicit format | Record interpretation |
| --- | --- | --- |
| m6A (required) | `bed12_fibertools` | Fibertools blocks; terminal sentinel blocks removed |
| Nucleosome | `bed12_fibertools` or `bed13_fiberhmm` | Fibertools sentinels removed; FiberHMM blocks retained |
| TF footprint | `bed15_fiberhmm_tf` or `bed6_per_sample` | FiberHMM per-block score threshold; BED6 already filtered |
| FIRE | `bed11_fire` | Per-read FIRE intervals |

Formats are validated explicitly; provenance is never guessed. The m6A
track is fibertools BED12. FiberHMM-only modification output is not a
supported replacement for that required track. No allele-specific RDS or
phasing input is required.

GENCODE and plant GTF/GFF3 annotations are accepted when they contain
transcript/mRNA rows with gene relationships. Canonical transcripts are not
required. Distinct TSSs are retained, and shared TSSs retain their transcript
IDs. Default promoters cover 1,000 bp on each side: one-based
`[TSS - 1000, TSS + 999]`. `extract_canonical_tss()` is optional and reports
missing/ambiguous canonical selections. It defaults to the GENCODE
`tag=Ensembl_canonical`; other attribute/value pairs are supported.

Contig names and gene identifiers are preserved. Contigs must match tracks
exactly, including case, and annotations must use the same genome assembly.
Synthetic plant fixtures test both strands, nonhuman contigs, multiple
transcripts and absent canonical tags. This establishes input compatibility,
not biological validation of every plant caller. Promoters extending before
position 1 are rejected rather than silently changed.

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

## Running example data
```r
remotes::install_github("CrystalShann/plot_selected_Fiber-seq_loci")
library(selectFiberseqLociWorkflow)
e <- example_data()
result <- analyze_loci(e$samples, e$regions, clustering = FALSE)
```

