\name{analyze_loci}
\alias{analyze_loci}
\title{Analyze loci using native Bioconductor containers}
\description{Build sparse modification assays from full-span reads and
optionally cluster them using Manhattan distances and Leiden communities.}
\usage{
analyze_loci(samples, regions, clustering, window_size = 0,
    k_neighbors = 10, resolution = 1, seed = 1, nuc_min = 130,
    nuc_max = 160, min_tf_score = 50, smoothing_bp = 10)
}
\arguments{
\item{samples}{A sample DataFrame, data.frame or TSV path; see
\code{\link{make_sample_sheet}}.}
\item{regions}{A named GRanges, normally from \code{\link{prepare_regions}}.}
\item{clustering}{Required explicit TRUE or FALSE. Saved in result metadata.}
\item{window_size}{Clustering feature width in bases; zero uses observed sites.}
\item{k_neighbors}{Positive integer number of Manhattan nearest neighbors.}
\item{resolution}{Nonnegative Leiden resolution.}
\item{seed}{Nonnegative integer seed for Leiden clustering.}
\item{nuc_min, nuc_max}{Inclusive nucleosome size limits in bases.}
\item{min_tf_score}{Minimum FiberHMM TF per-block score. BED6 is prefiltered.}
\item{smoothing_bp}{Positive integer accessibility smoothing width in bases.}
}
\value{A named S4Vectors SimpleList of RangedSummarizedExperiment objects.
Each object has a sparse \code{m6a} assay with genomic bases in rows and reads
in columns. Its single-base GRanges rowRanges match assay rows; its colData
rownames and RID match assay columns. colData includes original_RID,
sample_name, sample_label, one-based alignment start/end, a GRanges
\code{read_range} column, and a GRangesList \code{footprints} column.
Footprint metadata include RID, track and class. Successful clustering adds
colData$cluster without changing assay or read order.}
\details{Only the longest overlapping alignment per read is considered; it
must span the entire locus. Reads with zero m6A calls remain in the assay and
denominators. Fibertools sentinel blocks are removed; FiberHMM blocks are not.
TF size classes are [10,30), [40,60), [60,80) bp. A supplied footprint track
with no retained calls triggers the warning ``no footprints called in the
selected region''. Empty loci return zero-column experiments.

Normal SummarizedExperiment subsetting keeps rowRanges, assay, colData,
alignment ranges and per-read footprints synchronized. Metadata contain the
original region, sample sheet, QC, parameters and optional clustering fit as
run provenance; these summaries are not recalculated by generic subsetting.
\code{locus_profile} recalculates denominators for the current columns.
Clustering requires three reads and two occupied features; otherwise metadata
record why it was skipped. Exact Manhattan distances use quadratic memory in
the number of reads. This package analyzes processed tracks, not raw reads.}
\examples{
e <- example_data()
x <- analyze_loci(e$samples, e$regions[1], clustering = FALSE)[[1]]
SummarizedExperiment::rowRanges(x)
SummarizedExperiment::colData(x)[, c("RID", "sample_name")]
small <- x[1:20, 1:3]
stopifnot(identical(colnames(small),
    rownames(SummarizedExperiment::colData(small))))
}
