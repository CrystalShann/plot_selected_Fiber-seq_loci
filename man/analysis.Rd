\name{analyze_loci}
\alias{analyze_loci}
\title{Analyze loci using native Bioconductor containers}
\description{Build sparse modification assays from full-span reads and
optionally cluster them using Manhattan distances and Leiden communities.}
\usage{
analyze_loci(samples, regions, clustering, window_size = 0,
    k_neighbors = 10, resolution = 1, seed = 1, nuc_min = 130,
    nuc_max = 160, min_tf_score = 50, smoothing_bp = 10,
    tf_length_ranges = list(c(10, 29), c(40, 59), c(60, 79)),
    autocorrelation = FALSE, autocorrelation_options = list())
}
\arguments{
\item{autocorrelation}{Logical, default FALSE. Compute independent per-read
ACFs using the same sample sheet and full-span reads. Works with either value
of clustering.}
\item{autocorrelation_options}{Named list of arguments to
\code{\link{locus_autocorrelation}}, excluding x. For example,
\code{list(max_lag = 500, clustering = TRUE, resolution = 0.4)}.
These settings affect only ACF analysis.}
\item{samples}{A sample DataFrame, data.frame or TSV path; see
\code{\link{make_sample_sheet}}.}
\item{regions}{A named GRanges, normally from \code{\link{prepare_regions}}.}
\item{clustering}{Required explicit TRUE or FALSE. Saved in result metadata.}
\item{window_size}{Nonnegative integer clustering bin width in bases. The
default, 0, uses individual observed sites. Positive widths average each
read's calls across observed sites within each occupied bin. Bins start at
the locus's first base. The stored assay remains single-base.}
\item{k_neighbors}{Positive integer number of Manhattan nearest neighbors,
default 10. The effective value is capped at the number of reads minus one.}
\item{resolution}{Nonnegative Leiden resolution, default 1. Larger values
generally favor finer partitions.}
\item{seed}{Nonnegative integer seed for Leiden clustering.}
\item{nuc_min, nuc_max}{Inclusive nucleosome size limits in bases.}
\item{min_tf_score}{Minimum FiberHMM TF per-block score. BED6 is prefiltered.}
\item{tf_length_ranges}{NULL to retain all positive TF footprint lengths, or a
nonempty list of nonoverlapping inclusive c(min, max) pairs in bp. Bounds must
be positive integers. For example, list(c(30, 39), c(80, 120)) selects two
custom classes and list(c(35, 35)) selects exactly 35 bp. Applies to both
FiberHMM BED15 and prefiltered BED6 TF inputs, independently of score filtering.
Defaults preserve the original 10--29, 40--59 and 60--79 bp classes.}
\item{smoothing_bp}{Positive integer accessibility smoothing width in bases.}
}
\value{A named S4Vectors SimpleList of RangedSummarizedExperiment objects.
Each object has a sparse \code{m6a} assay with genomic bases in rows and reads
in columns. Its single-base GRanges rowRanges match assay rows; its colData
rownames and RID match assay columns. colData includes original_RID,
sample_name, sample_label, one-based alignment start/end, a GRanges
\code{read_range} column, and a GRangesList \code{footprints} column.
Footprint metadata include RID, track and class. Successful clustering adds
colData$cluster without changing assay or read order. Enabled autocorrelation
adds an independent lag-by-read SummarizedExperiment in
\code{metadata(x)$autocorrelation}, with a separate status in
\code{metadata(x)$autocorrelation_status}. ACF errors produce a warning and
status while preserving the locus and its Manhattan clustering results.}
\details{Only the longest overlapping alignment per read is considered; it
must span the entire locus. Reads with zero m6A calls remain in the assay and
denominators. Fibertools sentinel blocks are removed; FiberHMM blocks are not.
TF length classes are set by \code{tf_length_ranges}; filtering uses original
footprint widths before clipping for display. A supplied footprint track
with no retained calls triggers the warning ``no footprints called in the
selected region''. Empty loci return zero-column experiments.

Normal SummarizedExperiment subsetting keeps rowRanges, assay, colData,
alignment ranges and per-read footprints synchronized. Metadata contain the
original region, sample sheet, QC, parameters, optional autocorrelation and clustering fit as
run provenance; these summaries are not recalculated by generic subsetting.
\code{locus_profile} recalculates denominators for the current columns.
The clustering controls \code{window_size}, \code{k_neighbors} and
\code{resolution} are ignored when \code{clustering = FALSE}. They do not
change the unclustered assay or pooled profile. \code{smoothing_bp} controls
profile smoothing independently of clustering bins. Settings are saved in
\code{metadata(x)$parameters}; run \code{analyze_loci} again to change a fit.
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
