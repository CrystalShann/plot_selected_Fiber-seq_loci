\name{plot_locus}
\alias{plot_locus}
\alias{locus_profile}
\alias{write_locus_results}
\title{Plot, summarize and export locus experiments}
\description{Generate per-read heatmaps, pooled accessibility profiles and
available footprint tracks directly from RangedSummarizedExperiment objects.}
\usage{
plot_locus(x, type = c("heatmap", "profile", "footprints"), clustered = FALSE)
locus_profile(x, by_cluster = FALSE)
write_locus_results(results, directory)
}
\arguments{
\item{x}{A locus RangedSummarizedExperiment from \code{\link{analyze_loci}}.}
\item{type}{Plot to construct.}
\item{clustered, by_cluster}{Use existing cluster assignments for ordering or
profile groups; does not run clustering. FALSE uses sample, read start, RID order.}
\item{results}{A locus experiment or named SimpleList of locus experiments.}
\item{directory}{New or empty output directory.}
}
\value{\code{plot_locus} returns a ComplexHeatmap Heatmap (heatmap), a ggplot
(profile), or a patchwork composition of ggplots (footprints).
\code{locus_profile} returns a DataFrame with one-based pos, base, group, cov,
met, frac, smooth_cov and smooth_frac. \code{write_locus_results} returns the
output directory invisibly.}
\details{The pooled denominator includes every retained read, including reads
with zero calls. Profiles use observed m6A positions and the tutorial sliding
window; all-zero loci use all positions. Empty loci have NA fractions.
Plots use the current assay columns and do not alter the experiment.
Heatmaps use black m6A calls, genomic coordinate ticks and a dashed TSS line.
The footprint composition stacks the accessibility profile, blue m6A calls,
and overlaid FIRE (orange), nucleosome (dark blue), and TF footprint tracks
(purple, magenta and cyan for the default size classes, with automatic colors
for custom TF length ranges).
Every read panel places sample annotations farthest left, followed by cluster
annotations when enabled. Clusters use \code{\link{cluster_colors}} throughout.
Sample colors use the Dynamic palette in sample-sheet label order.
Clustered read panels split rows by cluster; within each cluster, sample,
read start and read ID determine order. Coordinates increase left to right on
both strands, as in the original extraction tutorial. On the minus strand,
upstream of the TSS is to the right.

Each locus output directory contains locus.rds (the complete experiment),
m6a_matrix.rds, reads.tsv, qc.tsv and region.tsv. The unclustered directory
contains m6a_heatmap.pdf, aggregate_m6a_profile.pdf, footprint_tracks.pdf when
available, accessibility_profile.tsv, read_order.tsv and manifest.tsv.
Successful clustering adds a clustered directory with corresponding plots,
assignments.tsv and edges.tsv. Names unsafe for file paths are replaced by
underscores, with collisions rejected. Empty panels are recorded in manifest.tsv.
Region TSV coordinates are BED; all in-memory ranges and plot coordinates are
one-based inclusive.}
\examples{
e <- example_data()
x <- analyze_loci(e$samples, e$regions[1], clustering = FALSE)[[1]]
head(locus_profile(x))
plot_locus(x, "profile")
plot_locus(x, "footprints")
write_locus_results(x, tempfile("locus-results-"))
}
