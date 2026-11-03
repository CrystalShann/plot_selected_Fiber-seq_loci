\name{locus_autocorrelation}
\alias{locus_autocorrelation}
\alias{autocorrelation_profile}
\alias{plot_autocorrelation}
\alias{write_autocorrelation_results}
\title{Independent native R autocorrelation of selected-locus reads}
\description{Compute mean-centered per-read autocorrelations, optional ACF
clusters and repeat peaks; summarize, plot and export the results.}
\usage{
locus_autocorrelation(x, max_lag = NULL, clustering = FALSE,
    n_pcs = 50L, n_neighbors = 10L, resolution = 0.4, seed = 0L)
autocorrelation_profile(x,
    group_by = c("all", "sample", "cluster", "acf_cluster"))
plot_autocorrelation(x, type = c("profile", "heatmap"),
    group_by = c("all", "sample", "cluster", "acf_cluster"), lag_range = NULL)
write_autocorrelation_results(x, directory)
}
\arguments{
\item{x}{For locus_autocorrelation, a locus RangedSummarizedExperiment with
binary m6a calls at consecutive genomic bases. For the other functions, its
returned ACF SummarizedExperiment.}
\item{max_lag}{Largest lag in bp, including lag zero. NULL uses locus width
minus one. Must be less than locus width; at least two bases are required.}
\item{clustering}{Whether to cluster ACF profiles, default FALSE. Independent
of the Manhattan clustering option in analyze_loci. Requires uwot and irlba.}
\item{n_pcs}{Requested centered, unscaled PCA components, integer at least 2.}
\item{n_neighbors}{Requested correlation-distance neighbors including self,
integer at least 2. Capped at usable reads minus one.}
\item{resolution}{Positive resolution for Leiden on the ACF fuzzy graph.}
\item{seed}{Nonnegative integer seed. The caller's RNG state is restored.}
\item{group_by}{Pool reads, group by sample, existing Manhattan cluster, or
independent ACF cluster. Unassigned ACF reads are excluded only from ACF-cluster
views and remain in all/sample views and read exports.}
\item{type}{Mean profile or per-read heatmap.}
\item{lag_range}{Inclusive integer bounds for displayed lags. NULL omits zero
when positive lags are available.}
\item{directory}{New or empty output directory.}
}
\value{locus_autocorrelation returns a SummarizedExperiment with an acf assay
(lags by reads), rowData$lag_bp, aligned read metadata, validity flags, repeat
peaks and optional acf_cluster, acf_cluster_status and acf_umap1/2 columns.
Original Manhattan cluster labels are preserved. autocorrelation_profile
returns a DataFrame of means, SD, SEM and valid/invalid counts.
plot_autocorrelation returns a ggplot or ComplexHeatmap object.
write_autocorrelation_results invisibly returns the output path.}
\details{Uses stats::acf with a fixed lag-zero variance denominator and no
smoothing. Constant reads have NA profiles. SD and SEM require two valid reads.
Repeat peaks select the strongest positive local maximum at 140--250 bp;
max_lag must be at least 251 for the entire band, otherwise peaks are NA.
ACF clustering uses PCA, correlation neighbors, a uwot fuzzy graph and igraph
Leiden; UMAP is for display. Labels are zero-based and ordered by cluster size.
Fewer than four valid reads, identical profiles and undefined score correlations
have explicit statuses. Computation uses all consecutive bases, independently
of Manhattan binning. No Python, phasing files or external project paths are
required. Allele-specific LCL reporting is not part of this general locus API.
Metadata and graph fits describe the original run and are not recalculated on
subsetting. Call locus_autocorrelation on a subset to recompute its analysis.
Exports contain autocorrelation.rds, lags.tsv, reads.tsv, group summaries and
profile/heatmap PDFs, plus ACF edges when clustering was requested.
write_locus_results automatically writes these in an autocorrelation subdirectory.}
\examples{
e <- example_data()
x <- analyze_loci(e$samples, e$regions[1], clustering = FALSE,
    autocorrelation = TRUE)[[1]]
a <- S4Vectors::metadata(x)$autocorrelation
autocorrelation_profile(a, "sample")
plot_autocorrelation(a, group_by = "sample")
}
