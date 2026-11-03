library(selectFiberseqLociWorkflow)
library(SummarizedExperiment)
e <- example_data()
same_fit <- function(x, y) {
    a <- metadata(x)$clustering; b <- metadata(y)$clustering
    stopifnot(identical(igraph::as_data_frame(a$graph), igraph::as_data_frame(b$graph)))
    a$graph <- b$graph <- NULL
    identical(a, b)
}
# Exercise all combinations using the identical sample sheet.
none <- analyze_loci(e$samples, e$regions[1], clustering = FALSE)[[1]]
acf_only <- analyze_loci(e$samples, e$regions[1], clustering = FALSE,
    autocorrelation = TRUE)[[1]]
clustered <- analyze_loci(e$samples, e$regions[1], clustering = TRUE)[[1]]
both <- analyze_loci(e$samples, e$regions[1], clustering = TRUE,
    autocorrelation = TRUE)[[1]]
stopifnot(identical(assay(none), assay(acf_only)),
    identical(colData(none), colData(acf_only)),
    identical(assay(clustered), assay(both)),
    identical(colData(clustered), colData(both)),
    same_fit(clustered, both),
    identical(locus_profile(clustered), locus_profile(both)),
    identical(metadata(both)$autocorrelation_status, "complete"),
    is.null(metadata(clustered)$autocorrelation))
a <- metadata(both)$autocorrelation
stopifnot(identical(colnames(a), colnames(both)),
    identical(colData(a)$cluster, colData(both)$cluster))
# Verify numerical normalization against its explicit definition.
v <- as.numeric(assay(none)[, which(colData(a)$acf_valid)[1]])
z <- v - mean(v)
expected <- vapply(0:5, function(k)
    sum(z[seq_len(length(z)-k)] * z[seq.int(1+k,length(z))])/sum(z^2), numeric(1))
stopifnot(isTRUE(all.equal(unname(assay(a)[1:6, which(colData(a)$acf_valid)[1]]), expected)))
# Constant and empty inputs retain valid container shapes and explicit statuses.
constant <- none[, 1:2]
assay(constant) <- matrix(0, nrow(constant), ncol(constant), dimnames=dimnames(constant))
ac <- locus_autocorrelation(constant)
stopifnot(all(is.na(assay(ac))), !any(colData(ac)$acf_valid),
    all(is.na(autocorrelation_profile(ac)$mean_acf)))
empty <- locus_autocorrelation(none[, integer()])
stopifnot(ncol(empty) == 0L, nrow(empty) == nrow(none))
stopifnot(inherits(try(locus_autocorrelation(none, max_lag=nrow(none)), silent=TRUE), "try-error"))
failed <- suppressWarnings(analyze_loci(e$samples, e$regions[1], clustering=TRUE,
    autocorrelation=TRUE, autocorrelation_options=list(max_lag=1000000))[[1]])
stopifnot(identical(colData(failed), colData(clustered)),
    same_fit(failed, clustered),
    startsWith(metadata(failed)$autocorrelation_status, "failed:"))
if (requireNamespace("uwot", quietly=TRUE) && requireNamespace("irlba", quietly=TRUE)) {
    set.seed(14)
    rng <- .Random.seed
    fitted <- locus_autocorrelation(clustered, clustering=TRUE)
    stopifnot(identical(rng, .Random.seed),
        identical(colData(fitted)$cluster, colData(clustered)$cluster),
        "acf_cluster" %in% names(colData(fitted)))
    repeat_fit <- locus_autocorrelation(clustered, clustering=TRUE)
    stopifnot(identical(colData(fitted)$acf_cluster, colData(repeat_fit)$acf_cluster))
    autocorrelation_profile(fitted, "acf_cluster")
    constant_fit <- locus_autocorrelation(constant, clustering=TRUE)
    stopifnot(all(colData(constant_fit)$acf_cluster_status == "zero_variance"))
    plot_autocorrelation(constant_fit, group_by="acf_cluster")
}
# Export separation and round-trip read alignment.
out <- tempfile("acf-export-")
write_locus_results(S4Vectors::SimpleList(test=both), out)
restored <- readRDS(file.path(out,"test","autocorrelation","autocorrelation.rds"))
stopifnot(identical(assay(restored), assay(a)),
    identical(colData(restored), colData(a)), identical(rowData(restored), rowData(a)),
    same_fit(restored, a), file.exists(file.path(out,"test","clustered","assignments.tsv")))
unlink(out, recursive=TRUE)
unlink(e$directory, recursive=TRUE)
