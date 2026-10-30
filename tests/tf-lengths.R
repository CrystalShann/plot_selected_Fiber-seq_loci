library(selectFiberseqLociWorkflow)
library(SummarizedExperiment)
library(GenomicRanges)

plant <- system.file("extdata", "plant", package = "selectFiberseqLociWorkflow")
samples <- make_sample_sheet(file.path(plant, "samples.tsv"))
regions <- prepare_regions("Chr5:1901-2600")
# Exercise inclusive endpoints, former gaps, exact widths and lengths >79 bp.
widths <- c(1, 9, 10, 29, 30, 35, 39, 40, 59, 60, 79, 80, 120, 121)
starts <- 2000 + seq_along(widths) * 10
scratch <- tempfile("tf-length-tests-")
dir.create(scratch)
write_track <- function(rows, name) {
    bed <- file.path(scratch, paste0(name, ".bed"))
    write.table(rows, bed, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = FALSE)
    gz <- Rsamtools::bgzip(bed, overwrite = TRUE)
    Rsamtools::indexTabix(gz, format = "bed")
    gz
}
bed6 <- data.frame("Chr5", starts, starts + widths, "read1", 0, "+")
bed15 <- data.frame("Chr5", 0, 10000, "read1", 60, "+", 0, 10000,
    "0,0,0", length(widths), paste(widths, collapse = ","),
    paste(starts, collapse = ","), paste(rep(80, length(widths)), collapse = ","),
    paste(rep(0, length(widths)), collapse = ","),
    paste(rep(0, length(widths)), collapse = ","))
paths <- c(bed6_per_sample = write_track(bed6, "custom"),
    bed15_fiberhmm_tf = write_track(bed15, "fiberhmm"))
get_tf <- function(x) {
    fp <- unlist(colData(x)$footprints, use.names = FALSE)
    fp[S4Vectors::mcols(fp)$track == "tf"]
}
run <- function(...) suppressWarnings(analyze_loci(samples, regions,
    clustering = FALSE, ...))[[1]]
for (format in names(paths)) {
    samples$tf <- paths[[format]]
    samples$tf_format <- format
    default <- run()
    stopifnot(identical(sort(width(get_tf(default))),
        as.integer(c(10, 29, 40, 59, 60, 79))))
    custom <- list(c(80, 120), c(30, 39)) # unsorted, disjoint ranges are valid
    x <- run(tf_length_ranges = custom)
    stopifnot(identical(sort(width(get_tf(x))), as.integer(c(30, 35, 39, 80, 120))),
        identical(S4Vectors::metadata(x)$parameters$tf_length_ranges, custom),
        identical(assay(default), assay(x)),
        identical(locus_profile(default), locus_profile(x)))
    exact <- run(tf_length_ranges = list(c(35, 35)))
    stopifnot(identical(width(get_tf(exact)), 35L))
    all <- run(tf_length_ranges = NULL)
    stopifnot(identical(sort(width(get_tf(all))), as.integer(widths)))
    for (object in list(x, exact, all)) {
        plot <- plot_locus(object, "footprints")
        built <- ggplot2::ggplot_build(plot[[3]])
        scale <- built$plot$scales$get_scales("colour")
        colors <- scale$map(unique(mcols(get_tf(object))$class))
        stopifnot(all(!is.na(colors)), all(colors != "grey50"))
    }
    # Score filtering remains independent, and does not apply to BED6.
    scored <- run(tf_length_ranges = NULL, min_tf_score = 90)
    stopifnot(length(get_tf(scored)) ==
        if (format == "bed6_per_sample") length(widths) else 0L)
    empty <- run(tf_length_ranges = list(c(500, 600)))
    stopifnot(length(get_tf(empty)) == 0L)
    invisible(ggplot2::ggplot_build(plot_locus(empty, "footprints")[[3]]))
    output <- file.path(scratch, format)
    write_locus_results(x, output)
    restored <- readRDS(list.files(output, "^locus[.]rds$",
        recursive = TRUE, full.names = TRUE)[1])
    stopifnot(identical(colData(restored)$footprints, colData(x)$footprints),
        identical(S4Vectors::metadata(restored)$parameters$tf_length_ranges, custom))
}
for (bad in list(list(), c(10, 20), list(c(0, 10)), list(c(20, 10)),
    list(c(1.5, 10)), list(c(NA, 10)), list(c(1, Inf)), list("10"),
    list(c(10, 30), c(30, 50)), list(c(1, 100), c(5, 10)))) {
    err <- tryCatch(run(tf_length_ranges = bad), error = identity)
    stopifnot(inherits(err, "error"),
        grepl("tf_length_ranges", conditionMessage(err)))
}
cat("Configurable TF lengths, both formats, plots and exports passed\n")
