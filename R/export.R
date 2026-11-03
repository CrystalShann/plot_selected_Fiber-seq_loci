# Export locus experiments, summaries and figures.
write_locus_results <- function(results, directory) {
    if (methods::is(results, "RangedSummarizedExperiment"))
        results <- S4Vectors::SimpleList(stats::setNames(list(results),
            names(S4Vectors::metadata(results)$region)))
    if (!methods::is(results, "List") || !length(results))
        stop("results must be a nonempty SimpleList of locus experiments")
    ids <- names(results)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)))
        stop("Results require locus names")
    safe <- gsub("[^A-Za-z0-9_.-]", "_", ids)
    if (anyDuplicated(safe) || any(safe %in% c(".", "..")))
        stop("Locus names collide as filenames; choose unique safe names")
    if (dir.exists(directory) && length(list.files(directory,
        all.files = TRUE, no.. = TRUE))) stop("Choose a new or empty output directory")
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_along(results)) {
        x <- results[[i]]
        .check_locus(x)
        out <- file.path(directory, safe[i])
        dir.create(out, showWarnings = FALSE)
        saveRDS(x, file.path(out, "locus.rds"))
        saveRDS(SummarizedExperiment::assay(x, "m6a"),
            file.path(out, "m6a_matrix.rds"))
        cd <- SummarizedExperiment::colData(x)
        flat <- setdiff(names(cd), c("read_range", "footprints"))
        .write_table(as.data.frame(cd[, flat, drop = FALSE]),
            file.path(out, "reads.tsv"))
        md <- S4Vectors::metadata(x)
        if (!is.null(md$autocorrelation))
            write_autocorrelation_results(md$autocorrelation,
                file.path(out, "autocorrelation"))
        .write_table(as.data.frame(md$qc), file.path(out, "qc.tsv"))
        .write_table(.regions_table(md$region), file.path(out, "region.tsv"))
        .write_panels(x, file.path(out, "unclustered"), FALSE)
        if ("cluster" %in% names(cd)) {
            .write_panels(x, file.path(out, "clustered"), TRUE)
            .write_table(as.data.frame(cd[, c("RID", "cluster")]),
                file.path(out, "clustered", "assignments.tsv"))
            .write_table(igraph::as_data_frame(md$clustering$graph, what = "edges"),
                file.path(out, "clustered", "edges.tsv"))
        }
    }
    invisible(normalizePath(directory))
}

# Save view-specific plots, profiles and a manifest.
.write_panels <- function(x, directory, clustered) {
    dir.create(directory, showWarnings = FALSE)
    ord <- .read_order(x, clustered)
    .write_table(data.frame(RID = colnames(x)[ord]),
        file.path(directory, "read_order.tsv"))
    .write_table(as.data.frame(locus_profile(x, by_cluster = clustered)),
        file.path(directory, "accessibility_profile.tsv"))
    types <- c("heatmap", "profile", "footprints")
    files <- c("m6a_heatmap.pdf", "aggregate_m6a_profile.pdf", "footprint_tracks.pdf")
    md <- S4Vectors::metadata(x)
    status <- rep("written", length(types))
    detail <- rep("", length(types))
    if (!length(md$supplied_tracks)) {
        status[3] <- "omitted"
        detail[3] <- "No optional tracks supplied"
    }
    if (!ncol(x)) {
        status[] <- "empty"
        detail[] <- "No full-span reads"
    }
    if (length(md$empty_footprint_tracks))
        detail[3] <- "no footprints called in the selected region"
    for (i in which(status == "written"))
        .save_plot(plot_locus(x, types[i], clustered),
            file.path(directory, files[i]), c(10, 3.5, 14)[i])
    .write_table(data.frame(panel = types, status = status, detail = detail),
        file.path(directory, "manifest.tsv"))
}

# Write a heatmap or ggplot figure to PDF.
.save_plot <- function(plot, path, height) {
    grDevices::pdf(path, width = 10, height = height)
    on.exit(grDevices::dev.off())
    if (methods::is(plot, "Heatmap")) ComplexHeatmap::draw(plot) else print(plot)
    invisible(path)
}
