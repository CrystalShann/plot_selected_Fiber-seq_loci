# Independent native R ACF analysis of the same selected-locus reads.
.acf_check_locus <- function(x) {
    if (!methods::is(x, "RangedSummarizedExperiment") ||
        !"m6a" %in% SummarizedExperiment::assayNames(x))
        stop("x must be a RangedSummarizedExperiment with an m6a assay")
    methods::validObject(x)
    rr <- SummarizedExperiment::rowRanges(x)
    cd <- SummarizedExperiment::colData(x)
    if (!methods::is(rr, "GRanges") || any(BiocGenerics::width(rr) != 1L) ||
        length(unique(as.character(GenomicRanges::seqnames(rr)))) > 1L ||
        anyDuplicated(BiocGenerics::start(rr))) stop("Invalid single-locus ranges")
    if (ncol(x) && !identical(colnames(x), as.character(cd$RID)))
        stop("Read IDs and colData must match assay columns")
    invisible(TRUE)
}


.acf_read_order <- function(x, clustered) {
    cd <- SummarizedExperiment::colData(x)
    samples <- S4Vectors::metadata(x)$samples$sample_name
    if (clustered && !"cluster" %in% names(cd))
        stop("No clustering in this result")
    if (clustered) order(cd$cluster, match(cd$sample_name, samples),
        cd$start, cd$RID) else
        order(match(cd$sample_name, samples), cd$start, cd$RID)
}


.acf_sample_colors <- function(x) {
    labels <- unique(S4Vectors::metadata(x)$samples$sample_label)
    stats::setNames(grDevices::hcl.colors(length(labels), "Dynamic"), labels)
}
.acf_save_plot <- function(plot, path, height) {
    grDevices::pdf(path, width = 10, height = height)
    on.exit(grDevices::dev.off())
    if (methods::is(plot, "Heatmap")) ComplexHeatmap::draw(plot) else print(plot)
    invisible(path)
}
.acf_write_table <- function(x, path) {
    utils::write.table(x, file = path, sep = "\t", quote = TRUE,
        row.names = FALSE, na = "NA")
}
locus_autocorrelation <- function(x, max_lag = NULL, clustering = FALSE,
    n_pcs = 50L, n_neighbors = 10L, resolution = 0.4, seed = 0L) {
    .acf_check_locus(x)
    .acf_validate_options(list(max_lag = max_lag, clustering = clustering,
        n_pcs = n_pcs, n_neighbors = n_neighbors, resolution = resolution, seed = seed))
    ranges <- SummarizedExperiment::rowRanges(x)
    positions <- BiocGenerics::start(ranges)
    if (length(positions) < 2L || any(diff(positions) != 1L))
        stop("ACF requires at least two consecutive, increasing genomic positions")
    n <- nrow(x)
    if (is.null(max_lag)) max_lag <- n - 1L
    if (!is.numeric(max_lag) || length(max_lag) != 1L ||
        !is.finite(max_lag) || max_lag != floor(max_lag) ||
        max_lag < 0 || max_lag >= n)
        stop("max_lag must be an integer from 0 to locus width minus one")
    max_lag <- as.integer(max_lag)
    m <- SummarizedExperiment::assay(x, "m6a")
    if (any(!is.finite(m)) || any(m != 0 & m != 1))
        stop("ACF requires finite binary 0/1 calls")
    computed <- .acf_compute(t(as.matrix(m)), n_features = max_lag + 1L)
    valid <- computed$valid
    lags <- seq.int(0L, max_lag)
    values <- t(computed$profiles)
    dimnames(values) <- list(paste0("lag", lags), colnames(x))
    cd <- SummarizedExperiment::colData(x)
    cd$acf_valid <- unname(valid)
    cd$acf_status <- ifelse(valid, "valid", "constant")
    peaks <- vapply(seq_len(ncol(x)), function(j)
        .acf_repeat_peak(values[, j]), numeric(2L))
    cd$positive_local_peak_140_250_bp <- peaks[1L, ]
    cd$peak_acf <- peaks[2L, ]
    cd$peak_band_complete <- rep(max_lag >= 251L, ncol(x))
    md <- S4Vectors::metadata(x)
    md$autocorrelation <- NULL
    md$autocorrelation <- list(max_lag = max_lag, source_width = n,
        source_ranges = ranges,
        normalization = "sum(centered[t] * centered[t+k]) / sum(centered^2)",
        method = "stats::acf; mean-centered linear correlation; no smoothing",
        implementation = "native R autocorrelation",
        source_package_version = as.character(utils::packageVersion(
            "selectFiberseqLociWorkflow")))
    if (clustering) {
        if (max_lag < 2L) stop("ACF clustering requires at least three lag features")
        for (pkg in c("uwot", "irlba"))
            if (!requireNamespace(pkg, quietly = TRUE))
                stop("Install optional package '", pkg, "' for ACF clustering")
        fit <- withr::with_seed(seed, .acf_cluster_profiles(computed$profiles,
            valid, n_pcs, n_neighbors, resolution, seed))
        cd$acf_cluster <- fit$labels
        cd$acf_cluster_status <- fit$status
        cd$acf_umap1 <- fit$embedding[, 1L]
        cd$acf_umap2 <- fit$embedding[, 2L]
        fit$edges$source <- colnames(x)[fit$edges$source]
        fit$edges$target <- colnames(x)[fit$edges$target]
        md$autocorrelation$clustering <- fit
    }
    out <- SummarizedExperiment::SummarizedExperiment(assays = list(acf = values),
        rowData = S4Vectors::DataFrame(lag_bp = lags, row.names = rownames(values)),
        colData = cd, metadata = md)
    methods::validObject(out)
    out
}

.check_acf <- function(x) {
    if (!methods::is(x, "SummarizedExperiment") ||
        methods::is(x, "RangedSummarizedExperiment") ||
        !"acf" %in% SummarizedExperiment::assayNames(x))
        stop("x must be an ACF SummarizedExperiment from locus_autocorrelation()")
    methods::validObject(x)
    lags <- SummarizedExperiment::rowData(x)$lag_bp
    cd <- SummarizedExperiment::colData(x)
    if (!is.numeric(lags) || length(lags) != nrow(x) ||
        any(!is.finite(lags)) || any(lags < 0 | lags != floor(lags)) ||
        anyDuplicated(lags)) stop("Invalid lag_bp row metadata")
    if (!all(c("RID", "acf_valid") %in% names(cd)) ||
        !is.logical(cd$acf_valid) || anyNA(cd$acf_valid) ||
        (ncol(x) > 0L && !identical(colnames(x), as.character(cd$RID))))
        stop("Read IDs and ACF validity flags must align with assay columns")
    m <- SummarizedExperiment::assay(x, "acf")
    if (any(!is.finite(m[, cd$acf_valid, drop = FALSE])) ||
        any(!is.na(m[, !cd$acf_valid, drop = FALSE])))
        stop("Valid ACF columns must be finite; invalid columns must be NA")
    invisible(TRUE)
}

.acf_groups <- function(x, group_by) {
    cd <- SummarizedExperiment::colData(x)
    if (group_by == "all") return(list(All = seq_len(ncol(x))))
    field <- if (group_by == "sample") "sample_name" else "cluster"
    if (!field %in% names(cd)) stop("No ", field, " assignments in this result")
    values <- as.character(cd[[field]])
    if (anyNA(values) || any(!nzchar(values)))
        stop("Missing ", field, " assignments in this result")
    if (group_by == "sample") {
        levels <- S4Vectors::metadata(x)$samples$sample_name
        levels <- unique(c(levels[levels %in% values], values))
    } else {
        levels <- if (is.factor(cd$cluster)) levels(cd$cluster) else unique(values)
        levels <- levels[levels %in% values]
    }
    stats::setNames(lapply(levels, function(g) which(values == g)), levels)
}

autocorrelation_profile <- function(x, group_by = c("all", "sample", "cluster", "acf_cluster")) {
    .check_acf(x)
    group_by <- match.arg(group_by)
    if (group_by == "acf_cluster") {
        x <- .acf_cluster_view(x)
        group_by <- "cluster"
    }
    groups <- .acf_groups(x, group_by)
    empty <- S4Vectors::DataFrame(lag_bp = integer(), group = character(),
        mean_acf = numeric(), sd_acf = numeric(), sem_acf = numeric(),
        n_reads = integer(), n_valid = integer(), n_invalid = integer())
    if (!length(groups) || !nrow(x)) return(empty)
    m <- SummarizedExperiment::assay(x, "acf")
    valid <- SummarizedExperiment::colData(x)$acf_valid
    result <- lapply(names(groups), function(g) {
        all <- groups[[g]]
        selected <- all[valid[all]]
        n <- length(selected)
        z <- m[, selected, drop = FALSE]
        avg <- if (n) rowMeans(z) else rep(NA_real_, nrow(x))
        sd <- if (n > 1L) apply(z, 1L, stats::sd) else rep(NA_real_, nrow(x))
        S4Vectors::DataFrame(lag_bp = SummarizedExperiment::rowData(x)$lag_bp,
            group = rep(g, nrow(x)), mean_acf = unname(avg), sd_acf = unname(sd),
            sem_acf = unname(sd / sqrt(n)), n_reads = rep(length(all), nrow(x)),
            n_valid = rep(n, nrow(x)), n_invalid = rep(length(all) - n, nrow(x)))
    })
    do.call(rbind, result)
}

.acf_lags <- function(x, lag_range) {
    lags <- SummarizedExperiment::rowData(x)$lag_bp
    if (is.null(lag_range)) {
        selected <- if (any(lags > 0)) which(lags > 0) else seq_along(lags)
    } else {
        if (!is.numeric(lag_range) || length(lag_range) != 2L ||
            any(!is.finite(lag_range)) || any(lag_range < 0) ||
            any(lag_range != floor(lag_range)) || lag_range[1] > lag_range[2])
            stop("lag_range must be two nonnegative integers in increasing order")
        selected <- which(lags >= lag_range[1] & lags <= lag_range[2])
        if (!length(selected)) stop("lag_range contains no stored lags")
    }
    selected[order(lags[selected])]
}

.acf_group_colors <- function(x, group_by, groups) {
    if (group_by == "all") return(c(All = "blue"))
    if (group_by == "cluster")
        return(cluster_colors(SummarizedExperiment::colData(x)$cluster)[groups])
    samples <- S4Vectors::metadata(x)$samples
    colors <- .acf_sample_colors(x)
    stats::setNames(unname(colors[samples$sample_label[
        match(groups, samples$sample_name)]]), groups)
}

plot_autocorrelation <- function(x, type = c("profile", "heatmap"),
    group_by = c("all", "sample", "cluster", "acf_cluster"), lag_range = NULL) {
    .check_acf(x)
    type <- match.arg(type)
    group_by <- match.arg(group_by)
    if (group_by == "acf_cluster") {
        x <- .acf_cluster_view(x)
        group_by <- "cluster"
    }
    groups <- .acf_groups(x, group_by)
    selected <- .acf_lags(x, lag_range)
    if (!ncol(x) || !length(selected))
        return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(
            title = "No reads or lags available for autocorrelation"))
    if (type == "heatmap") return(.acf_heatmap(x[selected, ], group_by))
    d <- as.data.frame(autocorrelation_profile(x[selected, ], group_by))
    if (!any(d$n_valid > 0L))
        return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(
            title = "No valid ACFs: all retained reads are constant"))
    d$group <- factor(d$group, levels = names(groups))
    valid <- d[d$n_valid > 0L, , drop = FALSE]
    p <- ggplot2::ggplot(valid, ggplot2::aes(x = lag_bp, y = mean_acf,
        colour = group, group = group)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3)
    if (length(selected) > 1L) p <- p + ggplot2::geom_line(linewidth = 0.6) else
        p <- p + ggplot2::geom_point(size = 1.5)
    labels <- names(groups)
    if (group_by == "sample") {
        samples <- S4Vectors::metadata(x)$samples
        short <- samples$sample_label[match(labels, samples$sample_name)]
        if (!anyDuplicated(short)) labels <- short
    }
    p + ggplot2::scale_colour_manual(values = .acf_group_colors(x, group_by,
            names(groups)), drop = FALSE,
            breaks = names(groups), labels = labels,
            name = if (group_by == "cluster") "Cluster" else "Sample",
            guide = if (group_by == "all") "none" else
                ggplot2::guide_legend(ncol = if (group_by == "sample") 2 else 1)) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::labs(x = "Lag (bp)", y = "Mean autocorrelation",
            title = "Mean per-read autocorrelation",
            subtitle = sprintf("%d valid / %d retained reads; constant reads excluded from means; no smoothing",
                sum(SummarizedExperiment::colData(x)$acf_valid), ncol(x)))
}

.acf_heatmap <- function(x, group_by) {
    clustered <- group_by == "cluster"
    ord <- .acf_read_order(x, clustered)
    cd <- SummarizedExperiment::colData(x)[ord, , drop = FALSE]
    annotation <- data.frame(sample = factor(cd$sample_label,
        levels = names(.acf_sample_colors(x))))
    colors <- list(sample = .acf_sample_colors(x))
    if (clustered) {
        annotation$cluster <- cd$cluster
        colors$cluster <- cluster_colors(cd$cluster)
    }
    m <- t(SummarizedExperiment::assay(x, "acf")[, ord, drop = FALSE])
    lags <- SummarizedExperiment::rowData(x)$lag_bp
    ticks <- unique(round(seq(1L, length(lags), length.out = 5L)))
    ComplexHeatmap::Heatmap(m, name = "ACF",
        col = circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
        na_col = "grey85", cluster_rows = FALSE, cluster_columns = FALSE,
        cluster_row_slices = FALSE, row_split = if (clustered) cd$cluster else NULL,
        row_gap = grid::unit(0.6, "mm"), row_title_rot = 0,
        row_title_gp = grid::gpar(fontsize = 8),
        show_row_names = FALSE, show_column_names = FALSE,
        use_raster = TRUE, raster_quality = 2, raster_resize_mat = FALSE,
        raster_by_magick = FALSE,
        left_annotation = ComplexHeatmap::rowAnnotation(df = annotation,
            col = colors, annotation_legend_param = list(sample = list(ncol = 2))),
        bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
            lag_bp = ComplexHeatmap::anno_mark(at = ticks,
                labels = as.character(lags[ticks]), which = "column", side = "bottom",
                labels_gp = grid::gpar(fontsize = 8))),
        column_title = "Per-read autocorrelation; grey = undefined (constant read)",
        column_title_gp = grid::gpar(fontsize = 11))
}



write_autocorrelation_results <- function(x, directory) {
    .check_acf(x)
    if (length(directory) != 1L || is.na(directory) || !nzchar(directory))
        stop("directory must be a nonempty path")
    if (file.exists(directory) && !dir.exists(directory))
        stop("directory refers to an existing file")
    if (dir.exists(directory) && length(list.files(directory,
        all.files = TRUE, no.. = TRUE))) stop("Choose a new or empty output directory")
    if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE))
        stop("Could not create output directory")
    saveRDS(x, file.path(directory, "autocorrelation.rds"))
    .acf_write_table(as.data.frame(SummarizedExperiment::rowData(x)),
        file.path(directory, "lags.tsv"))
    cd <- SummarizedExperiment::colData(x)
    flat <- names(cd)[vapply(cd, function(v) is.atomic(v) || is.factor(v), logical(1))]
    .acf_write_table(as.data.frame(cd[, flat, drop = FALSE]),
        file.path(directory, "reads.tsv"))
    groups <- c("all", "sample", if ("cluster" %in% names(cd)) "cluster",
        if ("acf_cluster" %in% names(cd)) "acf_cluster")
    fit <- S4Vectors::metadata(x)$autocorrelation$clustering
    if (!is.null(fit)) .acf_write_table(fit$edges, file.path(directory, "acf_edges.tsv"))
    for (group in groups) {
        .acf_write_table(as.data.frame(autocorrelation_profile(x, group)),
            file.path(directory, paste0("acf_profile_", group, ".tsv")))
        for (type in c("profile", "heatmap")) {
            .acf_save_plot(plot_autocorrelation(x, type, group),
                file.path(directory, paste0("acf_", type, "_", group, ".pdf")),
                if (type == "heatmap") 10 else 4)
        }
    }
    invisible(normalizePath(directory))
}

.acf_cluster_view <- function(x) {
    cd <- SummarizedExperiment::colData(x)
    if (!"acf_cluster" %in% names(cd)) stop("No ACF clustering in this result")
    keep <- which(!is.na(cd$acf_cluster) & nzchar(cd$acf_cluster))
    x <- x[, keep, drop = FALSE]
    SummarizedExperiment::colData(x)$cluster <- factor(cd$acf_cluster[keep])
    x
}

.acf_validate_options <- function(options) {
    allowed <- c("max_lag", "clustering", "n_pcs", "n_neighbors", "resolution", "seed")
    if (!is.list(options) || (length(options) && (is.null(names(options)) ||
        anyNA(names(options)) || any(!names(options) %in% allowed) || anyDuplicated(names(options)))))
        stop("autocorrelation_options must be a uniquely named list of locus_autocorrelation arguments")
    for (name in intersect(names(options), c("max_lag", "n_pcs", "n_neighbors", "seed", "resolution"))) {
        v <- options[[name]]
        if (name == "max_lag" && is.null(v)) next
        lower <- if (name %in% c("n_pcs", "n_neighbors")) 2 else 0
        if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < lower ||
            (name != "resolution" && (v != floor(v) || v > .Machine$integer.max)) ||
            (name == "resolution" && v == 0)) stop("Invalid ACF option: ", name)
    }
    if ("clustering" %in% names(options) && (!is.logical(options$clustering) ||
        length(options$clustering) != 1L || is.na(options$clustering)))
        stop("ACF clustering must be TRUE or FALSE")
    invisible(TRUE)
}

.add_locus_autocorrelation <- function(x, options) {
    md <- S4Vectors::metadata(x)
    # A failure in optional ACF analysis must not discard the locus or its clusters.
    result <- tryCatch(do.call(locus_autocorrelation, c(list(x = x), options)),
        error = function(e) e)
    if (inherits(result, "error")) {
        md$autocorrelation_status <- paste0("failed: ", conditionMessage(result))
        warning(md$autocorrelation_status, call. = FALSE)
    } else {
        md$autocorrelation <- result
        md$autocorrelation_status <- "complete"
    }
    S4Vectors::metadata(x) <- md
    x
}

utils::globalVariables(c("lag_bp", "mean_acf"))
