.check_locus <- function(x) {
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

locus_profile <- function(x, by_cluster = FALSE) {
    .check_locus(x)
    m <- SummarizedExperiment::assay(x, "m6a")
    cd <- SummarizedExperiment::colData(x)
    pos <- BiocGenerics::start(SummarizedExperiment::rowRanges(x))
    groups <- list(All = seq_len(ncol(x)))
    if (by_cluster) {
        if (!"cluster" %in% names(cd)) stop("No clustering in this result")
        groups <- split(seq_len(ncol(x)), droplevels(cd$cluster))
    }
    observed <- Matrix::rowSums(m) > 0
    if (!any(observed)) observed[] <- TRUE
    p <- pos[observed]
    out <- lapply(names(groups), function(g) {
        idx <- groups[[g]]
        met <- as.numeric(Matrix::rowSums(m[observed, idx, drop = FALSE]))
        frac <- if (length(idx)) met / length(idx) else rep(NA_real_, length(p))
        S4Vectors::DataFrame(pos = p, base = rep("A", length(p)),
            group = rep(g, length(p)), cov = rep(length(idx), length(p)),
            met = met, frac = frac,
            smooth_cov = rep(length(idx), length(p)),
            smooth_frac = smooth_mean_slidewindow(frac, p,
                S4Vectors::metadata(x)$parameters$smoothing_bp))
    })
    do.call(rbind, out)
}

.read_order <- function(x, clustered) {
    cd <- SummarizedExperiment::colData(x)
    samples <- S4Vectors::metadata(x)$samples$sample_name
    if (clustered && !"cluster" %in% names(cd))
        stop("No clustering in this result")
    if (clustered) order(cd$cluster, match(cd$sample_name, samples),
        cd$start, cd$RID) else
        order(match(cd$sample_name, samples), cd$start, cd$RID)
}

plot_locus <- function(x, type = c("heatmap", "profile", "footprints"),
    clustered = FALSE) {
    .check_locus(x)
    type <- match.arg(type)
    if (!is.logical(clustered) || length(clustered) != 1L || is.na(clustered))
        stop("clustered must be TRUE or FALSE")
    if (!nrow(x) || !ncol(x)) {
        return(ggplot2::ggplot() + ggplot2::theme_void() +
            ggplot2::labs(title = "No full-span reads or positions"))
    }
    if (type == "heatmap") return(.plot_heatmap(x, clustered))
    if (type == "footprints") return(.plot_footprints(x, clustered))
    .plot_profile(x, clustered)
}

.plot_profile <- function(x, clustered, stacked = FALSE) {
    d <- as.data.frame(locus_profile(x, by_cluster = clustered))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = pos, y = smooth_frac))
    if (clustered) {
        p <- p + ggplot2::geom_line(ggplot2::aes(colour = group), linewidth = 0.6) +
            ggplot2::scale_colour_manual(
                values = cluster_colors(SummarizedExperiment::colData(x)$cluster),
                name = "Cluster")
        subtitle <- sprintf("%d full-span reads; %d bp sliding mean within each cluster; dashed = TSS",
            ncol(x), S4Vectors::metadata(x)$parameters$smoothing_bp)
    } else {
        if (!stacked) p <- p + ggplot2::geom_line(ggplot2::aes(y = frac),
            colour = "grey75", linewidth = 0.3)
        p <- p + ggplot2::geom_line(colour = "blue",
            linewidth = if (stacked) 0.5 else 0.6)
        subtitle <- sprintf("%d full-span reads; grey = per position, blue = %d bp sliding mean; dashed = TSS",
            ncol(x), S4Vectors::metadata(x)$parameters$smoothing_bp)
    }
    p <- p + theme_fiberseq("bw", base_size = if (stacked) 11 else 10) +
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
        ggplot2::coord_cartesian(xlim = .plot_bounds(x)) +
        ggplot2::labs(x = .axis_label(x), y = "m6A fraction of reads",
            title = paste0(.locus_title(x), ": aggregated accessibility"),
            subtitle = subtitle)
    .add_tss(p, x, linewidth = 0.5)
}

.plot_heatmap <- function(x, clustered) {
    ord <- .read_order(x, clustered)
    cd <- SummarizedExperiment::colData(x)[ord, , drop = FALSE]
    m <- as.matrix(Matrix::t(SummarizedExperiment::assay(x, "m6a")[,
        ord, drop = FALSE]))
    colors <- list(sample = .sample_colors(x))
    annotation <- data.frame(sample = factor(cd$sample_label,
        levels = names(colors$sample)))
    if (clustered) {
        annotation$cluster <- cd$cluster
        colors$cluster <- cluster_colors(cd$cluster)
    }
    positions <- BiocGenerics::start(SummarizedExperiment::rowRanges(x))
    ticks <- unique(round(seq(1, ncol(m), length.out = 5L)))
    bottom <- list(coordinate = ComplexHeatmap::anno_mark(at = ticks,
        labels = format(positions[ticks], big.mark = ",", scientific = FALSE),
        which = "column", side = "bottom", labels_gp = grid::gpar(fontsize = 8)))
    tss_col <- match(.plot_tss(x), positions)
    if (!is.na(tss_col)) bottom$TSS <- ComplexHeatmap::anno_mark(at = tss_col,
        labels = "TSS", which = "column", side = "bottom",
        labels_gp = grid::gpar(fontsize = 8, col = "#D55E00"))
    bottom$annotation_name_gp <- grid::gpar(fontsize = 8)
    ## Draw once per cluster slice, also when callers use draw() themselves.
    tss_line <- function(j, i, x, y, width, height, fill) {
        if (!is.na(tss_col)) grid::grid.lines(
            x = rep((tss_col - 0.5) / ncol(m), 2), y = c(0, 1),
            gp = grid::gpar(col = "#D55E00", lty = "dashed", lwd = 1.2))
    }
    ComplexHeatmap::Heatmap(m, name = "m6A",
        col = c("0" = "white", "1" = "black"),
        heatmap_legend_param = list(at = c(0, 1), labels = c("no m6A call", "m6A")),
        cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
        row_split = if (clustered) cd$cluster else NULL,
        row_gap = grid::unit(0.6, "mm"), row_title_rot = 0,
        row_title_gp = grid::gpar(fontsize = 8),
        show_row_names = FALSE, show_column_names = FALSE,
        use_raster = TRUE, raster_quality = 2, raster_resize_mat = FALSE,
        raster_by_magick = FALSE, layer_fun = tss_line,
        left_annotation = ComplexHeatmap::rowAnnotation(df = annotation, col = colors,
            annotation_legend_param = list(sample = list(ncol = 2))),
        bottom_annotation = do.call(ComplexHeatmap::HeatmapAnnotation, bottom),
        column_title = paste0(.locus_title(x), ", ", .axis_label(x), ":",
            paste(.plot_bounds(x), collapse = "-"), "\n",
            .plot_subtitle(x, clustered)),
        column_title_gp = grid::gpar(fontsize = 11))
}

.plot_footprints <- function(x, clustered) {
    ord <- .read_order(x, clustered)
    cd <- SummarizedExperiment::colData(x)[ord, , drop = FALSE]
    tracks <- unlist(cd$footprints, use.names = FALSE)
    ids <- rev(cd$RID)
    cd$RID <- factor(cd$RID, levels = ids)
    bounds <- .plot_bounds(x)
    group <- if (clustered) cd$cluster else rep("All", nrow(cd))
    scaffold <- data.frame(RID = cd$RID, cluster = group,
        start = rep(bounds[1], nrow(cd)), end = rep(bounds[2], nrow(cd)))
    calls <- Matrix::summary(SummarizedExperiment::assay(x, "m6a")[, ord, drop = FALSE])
    positions <- BiocGenerics::start(SummarizedExperiment::rowRanges(x))
    calls <- data.frame(pos = positions[calls$i], RID = cd$RID[calls$j],
        cluster = group[calls$j])
    reads <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = scaffold, ggplot2::aes(
            x = start, xend = end, y = RID, yend = RID),
            colour = "grey50", linewidth = 0.5) +
        ggplot2::geom_tile(data = calls, ggplot2::aes(x = pos, y = RID),
            fill = "blue", width = 1, height = 0.9)
    reads <- .read_panel(reads, x, cd, clustered, "m6A calls") +
        ggplot2::labs(x = NULL)
    footprints <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = scaffold, ggplot2::aes(
            x = start, xend = end, y = RID, yend = RID),
            colour = "grey70", linewidth = 0.5)
    if (length(tracks)) {
        d <- data.frame(start = pmax(bounds[1], BiocGenerics::start(tracks)),
            end = pmin(bounds[2], BiocGenerics::end(tracks)),
            as.data.frame(S4Vectors::mcols(tracks)))
        d$RID <- factor(d$RID, levels = ids)
        d$cluster <- group[match(d$RID, cd$RID)]
        ## Match plot_region_panels(): FIRE, nucleosomes, then TF on each read.
        for (kind in c("fire", "nucleosome", "tf")) {
            footprints <- footprints + ggplot2::geom_segment(
                data = d[d$track == kind & d$end >= d$start, , drop = FALSE],
                ggplot2::aes(x = start, xend = end, y = RID, yend = RID,
                    colour = class), linewidth = 0.5)
        }
        p <- S4Vectors::metadata(x)$parameters
        nuc_label <- paste0(p$nuc_min, "-", p$nuc_max, " bp nucleosome")
        colors <- c(fire = "#FF8C00", linker = "darkgray",
            nucleosome = "darkblue", "size10-30" = "purple",
            "size40-60" = "magenta", "size60-80" = "#41c4e1")
        colors[nuc_label] <- "darkblue"
        footprints <- footprints + ggplot2::scale_colour_manual(values = colors,
            name = "class", guide = ggplot2::guide_legend(order = 2))
    }
    footprints <- .read_panel(footprints, x, cd, clustered,
        "FIRE / nucleosomes / TF footprints")
    profile <- .plot_profile(x, clustered, stacked = TRUE) +
        ggplot2::labs(x = NULL, y = "m6A fraction", title = .locus_title(x),
            subtitle = .plot_subtitle(x, clustered))
    ## Align genomic positions across all three panels; only read panels
    ## contain sample and cluster strips in this left margin.
    profile$coordinates <- ggplot2::coord_cartesian(
        xlim = .annotation_layout(x, clustered)$limits, expand = FALSE)
    patchwork::wrap_plots(profile, reads, footprints, ncol = 1,
        heights = c(2, 6, 6), guides = "collect")
}

## Names used by ggplot2 and unchanged shared dplyr expressions.
utils::globalVariables(c("pos", "frac", "smooth_frac", "group", "start",
    "end", "RID", "class", "track", "read_start", "read_end"))
