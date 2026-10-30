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
    if (!nrow(x) || !ncol(x)) {
        return(ggplot2::ggplot() + ggplot2::theme_void() +
            ggplot2::labs(title = "No full-span reads or positions"))
    }
    if (type == "heatmap") return(.plot_heatmap(x, clustered))
    if (type == "footprints") return(.plot_footprints(x, clustered))
    d <- as.data.frame(locus_profile(x, by_cluster = clustered))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = pos, y = smooth_frac,
        colour = group)) + ggplot2::geom_line() + ggplot2::theme_bw() +
        ggplot2::labs(x = .axis_label(x), y = "m6A fraction of full-span reads")
    if (!clustered) p <- p + ggplot2::geom_line(
        ggplot2::aes(y = frac), colour = "grey70", linewidth = 0.3) +
        ggplot2::geom_line() + ggplot2::guides(colour = "none")
    .add_tss(p, x)
}

.axis_label <- function(x) paste(unique(as.character(
    GenomicRanges::seqnames(SummarizedExperiment::rowRanges(x)))), "(1-based)")

.add_tss <- function(p, x) {
    tss <- S4Vectors::mcols(S4Vectors::metadata(x)$region)$tss
    if (length(tss) == 1L && !is.na(tss)) p <- p + ggplot2::geom_vline(
        xintercept = tss, linetype = "dashed", colour = "#D55E00")
    p
}

.plot_heatmap <- function(x, clustered) {
    ord <- .read_order(x, clustered)
    cd <- SummarizedExperiment::colData(x)[ord, , drop = FALSE]
    m <- as.matrix(Matrix::t(SummarizedExperiment::assay(x, "m6a")[,
        ord, drop = FALSE]))
    annotation <- data.frame(sample = cd$sample_label)
    labels <- unique(S4Vectors::metadata(x)$samples$sample_label)
    colors <- list(sample = stats::setNames(
        grDevices::hcl.colors(length(labels), "Dark 3"), labels))
    if (clustered) {
        annotation$cluster <- cd$cluster
        clusters <- levels(cd$cluster)
        colors$cluster <- stats::setNames(
            grDevices::hcl.colors(length(clusters), "Set 2"), clusters)
    }
    positions <- BiocGenerics::start(SummarizedExperiment::rowRanges(x))
    ComplexHeatmap::Heatmap(m, name = "m6A",
        col = c("0" = "white", "1" = "black"), cluster_rows = FALSE,
        cluster_columns = FALSE, show_row_names = FALSE,
        show_column_names = FALSE, use_raster = FALSE,
        left_annotation = ComplexHeatmap::rowAnnotation(df = annotation, col = colors),
        column_title = paste(names(S4Vectors::metadata(x)$region),
            .axis_label(x), paste(range(positions), collapse = "-")))
}

.plot_footprints <- function(x, clustered) {
    ord <- .read_order(x, clustered)
    cd <- SummarizedExperiment::colData(x)[ord, , drop = FALSE]
    tracks <- unlist(cd$footprints, use.names = FALSE)
    ids <- rev(cd$RID)
    rr <- SummarizedExperiment::rowRanges(x)
    bounds <- range(BiocGenerics::start(rr))
    ## Build empty panels even when optional calls are absent.
    d <- data.frame(start = numeric(), end = numeric(), RID = character(),
        class = character(), track = character())
    if (length(tracks)) d <- data.frame(
        start = BiocGenerics::start(tracks), end = BiocGenerics::end(tracks),
        as.data.frame(S4Vectors::mcols(tracks)))
    d$RID <- factor(d$RID, levels = ids)
    scaffold <- data.frame(RID = factor(ids, levels = ids),
        start = rep(bounds[1], length(ids)), end = rep(bounds[2], length(ids)))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = start, xend = end,
        y = RID, yend = RID, colour = class)) +
        ggplot2::geom_segment(data = scaffold, colour = "grey90",
            inherit.aes = FALSE, ggplot2::aes(x = start, xend = end,
                y = RID, yend = RID)) +
        ggplot2::geom_segment(linewidth = 1.8) +
        ggplot2::scale_y_discrete(drop = FALSE) +
        ggplot2::coord_cartesian(xlim = bounds) + ggplot2::theme_bw() +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()) +
        ggplot2::labs(x = .axis_label(x), y = "Full-span reads", colour = "Track")
    if (nrow(d)) p <- p + ggplot2::facet_wrap(~track, ncol = 1)
    .add_tss(p, x)
}

## Names used by ggplot2 and unchanged shared dplyr expressions.
utils::globalVariables(c("pos", "frac", "smooth_frac", "group", "start",
    "end", "RID", "class", "track", "read_start", "read_end"))
