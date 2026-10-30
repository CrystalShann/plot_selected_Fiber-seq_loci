cluster_colors <- function(clusters) {
    labels <- if (is.factor(clusters)) levels(clusters) else
        unique(as.character(clusters))
    if (anyNA(labels) || any(!nzchar(labels)))
        stop("Cluster labels must be nonempty and nonmissing")
    if (!length(labels)) return(stats::setNames(character(), character()))
    if (all(grepl("^cluster[1-9][0-9]*$", labels)))
        cluster_id_palette(labels) else cluster_palette(labels)
}

.sample_colors <- function(x) {
    labels <- unique(S4Vectors::metadata(x)$samples$sample_label)
    stats::setNames(grDevices::hcl.colors(length(labels), "Dynamic"), labels)
}

.axis_label <- function(x) paste(unique(as.character(
    GenomicRanges::seqnames(SummarizedExperiment::rowRanges(x)))), "(1-based)")

.plot_bounds <- function(x) range(BiocGenerics::start(
    SummarizedExperiment::rowRanges(x)))

.plot_tss <- function(x) {
    tss <- S4Vectors::mcols(S4Vectors::metadata(x)$region)$tss
    if (length(tss) == 1L && !is.na(tss)) tss else NA_integer_
}

.add_tss <- function(p, x, linewidth = 0.4) {
    tss <- .plot_tss(x)
    if (!is.na(tss)) p <- p + ggplot2::geom_vline(xintercept = tss,
        linetype = "dashed", colour = "#D55E00", linewidth = linewidth)
    p
}

.locus_title <- function(x) {
    r <- S4Vectors::metadata(x)$region
    paste0(names(r), if (!is.na(.plot_tss(x))) " promoter" else " locus")
}

.plot_subtitle <- function(x, clustered) {
    sprintf("%d full-span reads; ordered by %ssample, read start and read ID; dashed = TSS",
        ncol(x), if (clustered) "cluster, then " else "")
}

## The original tutorial keeps genomic coordinates on both strands.
## Annotation strips occupy a separate margin left of the genomic window.
.annotation_layout <- function(x, clustered) {
    bounds <- .plot_bounds(x)
    span <- max(1, diff(bounds))
    centers <- bounds[1] - span * if (clustered) c(0.07, 0.035) else 0.035
    names(centers) <- if (clustered) c("sample", "cluster") else "sample"
    list(centers = centers, width = span * 0.025,
        limits = c(min(centers) - span * 0.02, bounds[2]))
}

.read_annotations <- function(p, x, cd, clustered) {
    layout <- .annotation_layout(x, clustered)
    data <- data.frame(RID = cd$RID, cluster = if (clustered) cd$cluster else "All",
        annotation_x = rep(layout$centers["sample"], nrow(cd)),
        annotation_key = paste0("Sample: ", cd$sample_label))
    colors <- .sample_colors(x)
    names(colors) <- paste0("Sample: ", names(colors))
    if (clustered) {
        extra <- data
        extra$annotation_x <- layout$centers["cluster"]
        extra$annotation_key <- paste0("Cluster: ", cd$cluster)
        data <- rbind(data, extra)
        cc <- cluster_colors(cd$cluster)
        colors <- c(colors, stats::setNames(cc, paste0("Cluster: ", names(cc))))
    }
    data$RID <- factor(data$RID, levels = levels(cd$RID))
    data$annotation_key <- factor(data$annotation_key, levels = names(colors))
    p + ggplot2::geom_tile(data = data, inherit.aes = FALSE,
        ggplot2::aes(x = annotation_x, y = RID, fill = annotation_key),
        width = layout$width, height = 1) +
        ggplot2::scale_fill_manual(values = colors, breaks = names(colors),
            labels = sub("^(Sample|Cluster): ", "", names(colors)),
            name = if (clustered) "Sample / cluster" else "Sample",
            guide = ggplot2::guide_legend(ncol = 2, order = 1))
}

.read_panel <- function(p, x, cd, clustered, ylabel) {
    if (clustered) p <- p + ggplot2::facet_grid(cluster ~ .,
        scales = "free_y", space = "free_y")
    bounds <- .plot_bounds(x)
    ticks <- pretty(bounds, n = 5)
    ticks <- ticks[ticks >= bounds[1] & ticks <= bounds[2]]
    p <- .read_annotations(p, x, cd, clustered) +
        ggplot2::scale_x_continuous(breaks = ticks,
            labels = function(v) format(v, big.mark = ",", scientific = FALSE)) +
        ggplot2::scale_y_discrete(drop = TRUE) +
        ggplot2::coord_cartesian(xlim = .annotation_layout(x, clustered)$limits,
            expand = FALSE) +
        theme_fiberseq("bw", base_size = 11, overrides = list(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            panel.spacing.y = grid::unit(0.6, "mm"))) +
        ggplot2::labs(x = .axis_label(x), y = ylabel)
    .add_tss(p, x)
}

utils::globalVariables(c("annotation_x", "annotation_key", "cluster"))
