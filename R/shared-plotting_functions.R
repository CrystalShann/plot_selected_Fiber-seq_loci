# Unchanged selections from code/parsing_functions/plotting_functions.r.
# Source ranges and SHA256 checksums are recorded in inst/SHARED_*.

LEIDEN_CLUSTER_COLORS <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
  "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
  "green1", "yellow4", "yellow3", "darkorange4", "brown"
)

cluster_palette <- function(levels) {
  n <- length(levels)
  cols <- if (n <= length(LEIDEN_CLUSTER_COLORS)) {
    LEIDEN_CLUSTER_COLORS[seq_len(n)]
  } else {
    grDevices::colorRampPalette(LEIDEN_CLUSTER_COLORS)(n)
  }
  stats::setNames(cols, levels)
}

cluster_id_palette <- function(levels) {
  idx <- suppressWarnings(as.integer(sub("^cluster", "", levels)))
  if (anyNA(idx) || any(!grepl("^cluster[0-9]+$", levels))) idx <- seq_along(levels)
  cols <- if (max(idx) <= length(LEIDEN_CLUSTER_COLORS)) LEIDEN_CLUSTER_COLORS[idx] else
    grDevices::colorRampPalette(LEIDEN_CLUSTER_COLORS)(max(idx))[idx]
  stats::setNames(cols, levels)
}

theme_fiberseq <- function(base = c("classic", "bw", "cowplot"),
                           base_size = 9, base_family = "", overrides = list()) {
  base <- match.arg(base)
  value <- switch(base,
    classic = ggplot2::theme_classic(base_size = base_size, base_family = base_family),
    bw = ggplot2::theme_bw(base_size = base_size, base_family = base_family),
    cowplot = cowplot::theme_cowplot(font_size = base_size, font_family = base_family))
  if (length(overrides)) value <- value + do.call(ggplot2::theme, overrides)
  value
}
