# Adapted from native R auto_correlation steps 02, 03 and repeat_peak in 04.
# Mean-centered, variance-normalized linear ACF; columns are consecutive bases.
# Lag k uses sum(x_centered[t] * x_centered[t+k]) / sum(x_centered^2).
.acf_compute <- function(binary, n_features = NULL, positions = NULL, window = NULL) {
  x <- as.matrix(binary)
  if (length(dim(binary)) != 2L || ncol(x) < 2L)
    stop("Expected reads x consecutive genomic bases")
  if (!is.null(window)) {
    if (is.null(positions)) stop("Window selection requires a position for each input column")
    if (length(positions) != ncol(x) || length(window) != 2L ||
        anyNA(c(positions, window)) || !all(diff(positions) == 1) ||
        window[1L] > window[2L] || !all(window %in% positions))
      stop("Expected consecutive 1-bp positions and covered inclusive window bounds")
    x <- x[, positions >= window[1L] & positions <= window[2L], drop = FALSE]
    if (ncol(x) < 2L) stop("Autocorrelation requires at least two bases")
  }
  if (is.null(n_features)) n_features <- ncol(x)
  if (length(n_features) != 1L || !is.numeric(n_features) || !is.finite(n_features) ||
      n_features != floor(n_features) || n_features < 1L || n_features > ncol(x))
    stop("n_features must be an integer in [1, read width]")
  if (!is.numeric(x) && !is.logical(x)) stop("Input must contain only fully covered binary 0/1 calls")
  if (any(!is.finite(x)) || !all(x %in% c(0, 1)))
    stop("Input must contain only fully covered binary 0/1 calls")
  valid <- unname(rowSums(x) > 0 & rowSums(x) < ncol(x))
  profiles <- matrix(NA_real_, nrow(x), n_features)
  for (i in which(valid)) {
    # stats::acf uses the fixed lag-zero denominator, without overlap correction.
    profiles[i, ] <- as.numeric(stats::acf(x[i, ], lag.max = n_features - 1L,
      type = "correlation", plot = FALSE, demean = TRUE)$acf)
  }
  list(profiles = profiles, valid = valid)
}

# Centered PCA -> exact correlation neighbors -> fuzzy graph -> Leiden.
# UMAP is for display only. All computation uses native R packages.
.acf_correlation_neighbors <- function(scores, k, block_size = 256L) {
  centered <- scores - rowMeans(scores)
  unit <- centered / sqrt(rowSums(centered^2))
  n <- nrow(unit)
  stopifnot(all(is.finite(unit)), k >= 2L, k < n)
  indices <- matrix(0L, n, k)
  distances <- matrix(0, n, k)
  # Block multiplication avoids allocating a dense n x n distance matrix.
  for (first in seq.int(1L, n, by = block_size)) {
    rows <- seq.int(first, min(n, first + block_size - 1L))
    d <- 1 - tcrossprod(unit[rows, , drop = FALSE], unit)
    d[] <- pmax(0, pmin(2, d))
    for (j in seq_along(rows)) {
      i <- rows[j]
      d[j, i] <- Inf
      other <- order(d[j, ], seq_len(n))[seq_len(k - 1L)]
      indices[i, ] <- c(i, other)  # uwot counts self among k neighbors
      distances[i, ] <- c(0, d[j, other])
    }
  }
  list(idx = indices, dist = distances)
}

.acf_cluster_profiles <- function(profiles, valid, n_pcs = 50L, n_neighbors = 10L,
                             resolution = 0.4, seed = 0L) {
  stopifnot(is.matrix(profiles), is.logical(valid), length(valid) == nrow(profiles),
    !anyNA(valid), all(is.finite(profiles[valid, , drop = FALSE])),
    length(n_pcs) == 1L, is.finite(n_pcs), n_pcs >= 2L, n_pcs == floor(n_pcs),
    length(n_neighbors) == 1L, is.finite(n_neighbors), n_neighbors >= 2L,
    n_neighbors == floor(n_neighbors), length(resolution) == 1L,
    is.finite(resolution), resolution > 0, length(seed) == 1L, is.finite(seed))
  n <- nrow(profiles)
  result <- list(labels = rep("", n),
    status = ifelse(valid, "insufficient_reads", "zero_variance"),
    embedding = matrix(NA_real_, n, 2L),
    edges = data.frame(source = integer(), target = integer(), weight = numeric()),
    info = list(method = "R PCA -> correlation neighbors -> uwot fuzzy graph -> igraph Leiden; UMAP for display",
      requested_pcs = n_pcs, requested_neighbors = n_neighbors,
      resolution = resolution, seed = seed, n_clustered = 0L, n_clusters = 0L))
  indices <- which(valid)
  if (length(indices) < 4L) return(result)
  values <- profiles[indices, , drop = FALSE]
  if (all(abs(sweep(values, 2L, values[1L, ], "-")) <= 1e-12)) {
    result$status[indices] <- "identical_profiles"
    return(result)
  }
  pcs <- min(n_pcs, nrow(values) - 1L, ncol(values) - 1L)
  if (pcs < 2L) stop("At least three lag features are required for clustering")
  # Retain lag zero, center columns, and do not variance-scale the features.
  set.seed(seed)
  pca <- if (pcs < min(dim(values))/2) {
    irlba::prcomp_irlba(values, n = pcs, center = TRUE, scale. = FALSE, tol = 1e-8)
  } else {
    stats::prcomp(values, center = TRUE, scale. = FALSE, rank. = pcs)
  }
  total_variance <- sum(apply(values, 2L, stats::var))
  scores <- pca$x[, seq_len(pcs), drop = FALSE]
  # Correlation between score rows depends on component signs. Use a stable
  # loading-based convention: each component's largest absolute loading is positive.
  signs <- vapply(seq_len(pcs), function(j)
    sign(pca$rotation[which.max(abs(pca$rotation[,j])),j]), numeric(1))
  scores <- sweep(scores, 2L, signs, "*")
  usable <- apply(scores, 1L, function(x) all(is.finite(x)) && stats::sd(x) > 1e-12)
  result$status[indices[!usable]] <- "undefined_pca_correlation"
  indices <- indices[usable]
  scores <- scores[usable, , drop = FALSE]
  if (length(indices) < 4L) return(result)
  k <- min(n_neighbors, length(indices) - 1L)
  nn <- .acf_correlation_neighbors(scores, k)
  set.seed(seed)
  embedding <- uwot::umap(scores, nn_method = nn, n_neighbors = k,
    n_components = 2L, init = "random", min_dist = 0.5, spread = 1,
    n_threads = 1L, n_sgd_threads = 1L, fast_sgd = FALSE,
    ret_extra = "fgraph", verbose = FALSE)
  adjacency <- embedding$fgraph
  stopifnot(all(is.finite(adjacency@x)), all(is.finite(embedding$embedding)))
  graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected",
    weighted = TRUE, diag = FALSE)
  set.seed(seed)
  partition <- igraph::cluster_leiden(graph, objective_function = "modularity",
    resolution = resolution, weights = igraph::E(graph)$weight, n_iterations = -1L)
  membership <- igraph::membership(partition)
  # Stable labels by decreasing size; report loader adds 1 for display.
  sizes <- table(membership)
  groups <- names(sizes)[order(-as.integer(sizes), as.integer(names(sizes)))]
  result$labels[indices] <- as.character(match(as.character(membership), groups) - 1L)
  result$status[indices] <- "clustered"
  result$embedding[indices, ] <- embedding$embedding
  edges <- as.data.frame(igraph::as_data_frame(graph, what = "edges"))
  result$edges <- data.frame(source = indices[as.integer(edges$from)],
    target = indices[as.integer(edges$to)], weight = edges$weight)
  result$info <- c(result$info[!names(result$info) %in% c("n_clustered", "n_clusters")],
    list(actual_pcs = pcs, actual_neighbors = k, n_clustered = length(indices),
      n_clusters = length(groups), pca_variance_explained = sum(pca$sdev[seq_len(pcs)]^2)/total_variance,
      neighbor_metric = "correlation", neighbor_representation = "PCA scores",
      leiden_flavor = "igraph", leiden_directed = FALSE,
      package_versions = c(igraph = as.character(utils::packageVersion("igraph")),
                           uwot = as.character(utils::packageVersion("uwot")))))
  result
}

# Match scipy.signal.find_peaks: midpoint (rounded down) for a flat local peak;
# endpoints excluded; equal-height peaks resolve to the first/lower lag.
.acf_repeat_peak <- function(profile, lower = 140L, upper = 250L) {
  missing <- c(lag = NA_real_, height = NA_real_)
  if (length(profile) <= upper + 1L || any(!is.finite(profile))) return(missing)
  runs <- rle(profile)
  if (length(runs$values) < 3L) return(missing)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  middle <- seq.int(2L, length(runs$values) - 1L)
  candidates <- middle[runs$values[middle] > runs$values[middle - 1L] &
                         runs$values[middle] > runs$values[middle + 1L]]
  indices <- floor((starts[candidates] + ends[candidates]) / 2)
  lags <- indices - 1L
  indices <- indices[lags >= lower & lags <= upper & profile[indices] > 0]
  if (!length(indices)) return(missing)
  peak <- indices[which.max(profile[indices])]
  c(lag = peak - 1L, height = profile[peak])
}
