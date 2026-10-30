bin_read_matrix <- function(met_mat, region_start, region_end, window_size = 0) {
  M   <- as.matrix(met_mat)
  pos <- as.integer(colnames(M))

  if (window_size == 0) {
    anno <- data.frame(feature = seq_along(pos), win_start = pos, win_end = pos,
                       mid = pos, n_sites = 1L)
    colnames(M) <- as.character(pos)
    return(list(mat = M, window_anno = anno, window_size = 0))
  }

  starts <- seq(region_start, region_end, by = window_size)
  anno <- data.frame(feature   = seq_along(starts),
                     win_start = starts,
                     win_end   = pmin(starts + window_size - 1, region_end))
  anno$mid <- as.integer(floor((anno$win_start + anno$win_end) / 2))

  win_of_pos <- findInterval(pos, anno$win_start)
  stopifnot(all(win_of_pos >= 1), all(win_of_pos <= nrow(anno)))
  # Counts how many m6A sites fall into each window
  anno$n_sites <- tabulate(win_of_pos, nbins = nrow(anno))

  # sum of each read's site calls per occupied window, then divide by the
  # number of sites in that window -> mean call per window
  win_sums <- t(rowsum(t(M), group = win_of_pos))
  occupied <- as.integer(colnames(win_sums))
  mat <- sweep(win_sums, 2, anno$n_sites[occupied], "/")
  colnames(mat) <- as.character(anno$mid[occupied])
  rownames(mat) <- rownames(M)

  list(mat = mat, window_anno = anno, window_size = window_size)
}

manhattan_knn_graph <- function(mat, k_neighbors = 50, sigma = NULL, verbose = TRUE) {
  n <- nrow(mat)
  if (n < 3) stop("fewer than 3 reads to cluster")
  # determine how many (k) neighbors can be used
  k_eff <- min(k_neighbors, n - 1)
  if (verbose && k_eff < k_neighbors)
    cat(sprintf("KNN graph: k reduced from %d to %d (only %d reads)\n",
                k_neighbors, k_eff, n))

  # calculate the Manhattan distance between every pair of reads
  D <- as.matrix(dist(mat, method = "manhattan"))
  # diagonal contains each read compared with itself, but a read can not count itself as its nearest neighbor
  # change 0 to inf
  diag(D) <- Inf

  # k nearest neighbours of every read
  nn <- t(apply(D, 1, function(d) order(d)[seq_len(k_eff)]))
  # turn the neighbor table into pairs
  from <- rep(seq_len(n), each = k_eff)
  # gets the corresponding nearest neighbor read ID
  to   <- as.vector(t(nn))
  # get distance for each pair
  d_nn <- D[cbind(from, to)]

  
  # convert distance into similarity - use the average distance among nearest neighbor pairs
  if (is.null(sigma)) sigma <- mean(d_nn)
  if (!is.finite(sigma) || sigma <= 0) {
    warning("all KNN distances are 0 (identical reads); using sigma = 1")
    sigma <- 1
  }
  affinity <- exp(-d_nn / sigma)

  # remove duplicate direction edges
  # ex: 1-->3 and 3-->1 is taken as the same representation
  # output is unique undirected pairs --> union KNN graph 
    # if A considers B a neighbor or B considers A a neighbor, A and B get connected
    # does not require both reads to choose each other
  a <- pmin(from, to); b <- pmax(from, to)
  keep <- !duplicated(a * n + b)

  # create the edge table using read ID
  edges <- data.frame(from = rownames(mat)[a[keep]],
                      to   = rownames(mat)[b[keep]],
                      weight = affinity[keep],
                      manhattan_distance = d_nn[keep], stringsAsFactors = FALSE)
  # convert edge table into igraph object
  g <- igraph::graph_from_data_frame(
    edges, directed = FALSE,
    vertices = data.frame(name = rownames(mat), stringsAsFactors = FALSE))

  if (verbose)
    cat(sprintf("KNN graph: %d reads, %d edges, sigma = %.4g, mean KNN distance = %.4g\n",
                n, igraph::ecount(g), sigma, mean(d_nn)))

  list(graph = g, sigma = sigma, k_eff = k_eff,
       mean_knn_dist = mean(d_nn), dist = D)
}

relabel_clusters_by_size <- function(membership) {
  ord <- names(sort(table(membership), decreasing = TRUE))
  factor(paste0("cluster", match(as.character(membership), ord)),
         levels = paste0("cluster", seq_along(ord)))
}

cluster_feature_profiles <- function(feat_mat, cluster, na.rm = TRUE) {
  profiles <- t(sapply(levels(cluster), function(cl)
    colMeans(feat_mat[cluster == cl, , drop = FALSE], na.rm = na.rm)))
  colnames(profiles) <- colnames(feat_mat)
  profiles
}

