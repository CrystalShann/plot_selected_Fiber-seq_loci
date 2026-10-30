## Native R analysis. The assay is positions x reads throughout the public API.
analyze_loci <- function(samples, regions, clustering,
    window_size = 0, k_neighbors = 10, resolution = 1, seed = 1,
    nuc_min = 130, nuc_max = 160, min_tf_score = 50, smoothing_bp = 10) {
    if (missing(clustering) || !is.logical(clustering) ||
        length(clustering) != 1L || is.na(clustering))
        stop("Choose clustering = TRUE or FALSE explicitly")
    values <- list(window_size, k_neighbors, resolution, seed, nuc_min,
        nuc_max, min_tf_score, smoothing_bp)
    valid <- vapply(values, function(v) is.numeric(v) && length(v) == 1L &&
        is.finite(v) && v >= 0, logical(1))
    if (!all(valid)) stop("Parameters must be finite nonnegative scalars")
    ints <- c(window_size, k_neighbors, seed, smoothing_bp)
    if (any(ints != floor(ints)) || k_neighbors < 1 || smoothing_bp < 1 ||
        nuc_min > nuc_max) stop("Invalid integer parameters or size bounds")
    samples <- make_sample_sheet(samples)
    regions <- prepare_regions(regions)
    validate_inputs(samples, regions)
    parameters <- list(clustering = clustering, window_size = window_size,
        k_neighbors = k_neighbors, resolution = resolution, seed = seed,
        nuc_min = nuc_min, nuc_max = nuc_max, min_tf_score = min_tf_score,
        smoothing_bp = smoothing_bp)
    ans <- lapply(seq_along(regions), function(i)
        .analyze_region(samples, regions[i], parameters))
    names(ans) <- names(regions)
    S4Vectors::SimpleList(ans)
}

.read_m6a_sample <- function(samples, i, region) {
    chr <- as.character(GenomicRanges::seqnames(region))
    path <- .expand_path(samples$m6a[i], samples$sample_name[i], chr)
    bed <- read_ft_bed12(path, region, longest_alignment = TRUE)
    overlap <- nrow(bed)
    if (overlap) bed <- bed[bed$start < BiocGenerics::start(region) &
        bed$end >= BiocGenerics::end(region), , drop = FALSE]
    meta <- data.frame(RID = character(), original_RID = character(),
        chr = character(), start = integer(), end = integer(),
        strand = character(), sample_name = character(), score = numeric())
    calls <- data.frame(RID = character(), pos = integer())
    if (nrow(bed)) {
        original <- as.character(bed$RID)
        bed$RID <- paste(samples$sample_name[i], original, sep = "::")
        meta <- data.frame(RID = bed$RID, original_RID = original,
            chr = bed$chr, start = bed$start + 1L, end = bed$end,
            strand = bed$strand, sample_name = samples$sample_name[i],
            score = bed$score)
        b <- convert_ft_bed12_to_bed6(bed, format = "bed12_fibertools",
            source = path, validate_blocks = FALSE)
        if (any(b$end - b$start != 1L)) stop("Non-single-base m6A block")
        keep <- b$end >= BiocGenerics::start(region) &
            b$end <= BiocGenerics::end(region)
        calls <- unique(data.frame(RID = b$RID[keep], pos = b$end[keep]))
    }
    list(meta = meta, calls = calls, qc = data.frame(
        sample_name = samples$sample_name[i], overlapping_reads = overlap,
        full_span_reads = nrow(meta)))
}

.analyze_region <- function(samples, region, parameters) {
    parts <- lapply(seq_len(nrow(samples)), function(i)
        .read_m6a_sample(samples, i, region))
    reads <- do.call(rbind, lapply(parts, `[[`, "meta"))
    calls <- do.call(rbind, lapply(parts, `[[`, "calls"))
    reads <- reads[order(match(reads$sample_name, samples$sample_name),
        reads$start, reads$RID), , drop = FALSE]
    pos <- seq.int(BiocGenerics::start(region), BiocGenerics::end(region))
    chr <- as.character(GenomicRanges::seqnames(region))
    matrix <- Matrix::sparseMatrix(i = match(calls$pos, pos),
        j = match(calls$RID, reads$RID), x = rep(1, nrow(calls)),
        dims = c(length(pos), nrow(reads)),
        dimnames = list(paste(chr, pos, sep = ":"), reads$RID))
    rr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, width = 1L),
        seqinfo = GenomeInfoDb::seqinfo(region))
    names(rr) <- rownames(matrix)
    cd <- S4Vectors::DataFrame(reads, row.names = reads$RID)
    cd$sample_label <- samples$sample_label[match(cd$sample_name,
        samples$sample_name)]
    cd$read_range <- GenomicRanges::GRanges(reads$chr,
        IRanges::IRanges(reads$start, reads$end),
        strand = sub("^\\.$", "*", reads$strand),
        seqinfo = GenomeInfoDb::seqinfo(region))
    tracks <- .read_tracks(samples, reads, region, parameters)
    cd$footprints <- tracks$by_read
    result <- SummarizedExperiment::SummarizedExperiment(
        assays = list(m6a = matrix), rowRanges = rr, colData = cd,
        metadata = list(region = region, samples = samples,
            qc = S4Vectors::DataFrame(do.call(rbind, lapply(parts, `[[`, "qc"))),
            parameters = parameters, supplied_tracks = tracks$supplied,
            empty_footprint_tracks = tracks$empty,
            package_version = as.character(utils::packageVersion(
                "selectFiberseqLociWorkflow"))))
    if (parameters$clustering) result <- .cluster_locus(result, parameters)
    methods::validObject(result)
    result
}

.cluster_locus <- function(x, p) {
    mat <- SummarizedExperiment::assay(x, "m6a")
    observed <- Matrix::rowSums(mat) > 0
    positions <- BiocGenerics::start(SummarizedExperiment::rowRanges(x))
    features <- if (p$window_size == 0) sum(observed) else
        length(unique((positions[observed] - min(positions)) %/% p$window_size))
    md <- S4Vectors::metadata(x)
    if (ncol(x) < 3L || features < 2L) {
        md$clustering_status <- "skipped: requires >=3 reads and >=2 features"
    } else {
        m <- Matrix::t(mat[observed, , drop = FALSE])
        colnames(m) <- as.character(positions[observed])
        cd <- SummarizedExperiment::colData(x)
        fit <- .fit_clusters(m, positions, p)
        cd$cluster <- fit$assignments$cluster[match(cd$RID, fit$assignments$RID)]
        SummarizedExperiment::colData(x) <- cd
        md$clustering_status <- "complete"
        ## The graph is provenance for the original run, not live subset data.
        md$clustering <- fit
    }
    S4Vectors::metadata(x) <- md
    x
}

## Keep reproducibility local: the caller's RNG state is restored by withr.
.fit_clusters <- function(m, positions, p) {
    binned <- bin_read_matrix(m, min(positions), max(positions), p$window_size)
    knn <- manhattan_knn_graph(binned$mat, k_neighbors = p$k_neighbors,
        verbose = FALSE)
    partition <- withr::with_seed(p$seed, igraph::cluster_leiden(knn$graph,
        objective_function = "modularity", weights = igraph::E(knn$graph)$weight,
        resolution = p$resolution, n_iterations = -1))
    membership <- igraph::membership(partition)[rownames(binned$mat)]
    cluster <- relabel_clusters_by_size(membership)
    assignments <- data.frame(RID = rownames(binned$mat), cluster = cluster)
    list(assignments = assignments, graph = knn$graph,
        profiles = cluster_feature_profiles(binned$mat, cluster),
        feat_mat = binned$mat, window_anno = binned$window_anno,
        n_clusters = nlevels(cluster), params = c(p,
            list(k_eff = knn$k_eff, sigma = knn$sigma,
                mean_knn_dist = knn$mean_knn_dist, quality = partition$quality)))
}
