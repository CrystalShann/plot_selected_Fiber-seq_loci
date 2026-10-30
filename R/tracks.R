.read_tracks <- function(samples, reads, region, parameters) {
    kinds <- intersect(c("nucleosome", "tf", "fire"), names(samples))
    kinds <- kinds[vapply(kinds, function(k)
        any(!is.na(samples[[k]]) & nzchar(samples[[k]])), logical(1))]
    tracks <- lapply(kinds, function(kind) {
        pieces <- lapply(seq_len(nrow(samples)), function(i)
            .read_track_sample(samples, i, reads, region, kind, parameters))
        if (length(pieces)) do.call(c, pieces) else GenomicRanges::GRanges()
    })
    names(tracks) <- kinds
    empty <- intersect(kinds[lengths(tracks) == 0L], c("nucleosome", "tf"))
    if (length(empty))
        warning("no footprints called in the selected region", call. = FALSE)
    all <- if (length(tracks)) unlist(GenomicRanges::GRangesList(tracks),
        use.names = FALSE) else GenomicRanges::GRanges()
    if (length(all)) {
        by_read <- S4Vectors::split(all, factor(S4Vectors::mcols(all)$RID,
            levels = reads$RID))
    } else {
        by_read <- GenomicRanges::GRangesList(stats::setNames(
            rep(list(GenomicRanges::GRanges()), nrow(reads)), reads$RID))
    }
    list(by_read = by_read, supplied = kinds, empty = empty)
}

.read_track_sample <- function(samples, i, reads, region, kind, p) {
    path <- samples[[kind]][i]
    empty <- GenomicRanges::GRanges()
    if (is.na(path) || !nzchar(path)) return(empty)
    chr <- as.character(GenomicRanges::seqnames(region))
    path <- .expand_path(path, samples$sample_name[i], chr)
    if (!chr %in% Rsamtools::seqnamesTabix(path)) return(empty)
    format <- samples[[paste0(kind, "_format")]][i]
    rows <- read_tabix_region(path, region)
    if (!nrow(rows)) return(empty)
    expected <- c(bed12_fibertools = 12L, bed13_fiberhmm = 13L,
        bed15_fiberhmm_tf = 15L, bed6_per_sample = 6L, bed11_fire = 11L)
    if (ncol(rows) != expected[[format]])
        stop("Wrong column count for ", format, ": ", path)
    selected <- reads$original_RID[reads$sample_name == samples$sample_name[i]]
    rows <- rows[as.character(rows[[4]]) %in% selected, , drop = FALSE]
    if (!nrow(rows)) return(empty)
    if (format %in% c("bed12_fibertools", "bed13_fiberhmm", "bed15_fiberhmm_tf")) {
        b <- convert_ft_bed12_to_bed6(rows, format = format, source = path,
            keep_block_scores = format != "bed12_fibertools",
            validate_blocks = format != "bed12_fibertools")
    } else {
        b <- data.frame(start = as.numeric(rows[[2]]),
            end = as.numeric(rows[[3]]), RID = as.character(rows[[4]]))
    }
    if (anyNA(b$start) || anyNA(b$end) ||
        any(b$start < 0 | b$end <= b$start | b$start != floor(b$start) |
            b$end != floor(b$end))) stop("Invalid track coordinates: ", path)
    b <- b[b$end >= BiocGenerics::start(region) &
        b$start < BiocGenerics::end(region), , drop = FALSE]
    if (!nrow(b)) return(empty)
    b$class <- kind
    size <- b$end - b$start
    if (kind == "nucleosome") {
        b <- b[size >= p$nuc_min & size <= p$nuc_max, , drop = FALSE]
        b$class <- rep(paste0(p$nuc_min, "-", p$nuc_max,
            " bp nucleosome"), nrow(b))
    } else if (kind == "tf") {
        if (format == "bed15_fiberhmm_tf") {
            if (anyNA(b$block_score)) stop("Invalid TF block scores: ", path)
            b <- b[b$block_score >= p$min_tf_score, , drop = FALSE]
        }
        size <- b$end - b$start
        b$class <- ifelse(size >= 10 & size < 30, "size10-30",
            ifelse(size >= 40 & size < 60, "size40-60",
                ifelse(size >= 60 & size < 80, "size60-80", NA_character_)))
        b <- b[!is.na(b$class), , drop = FALSE]
    }
    if (!nrow(b)) return(empty)
    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(b$start + 1L, b$end),
        seqinfo = GenomeInfoDb::seqinfo(region))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        RID = paste(samples$sample_name[i], b$RID, sep = "::"),
        track = rep(kind, nrow(b)), class = b$class)
    gr
}
