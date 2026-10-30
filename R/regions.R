# All region tables use BED coordinates 
# Optionally save a table as a TSV file
.write_table <- function(x, output) {
  if (!is.null(output)) {
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(x, output, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  }
  x
}
# Read a headered TSV 
.read_table <- function(path) utils::read.delim(path, stringsAsFactors = FALSE,
                                              check.names = FALSE, na.strings = "")
## rtracklayer handles GTF quoting and GFF3 attribute encoding.
# Extract transcript coordinates and gene relationships
.annotation <- function(path, canonical_attribute = NULL, canonical_value = NULL) {
    gr <- rtracklayer::import(path)
    mc <- S4Vectors::mcols(gr)
    gff <- "ID" %in% names(mc)
 
    attr <- function(key) {
        if (!key %in% names(mc)) return(rep("", length(gr)))
        vapply(as.list(mc[[key]]), paste, character(1), collapse = ",")
    }
    type <- as.character(mc$type)
    tx <- type %in% c("mRNA", "transcript", "ncRNA", "lnc_RNA", "tRNA", "rRNA")
    if (!any(tx)) stop("Annotation requires transcript/mRNA feature rows")
    id <- attr(if (gff) "ID" else "transcript_id")
    gene <- attr(if (gff) "Parent" else "gene_id")
    symbol <- attr(if (gff) "Name" else "gene_name")
    if (gff) {
        genes <- which(type == "gene")
        symbol <- unname(stats::setNames(symbol[genes], id[genes])[gene])
    }
    d <- data.frame(chr = as.character(GenomicRanges::seqnames(gr)),
        start1 = BiocGenerics::start(gr), end1 = BiocGenerics::end(gr),
        strand = as.character(BiocGenerics::strand(gr)), gene_id = gene,
        gene_name = symbol, transcript_id = id, stringsAsFactors = FALSE)[tx, ]
    unnamed <- is.na(d$gene_name) | !nzchar(d$gene_name)
    d$gene_name[unnamed] <- d$gene_id[unnamed]
    if (anyNA(d$transcript_id) || any(!nzchar(d$transcript_id)) ||
        anyNA(d$gene_id) || any(!nzchar(d$gene_id)) ||
        anyDuplicated(d$transcript_id))
        stop("Transcript IDs must be unique and have gene IDs/Parents")
    if (any(!d$strand %in% c("+", "-"))) stop("Transcript strand must be + or -")
    d$tss <- ifelse(d$strand == "+", d$start1, d$end1)
    if (!is.null(canonical_attribute)) {
        values <- strsplit(attr(canonical_attribute)[tx], ",", fixed = TRUE)
        # Flag transcripts carrying the requested canonical tag.
        d$canonical <- vapply(values, function(x) canonical_value %in% x,
            logical(1))
    }
    d
}
# Select requested genes or transcripts
.select_annotation <- function(d, genes, transcripts) {
  if (!is.null(genes)) {
    missing <- setdiff(genes, c(d$gene_id, d$gene_name))
    if (length(missing)) stop("Unknown genes: ", paste(missing, collapse = ", "))
    d <- d[d$gene_id %in% genes | d$gene_name %in% genes, ]
  }
  if (!is.null(transcripts)) {
    missing <- setdiff(transcripts, d$transcript_id)
    if (length(missing)) stop("Unknown transcripts: ", paste(missing, collapse = ", "))
    d <- d[d$transcript_id %in% transcripts, ]
  }
  if (!nrow(d)) stop("No selected transcripts")
  d
}
# Build promoter windows around distinct gene TSSs
.tss_regions <- function(d, flank_bp) {
  if (length(flank_bp) != 1 || is.na(flank_bp) || flank_bp < 0 || flank_bp != floor(flank_bp)) stop("Invalid flank_bp")
  key <- paste(d$gene_id, d$chr, d$tss, d$strand, sep = "|")
  # Merge transcripts sharing a gene, strand and TSS
  rows <- lapply(unique(key), function(k) {
    z <- d[key == k, ]; r <- z[1, ]
    data.frame(region_id = paste(r$gene_id, r$chr, r$tss, r$strand, sep = "_"),
      chr = r$chr, start = r$tss - 1 - flank_bp,
      end = r$tss - 1 + max(1, flank_bp), strand = r$strand, tss = r$tss,
      gene_id = r$gene_id, gene_name = r$gene_name,
      transcript_ids = paste(sort(unique(z$transcript_id)), collapse = ","))
  })
  do.call(rbind, rows)
}

# Validate coordinates and standardize region metadata
.normalize_regions <- function(d) {
  if (!all(c("chr", "start", "end") %in% names(d))) stop("Regions require chr, start, end")
  d$start <- as.numeric(d$start); d$end <- as.numeric(d$end)
  if (!nrow(d) || anyNA(d[, c("chr", "start", "end")]) || any(!nzchar(d$chr)) ||
      any(!is.finite(d$start) | !is.finite(d$end) | d$start != floor(d$start) |
          d$end != floor(d$end) | d$start < 0 | d$end <= d$start)) stop("Invalid region coordinates (or promoter extends before position 1)")
  if (!"region_id" %in% names(d)) d$region_id <- paste(d$chr, d$start, d$end, sep = "_")
  unnamed <- is.na(d$region_id) | d$region_id %in% c("", ".")
  d$region_id[unnamed] <- paste(d$chr[unnamed], d$start[unnamed], d$end[unnamed], sep = "_")
  d$region_id <- as.character(d$region_id)
  if (anyNA(d$region_id) || any(!nzchar(d$region_id)) || anyDuplicated(d$region_id) ||
      any(d$region_id %in% c(".", ".."))) stop("Region IDs must be unique safe names")
  if (!"strand" %in% names(d)) d$strand <- "*"
  if (!all(d$strand %in% c("+", "-", "*", "."))) stop("Invalid strand")
  if (!"tss" %in% names(d)) d$tss <- NA_real_
  d$tss <- as.numeric(d$tss)
  anchors <- !is.na(d$tss)
  if (any(!is.finite(d$tss[anchors]) | d$tss[anchors] != floor(d$tss[anchors]) |
          d$tss[anchors] <= d$start[anchors] | d$tss[anchors] > d$end[anchors])) stop("TSS must be an integer inside its region")
  d$analysis_start <- d$start + 1; d$analysis_end <- d$end
  d$width <- d$end - d$start
  d
}

# Convert coordinates or annotations to named genomic ranges
prepare_regions <- function(regions = NULL, output = NULL, coordinate_system = NULL,
                            annotation = NULL, genes = NULL, transcripts = NULL,
                            flank_bp = 1000, bed_annotation = c("plain", "gencode_custom")) {
  if (methods::is(regions, "GRanges")) {
    if (!is.null(annotation)) stop("Supply regions or annotation, not both")
    gr <- .regions_gr(.regions_table(regions), output)
    GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::seqlevels(regions)
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(regions)
    return(gr)
  }
  bed_annotation <- match.arg(bed_annotation)
  if (!is.null(annotation)) {
    if (!is.null(regions)) stop("Supply regions or annotation, not both")
    d <- .tss_regions(.select_annotation(.annotation(annotation), genes, transcripts), flank_bp)
  } else if (is.data.frame(regions) || methods::is(regions, "DataFrame")) {
    if (is.null(coordinate_system) || !coordinate_system %in% c("bed", "one_based"))
      stop("Coordinate tables need coordinate_system = 'bed' or 'one_based'")
    d <- as.data.frame(regions)
    if (coordinate_system == "one_based") d$start <- as.numeric(d$start) - 1
  } else if (length(regions) == 1 && file.exists(regions)) {
    con <- if (grepl("\\.gz$", regions)) gzfile(regions, open = "rt") else file(regions, open = "rt")
    on.exit(close(con))
    d <- utils::read.delim(con, header = FALSE, quote = "", comment.char = "#", stringsAsFactors = FALSE)
    if (ncol(d) < 3) stop("BED requires at least three columns")
    names(d)[seq_len(3)] <- c("chr", "start", "end")
    if (ncol(d) >= 4) names(d)[4] <- "region_id"
    if (ncol(d) >= 6) names(d)[6] <- "strand"
    if (bed_annotation == "gencode_custom") {
      if (!"region_id" %in% names(d) || any(d$end - d$start != 1)) stop("Custom TSS BED requires named single-base rows")
      parts <- strsplit(d$region_id, ";", fixed = TRUE)
      if (any(lengths(parts) < 3)) stop("Custom GENCODE names require gene ID;annotation;gene name")
      d$gene_id <- vapply(parts, `[`, character(1), 1)
      d$gene_name <- vapply(parts, `[`, character(1), 3)
      if (!is.null(genes)) d <- d[d$gene_id %in% genes | d$gene_name %in% genes, ]
      d$tss <- d$start + 1
      d$region_id <- paste(d$gene_id, d$chr, d$tss, sep = "_")
      d$start <- d$start - flank_bp; d$end <- d$tss - 1 + max(1, flank_bp)
    }
  } else {
    if (!length(regions) || anyNA(regions) || any(!grepl("^[^:]+:[0-9]+-[0-9]+$", regions))) stop("Supply BED, coordinates, or annotation")
    d <- data.frame(chr = sub(":.*", "", regions),
      start = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", regions)) - 1,
      end = as.numeric(sub(".*-", "", regions)))
  }
  .regions_gr(d, output)
}

# Select canonical transcript TSSs or promoter windows
extract_canonical_tss <- function(annotation, output = NULL, genes = NULL,
                                  attribute = "tag", value = "Ensembl_canonical", flank_bp = 0) {
  d <- .select_annotation(.annotation(annotation, attribute, value), genes, NULL)
  counts <- table(factor(d$gene_id[d$canonical], levels = unique(d$gene_id)))
  if (any(counts != 1)) stop("Expected one canonical transcript per gene; missing/ambiguous: ",
                            paste(names(counts)[counts != 1], collapse = ", "))
  .regions_gr(.tss_regions(d[d$canonical, ], flank_bp), output)
}

# Convert a BED-coordinate table to genomic ranges.
.regions_gr <- function(d, output = NULL) {
    d <- .normalize_regions(d)
    .write_table(d, output)
    d$strand[d$strand == "."] <- "*"
    gr <- GenomicRanges::makeGRangesFromDataFrame(d,
        seqnames.field = "chr", starts.in.df.are.0based = TRUE,
        keep.extra.columns = TRUE)
    names(gr) <- d$region_id
    discard <- intersect(c("analysis_start", "analysis_end"),
        names(S4Vectors::mcols(gr)))
    S4Vectors::mcols(gr)[discard] <- NULL
    gr
}

# Convert genomic ranges to a normalized BED-coordinate table
.regions_table <- function(gr) {
    if (!methods::is(gr, "GRanges")) stop("regions must be a GRanges")
    d <- as.data.frame(S4Vectors::mcols(gr))
    d$chr <- as.character(GenomicRanges::seqnames(gr))
    d$start <- BiocGenerics::start(gr) - 1L
    d$end <- BiocGenerics::end(gr)
    d$strand <- as.character(BiocGenerics::strand(gr))
    if (!is.null(names(gr))) d$region_id <- names(gr)
    .normalize_regions(d)
}
