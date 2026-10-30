# Read indexed track records overlapping a region.
read_tabix_region <- function(ft_extracted_file, region_gr) {
  tabix_index <- Rsamtools::TabixFile(ft_extracted_file)
  compressed_records <- Rsamtools::scanTabix(tabix_index, param = region_gr)
  # scanTabix gives one character vector of lines per range; pool them
  raw_text <- unlist(compressed_records, use.names = FALSE)

  if (length(raw_text) == 0) {
    return(data.frame())
  }

  read.table(
    text = raw_text,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
}

# Return the column count for a declared format.
footprint_format_columns <- function(format) {
  columns <- c(bed12_fibertools = 12L, bed13_fiberhmm = 13L,
               bed15_fiberhmm_tf = 15L, bed6_per_sample = 6L, bed4_pooled = 4L)
  if (length(format) != 1L || is.na(format) || !format %in% names(columns))
    stop("Unsupported footprint format: ", paste(format, collapse = ", "))
  unname(columns[[format]])
}

# Report a footprint file column-count mismatch.
footprint_column_error <- function(path, format, actual) {
  stop("Footprint file ", path, " (configured format '", format, "'): expected ",
       footprint_format_columns(format), " columns, found ", actual, call. = FALSE)
}

# Keep the longest alignment for each read.
keep_longest_alignment <- function(df) {
  df <- as.data.frame(df)
  if (!nrow(df)) return(df)
  df <- df[order(df$RID, -(df$end - df$start)), , drop = FALSE]
  df[!duplicated(df$RID), , drop = FALSE]
}

# Load fibertools BED12 records for selected reads.
read_ft_bed12 <- function(bed_file, region = NULL, longest_alignment = FALSE) {
  if (is.null(region)) {
    df <- data.table::fread(bed_file, header = FALSE)
  } else {
    if (!file.exists(paste0(bed_file, ".tbi"))) stop("Missing tabix index for ", bed_file)
    df <- read_tabix_region(bed_file, region)
  }
  if (!nrow(df)) return(as.data.frame(df))
  if (ncol(df) != 12L) stop("Expected fibertools BED12 in ", bed_file)
  names(df) <- c('chr','start','end','RID','score','strand','x1','x2','rgb','blockCount','blockSizes','blockStarts')
  if (longest_alignment) {
    df <- keep_longest_alignment(df)
  }
  return(df)
}

# Expand footprint blocks into individual BED intervals.
convert_ft_bed12_to_bed6 <- function(bed12_df, include_read_start_end = FALSE,
                                   drop_sentinels = TRUE, keep_block_scores = FALSE,
                                   format = NULL, source = NULL,
                                   longest_alignment = FALSE, validate_blocks = FALSE) {
  # an explicit format fixes the expected column count and the sentinel rule
  if (!is.null(format)) {
    expected <- footprint_format_columns(format)
    if (!format %in% c("bed12_fibertools", "bed13_fiberhmm", "bed15_fiberhmm_tf"))
      stop("Not a BED12 block format: ", format)
    configured_sentinels <- format == "bed12_fibertools"
    if (!missing(drop_sentinels) && !identical(drop_sentinels, configured_sentinels))
      stop("drop_sentinels conflicts with configured format '", format, "'")
    drop_sentinels <- configured_sentinels
    if (ncol(bed12_df) != expected)
      footprint_column_error(if (is.null(source)) "<unknown>" else source,
                             format, ncol(bed12_df))
  }
  if (nrow(bed12_df) == 0) {
    return(data.frame())
  }
  # Existing calls retain the BED12-only contract. FiberHMM callers opt in by
  # keeping all blocks or requesting its per-block scores.
  if (is.null(format)) {
    legacy <- isTRUE(drop_sentinels) && !isTRUE(keep_block_scores)
    if (legacy && ncol(bed12_df) != 12L) stop("Expected 12 BED columns")
    if (!ncol(bed12_df) %in% c(12L, 13L)) stop("Expected 12 or 13 BED columns")
  }
  has_block_scores <- if (is.null(format)) ncol(bed12_df) == 13L else format != "bed12_fibertools"
  colnames(bed12_df) <- c('chr','start','end','RID','score','strand','read_start','read_end','rgb','blockCount','blockSizes','blockStarts',
                         if (has_block_scores) 'blockScores',
                         if (identical(format, "bed15_fiberhmm_tf")) c('blockEdgeLeft', 'blockEdgeRight'))
  if (identical(format, "bed12_fibertools")) {
    bed12_df$blockCount <- as.integer(bed12_df$blockCount)
    bed12_df$blockSizes <- as.character(bed12_df$blockSizes)
    bed12_df$blockStarts <- as.character(bed12_df$blockStarts)
  }
  if (longest_alignment) bed12_df <- keep_longest_alignment(bed12_df)

  # expand BED12 into one row per block
  block_sizes_list  <- strsplit(sub(",$", "", bed12_df$blockSizes),  ",", fixed = TRUE)
  block_starts_list <- strsplit(sub(",$", "", bed12_df$blockStarts), ",", fixed = TRUE)
  n_blocks <- lengths(block_sizes_list)
  offsets <- suppressWarnings(as.integer(unlist(block_starts_list)))
  sizes <- suppressWarnings(as.integer(unlist(block_sizes_list)))
  # blockCount, blockSizes and blockStarts must agree and be integers
  if (anyNA(bed12_df$blockCount) || any(n_blocks != bed12_df$blockCount) ||
      any(lengths(block_starts_list) != n_blocks) || anyNA(offsets) || anyNA(sizes))
    stop("Invalid BED12 blocks")
  # blockStarts are offsets from chromStart; repeat each read's fields per block
  starts <- rep(bed12_df$start, n_blocks) + offsets
  bed6_df <- data.frame(
    chr    = rep(bed12_df$chr,    n_blocks),
    start  = starts,
    end    = starts + sizes,
    RID    = rep(bed12_df$RID,    n_blocks),
    score  = rep(bed12_df$score,  n_blocks),
    strand = rep(bed12_df$strand, n_blocks),
    stringsAsFactors = FALSE
  )

  if (validate_blocks && (any(bed6_df$end <= bed6_df$start) ||
      any(bed6_df$start < rep(bed12_df$start, n_blocks))))
    stop("Invalid blocks in ", if (is.null(source)) "<unknown>" else source,
         " (configured format '", if (is.null(format)) "<unspecified>" else format, "')")

  if (keep_block_scores && has_block_scores) {
    block_scores_list <- strsplit(sub(",$", "", bed12_df$blockScores), ",", fixed = TRUE)
    if (any(lengths(block_scores_list) != n_blocks))
      stop("blockScores length does not match blockStarts")
    bed6_df$block_score <- as.numeric(unlist(block_scores_list, use.names = FALSE))
  }

  # Preserve the original RID grouping and ordering for every existing LCL
  # call. Keeping all blocks instead preserves alignment/block order.
  if (drop_sentinels) {
    bed6_df <- dplyr::ungroup(dplyr::slice(dplyr::group_by(bed6_df, RID),
                                         -c(1, dplyr::n())))
  }

  if (include_read_start_end) {
    if (drop_sentinels) {
      read_start_end_df <- dplyr::select(bed12_df, RID, read_start, read_end)
      bed6_df <- dplyr::left_join(bed6_df, read_start_end_df, by = "RID")
    } else {
      bed6_df$read_start <- rep(bed12_df$read_start, n_blocks)
      bed6_df$read_end <- rep(bed12_df$read_end, n_blocks)
    }
  }
  return(as.data.frame(bed6_df))
}

# Average values within centered genomic windows.
smooth_mean_slidewindow <- function(x, pos, n = 20) {
  # Average neighboring values around each position.
  sapply(pos, function(i) mean(x[which(pos >= (i - n / 2) & pos <= (i + n / 2))]))
}

