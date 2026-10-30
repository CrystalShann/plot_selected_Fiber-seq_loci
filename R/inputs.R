.formats <- list(m6a = "bed12_fibertools", nucleosome = c("bed12_fibertools", "bed13_fiberhmm"),
                 tf = c("bed15_fiberhmm_tf", "bed6_per_sample"), fire = "bed11_fire")
.absolute <- function(p, base = getwd()) {
  p <- path.expand(p)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", p)) return(p)
  file.path(normalizePath(base, mustWork = TRUE), p)
}
.expand_path <- function(p, sample, chr) {
  p <- gsub("{sample_name}", sample, p, fixed = TRUE)
  p <- gsub("{chr}", chr, p, fixed = TRUE)
  if (grepl("[{}]", p)) stop("Only {sample_name} and {chr} path placeholders are supported")
  p
}
make_sample_sheet <- function(samples, output = NULL, base_dir = getwd()) {
  if (is.character(samples) && length(samples) == 1L) {
    base_dir <- dirname(normalizePath(samples, mustWork = TRUE))
    samples <- .read_table(samples)
  }
  if (methods::is(samples, "DataFrame")) samples <- as.data.frame(samples)
  if (!is.data.frame(samples) || !nrow(samples) ||
      !all(c("sample_name", "m6a", "m6a_format") %in% names(samples)))
    stop("Samples require sample_name, m6a, and m6a_format")
  if (anyNA(samples$sample_name) || any(!grepl("^[A-Za-z0-9_.-]+$", samples$sample_name)) || anyDuplicated(samples$sample_name))
    stop("Sample names must be unique nonempty safe names")
  if (!"sample_label" %in% names(samples)) samples$sample_label <- samples$sample_name
  if (anyNA(samples$sample_label) || any(!nzchar(samples$sample_label))) stop("Sample labels must be nonempty")
  for (kind in names(.formats)) {
    if (!kind %in% names(samples)) next
    fmt <- paste0(kind, "_format")
    used <- !is.na(samples[[kind]]) & nzchar(samples[[kind]])
    if (kind == "m6a" && !all(used)) stop("Every sample requires m6a")
    if (any(used) && (!fmt %in% names(samples) || any(!samples[[fmt]][used] %in% .formats[[kind]])))
      stop("Unsupported/missing ", fmt, "; choose ", paste(.formats[[kind]], collapse = ", "))
    samples[[kind]][used] <- vapply(samples[[kind]][used], .absolute, character(1), base = base_dir)
  }
  .write_table(samples, output)
  S4Vectors::DataFrame(samples)
}

validate_inputs <- function(samples, regions) {
    samples <- make_sample_sheet(samples)
    regions <- prepare_regions(regions)
    expected <- c(bed12_fibertools = 12L, bed13_fiberhmm = 13L,
        bed15_fiberhmm_tf = 15L, bed6_per_sample = 6L, bed11_fire = 11L)
    checked <- character()
    for (i in seq_len(nrow(samples))) {
        for (chr in unique(as.character(GenomicRanges::seqnames(regions)))) {
            for (kind in intersect(names(.formats), names(samples))) {
                p <- samples[[kind]][i]
                if (is.na(p) || !nzchar(p)) next
                p <- .expand_path(p, samples$sample_name[i], chr)
                if (!file.exists(p) || !file.exists(paste0(p, ".tbi")))
                    stop("Missing bgzip input or .tbi: ", p)
                seqs <- Rsamtools::seqnamesTabix(p)
                ## Optional tracks may legitimately have no calls on a contig.
                if (kind == "m6a" && !chr %in% seqs)
                    stop("Contig ", chr, " absent from ", p,
                        "; contig names must match exactly")
                fmt <- samples[[paste0(kind, "_format")]][i]
                key <- paste(p, fmt)
                if (key %in% checked) next
                .validate_track_header(p, fmt, expected[[fmt]], kind)
                checked <- c(checked, key)
            }
        }
    }
    invisible(TRUE)
}

.validate_track_header <- function(path, format, expected, kind) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con))
    lines <- readLines(con, n = 100L, warn = FALSE)
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (!length(lines)) {
        if (kind == "m6a") stop("Empty m6A input: ", path)
        return(invisible(NULL))
    }
    if (any(lengths(strsplit(lines, "\t", fixed = TRUE)) != expected))
        stop("Column count does not match ", format, ": ", path)
}

example_data <- function(directory = tempfile("fiberseq-example-")) {
    if (dir.exists(directory) &&
        length(list.files(directory, all.files = TRUE, no.. = TRUE)))
        stop("Choose a new or empty example directory")
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    root <- system.file("extdata", "example", package = "selectFiberseqLociWorkflow")
    files <- list.files(root, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    if (!length(files) || !all(file.copy(files, directory, recursive = TRUE)))
        stop("Could not copy example data")
    directory <- normalizePath(directory)
    list(samples = make_sample_sheet(file.path(directory, "samples.tsv")),
        regions = prepare_regions(.read_table(file.path(directory, "regions.tsv")),
            coordinate_system = "bed"), directory = directory)
}
