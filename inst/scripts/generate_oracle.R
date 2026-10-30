## Regenerate independent tutorial expectations from unchanged shared sources.
args <- commandArgs(TRUE)
root <- normalizePath(args[1])
library(dplyr)
library(GenomicRanges)
library(Matrix)
source(args[2]) # unchanged shared parsing source identified in SHARED_SOURCES.json
base <- file.path(root, "inst", "extdata", "example")
samples <- read.delim(file.path(base, "samples.tsv"))
samples$fire_dir <- ""
regions <- read.delim(file.path(base, "regions.tsv"))
out <- file.path(root, "inst", "extdata", "expected")
dir.create(out, showWarnings = FALSE)
for (i in seq_len(nrow(regions))) {
    region <- regions[i, ]
    paths <- file.path(base, gsub("{chr}", region$chr, samples$m6a, fixed=TRUE))
    dat <- assemble_region_m6a(samples, region, full_span = TRUE,
        positions = seq(region$analysis_start, region$analysis_end),
        m6a_paths = paths)
    m <- dat$met_mat
    d <- Matrix::summary(m)
    calls <- data.frame(RID = rownames(m)[d$i], pos = as.integer(colnames(m))[d$j],
        base = rep("A", nrow(d)))
    profile <- pileup_reads(calls, dat$rids_df, region$analysis_start,
        region$analysis_end, smooth = TRUE, window_n = 10)
    saveRDS(list(matrix = m, profile = profile),
        file.path(out, paste0(region$region_id, ".rds")), version = 2)
}
