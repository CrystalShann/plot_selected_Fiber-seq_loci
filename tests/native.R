library(selectFiberseqLociWorkflow)
library(SummarizedExperiment)
library(GenomicRanges)
fail <- function(expr, pattern = NULL) {
    err <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(err, "error"))
    if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(err)))
}
collect_warnings <- function(expr) {
    messages <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = value, messages = messages)
}
check_alignment <- function(x) {
    stopifnot(is(x, "RangedSummarizedExperiment"), validObject(x),
        is(assay(x), "sparseMatrix"), is(rowRanges(x), "GRanges"),
        is(colData(x)$read_range, "GRanges"),
        is(colData(x)$footprints, "GRangesList"),
        identical(rownames(x), names(rowRanges(x))),
        identical(colnames(x), rownames(colData(x))),
        identical(colnames(x), colData(x)$RID),
        identical(start(colData(x)$read_range), as.integer(colData(x)$start)),
        identical(end(colData(x)$read_range), as.integer(colData(x)$end)))
    footprints <- colData(x)$footprints
    for (j in seq_along(footprints))
        stopifnot(all(mcols(footprints[[j]])$RID == colnames(x)[j]))
}

example <- example_data()
stopifnot(is(example$samples, "DataFrame"), is(example$regions, "GRanges"))
fail(analyze_loci(example$samples, example$regions), "clustering")
a <- analyze_loci(example$samples, example$regions, clustering = FALSE)
set.seed(91)
rng_before <- .Random.seed
b <- analyze_loci(example$samples, example$regions, clustering = TRUE)
stopifnot(identical(rng_before, .Random.seed))
for (i in seq_along(a)) {
    x <- a[[i]]; y <- b[[i]]
    check_alignment(x); check_alignment(y)
    stopifnot(identical(assay(x), assay(y)), identical(colnames(x), colnames(y)),
        identical(locus_profile(x), locus_profile(y)),
        !"cluster" %in% names(colData(x)), "cluster" %in% names(colData(y)))
    sub <- x[c(20, 3, 100), c(3, 1)]
    check_alignment(sub)
    stopifnot(identical(assay(sub), assay(x)[c(20, 3, 100), c(3, 1)]),
        identical(colData(sub)$footprints, colData(x)$footprints[c(3, 1)]))
    expected <- readRDS(system.file("extdata", "expected",
        paste0(names(a)[i], ".rds"), package = "selectFiberseqLociWorkflow"))
    m <- expected$matrix[colnames(x), , drop = FALSE]
    stopifnot(identical(unname(as.matrix(assay(x))), unname(t(as.matrix(m)))))
    profile <- as.data.frame(locus_profile(x))
    stopifnot(isTRUE(all.equal(profile, expected$profile, check.attributes = FALSE)),
        all(profile$cov == ncol(x)))
}
# Exports retain native containers and identical unclustered profiles.
stopifnot(identical(cluster_colors(c("cluster3", "cluster1")),
    c(cluster3 = "green4", cluster1 = "dodgerblue2")),
    identical(cluster_colors(factor("cluster3", levels = paste0("cluster", 1:3))),
        c(cluster1 = "dodgerblue2", cluster2 = "#E31A1C", cluster3 = "green4")),
    length(cluster_colors(character())) == 0L,
    length(cluster_colors(paste0("cluster", 1:30))) == 30L)
fail(cluster_colors(c("cluster1", NA)), "nonmissing")
for (clustered in c(FALSE, TRUE)) {
    h <- plot_locus(b[[1]], "heatmap", clustered = clustered)
    stopifnot(identical(names(h@left_annotation@anno_list),
        if (clustered) c("sample", "cluster") else "sample"),
        identical(names(h@bottom_annotation@anno_list), c("coordinate", "TSS")))
    composition <- plot_locus(b[[1]], "footprints", clustered = clustered)
    stopifnot(inherits(composition, "patchwork"), length(composition) == 3L)
    expected_order <- with(as.data.frame(colData(b[[1]])),
        if (clustered) order(cluster, match(sample_name, example$samples$sample_name),
            start, RID) else order(match(sample_name, example$samples$sample_name),
            start, RID))
    for (panel in c(2L, 3L)) {
        layers <- composition[[panel]]$layers
        bars <- layers[[which(vapply(layers, function(z)
            "annotation_key" %in% names(z$data), logical(1)))]]$data
        stopifnot(identical(levels(bars$RID), rev(colnames(b[[1]])[expected_order])))
        if (clustered) stopifnot(
            max(bars$annotation_x[grepl("^Sample:", bars$annotation_key)]) <
            min(bars$annotation_x[grepl("^Cluster:", bars$annotation_key)]))
        # Building also exercises facet/scale alignment with read metadata.
        ggplot2::ggplot_build(composition[[panel]])
    }
}
out <- tempfile("native-export-")
write_locus_results(b, out)
stopifnot(file.exists(file.path(out, "ALAS1", "clustered", "m6a_heatmap.pdf")))
check_alignment(readRDS(file.path(out, "ALAS1", "locus.rds")))

# Declared coordinate conversions, preserved identifiers and seqinfo.
r1 <- prepare_regions("Chr5:1-100")
r2 <- prepare_regions(data.frame(chr = "Chr5", start = 0, end = 100),
    coordinate_system = "bed")
stopifnot(identical(start(r1), start(r2)), identical(end(r1), end(r2)))
fail(prepare_regions(data.frame(chr = "Chr5", start = 0, end = 100)))
si <- GenomeInfoDb::Seqinfo("Chr5", seqlengths = 10000, genome = "fixture")
gr <- GRanges("Chr5", IRanges(1001, 3000), seqinfo = si)
names(gr) <- "plant:gene+1"
stopifnot(identical(seqinfo(prepare_regions(gr)), si),
    identical(names(prepare_regions(gr)), names(gr)))
plant <- system.file("extdata", "plant", package = "selectFiberseqLociWorkflow")
pg <- prepare_regions(annotation = file.path(plant, "annotation.gff3"))
stopifnot(length(pg) == 3L,
    identical(as.character(seqnames(pg)), c("Chr5", "Chr5", "scaffold_42")),
    identical(as.character(strand(pg)), c("+", "+", "-")),
    mcols(pg)$tss[3] == 5000,
    mcols(pg)$transcript_ids[1] == "AT5G00010.1,AT5G00010.3")
fail(extract_canonical_tss(file.path(plant, "annotation.gff3")), "missing/ambiguous")
canonical <- extract_canonical_tss(file.path(plant, "annotation.gff3"),
    genes = "AT5G00010", attribute = "canonical", value = "true")
stopifnot(start(canonical) == 2001L, width(canonical) == 1L)
gtf <- extract_canonical_tss(file.path(plant, "annotation.gtf"))
stopifnot(length(gtf) == 2L, all(width(gtf) == 1L))
ps <- make_sample_sheet(file.path(plant, "samples.tsv"))
p <- collect_warnings(analyze_loci(ps, pg, clustering = TRUE))$value
check_alignment(p[[1]])
stopifnot(ncol(p[[1]]) == 3L, any(Matrix::colSums(assay(p[[1]])) == 0))
fp <- unlist(colData(p[[1]])$footprints, use.names = FALSE)
stopifnot(all(width(fp[mcols(fp)$track == "tf"]) == 20L))
zero <- collect_warnings(analyze_loci(ps,
    prepare_regions("Chr5:5001-5100"), clustering = FALSE))$value[[1]]
stopifnot(ncol(zero) == 3L, sum(assay(zero)) == 0,
    all(locus_profile(zero)$frac == 0))
empty <- collect_warnings(analyze_loci(ps,
    prepare_regions("Chr5:11001-11100"), clustering = TRUE))$value[[1]]
stopifnot(ncol(empty) == 0L, nrow(empty) == 100L, validObject(empty),
    all(is.na(locus_profile(empty)$frac)))
write_locus_results(empty, tempfile("empty-export-"))
filtered <- collect_warnings(analyze_loci(ps, pg[1], clustering = FALSE,
    nuc_min = 1000, nuc_max = 2000, min_tf_score = 100))
stopifnot("no footprints called in the selected region" %in% filtered$messages,
    sum(lengths(colData(filtered$value[[1]])$footprints)) == 0L)
# Optional tracks omitted versus supplied but missing or wrongly declared.
naked <- analyze_loci(ps[, c("sample_name", "m6a", "m6a_format")],
    pg[1], clustering = FALSE)[[1]]
stopifnot(sum(lengths(colData(naked)$footprints)) == 0L)
# Empty optional tracks and zero-call reads still render read backbones.
for (result in list(naked, zero, filtered$value[[1]])) {
    plot <- plot_locus(result, "footprints")
    for (panel in seq_len(3L)) ggplot2::ggplot_build(plot[[panel]])
}
bad <- ps; bad$tf_format <- "bed6_per_sample"
fail(validate_inputs(bad, pg), "Column count")
bad <- ps; bad$tf <- tempfile("missing-track-")
fail(validate_inputs(bad, pg), "Missing bgzip")
fail(validate_inputs(ps, prepare_regions("chr5:1001-3000")), "match exactly")
cat("Native interoperability, tutorial equivalence and edge-case checks passed\n")
