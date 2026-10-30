\name{prepare_regions}
\alias{prepare_regions}
\alias{extract_canonical_tss}
\title{Prepare genomic regions as GRanges}
\description{Import coordinates or transcript annotations using Bioconductor
ranges. Coordinates inside returned objects are always one-based inclusive.}
\usage{
prepare_regions(regions = NULL, output = NULL, coordinate_system = NULL,
    annotation = NULL, genes = NULL, transcripts = NULL, flank_bp = 1000,
    bed_annotation = c("plain", "gencode_custom"))
extract_canonical_tss(annotation, output = NULL, genes = NULL,
    attribute = "tag", value = "Ensembl_canonical", flank_bp = 0)
}
\arguments{
\item{regions}{GRanges, a BED filename, coordinate strings such as
\code{"Chr5:1001-3000"}, or a data.frame/DataFrame with chr, start and end.}
\item{output}{Optional TSV filename; exported start/end follow BED conventions.}
\item{coordinate_system}{Required for tables: \code{"bed"} (zero-based,
half-open) or \code{"one_based"} (one-based, inclusive).}
\item{annotation}{GTF or GFF3 filename containing transcript or mRNA rows.}
\item{genes}{Optional gene identifiers or gene names to select.}
\item{transcripts}{Optional transcript identifiers to select.}
\item{flank_bp}{Promoter flank in bases. A positive F gives
[TSS - F, TSS + F - 1] in one-based coordinates; zero gives the TSS base.}
\item{bed_annotation}{Use plain BED by default. \code{"gencode_custom"}
accepts a single-base TSS BED whose name is gene ID;annotation;gene name.}
\item{attribute, value}{Annotation attribute and value identifying canonical
transcripts. Missing or ambiguous canonical selections raise an error.}
}
\value{A \link[GenomicRanges:GRanges-class]{GRanges} with unique locus names
and metadata including TSS, gene IDs and comma-separated transcript IDs when
available. Distinct TSSs within a gene are retained; shared TSSs retain all
associated transcript IDs. Contig names and identifiers are preserved exactly.}
\details{Ordinary extraction requires no canonical tags. GFF3 transcripts
must have unique IDs and gene Parents. Promoters extending before base 1 are
rejected. No chromosome-name or species restrictions are imposed. The genome
assembly and contig naming must agree with the processed tracks.}
\examples{
prepare_regions(c("Chr5:1001-3000", "scaffold_42:4001-6000"))
a <- system.file("extdata", "plant", "annotation.gff3",
    package = "selectFiberseqLociWorkflow")
prepare_regions(annotation = a)
extract_canonical_tss(a, genes = "AT5G00010", attribute = "canonical",
    value = "true")
}
