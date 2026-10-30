\name{make_sample_sheet}
\alias{make_sample_sheet}
\alias{validate_inputs}
\alias{example_data}
\title{Prepare and validate explicit-format sample inputs}
\description{Create a portable sample sheet and validate indexed tracks against
selected genomic ranges. Copy a small real example into a writable directory.}
\usage{
make_sample_sheet(samples, output = NULL, base_dir = getwd())
validate_inputs(samples, regions)
example_data(directory = tempfile("fiberseq-example-"))
}
\arguments{
\item{samples}{A data.frame, DataFrame, or TSV filename. Required columns are
sample_name, m6a and m6a_format. Optional columns are sample_label, nucleosome,
nucleosome_format, tf, tf_format, fire and fire_format.}
\item{output}{Optional sample-sheet TSV filename.}
\item{base_dir}{Base for relative paths in an in-memory table. A TSV filename
instead resolves relative paths against its own directory.}
\item{regions}{A GRanges of selected loci.}
\item{directory}{New or empty directory for a writable example copy.}
}
\details{Paths may contain \code{\{sample_name\}} and \code{\{chr\}}
placeholders. Every supplied track needs a bgzip file and adjacent .tbi index.
Explicit formats are m6a: bed12_fibertools; nucleosome: bed12_fibertools or
bed13_fiberhmm; tf: bed15_fiberhmm_tf or bed6_per_sample; fire: bed11_fire.
Optional tracks can be omitted or blank. Supplied missing files are errors.
The input format is never inferred from a filename or species. Empty optional
tracks are accepted; an absent m6A contig is an error.}
\value{\code{make_sample_sheet} returns an S4Vectors DataFrame with absolute
paths. \code{validate_inputs} returns TRUE invisibly, or raises an error.
\code{example_data} returns a list with samples (DataFrame), regions (GRanges)
and directory (the copied example location).}
\examples{
e <- example_data()
e$samples
validate_inputs(e$samples, e$regions)
}
