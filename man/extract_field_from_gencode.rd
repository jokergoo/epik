\name{extract_field_from_gencode}
\alias{extract_field_from_gencode}
\title{
Extract field from gencode GTF file
}
\description{
Extract field from gencode GTF file
}
\usage{
extract_field_from_gencode(file, level = "gene",
    primary_key = "gene_id", field = "gene_name")
}
\arguments{

  \item{file}{the input GTF file}
  \item{level}{level of the annotation (e.g. gene, transcript, exon, the third column in GTF file)}
  \item{primary_key}{primary field}
  \item{field}{field to be retrieved}

}
\details{
Although GTF file can be imported by e.g. \code{\link[GenomicFeatures]{makeTxDbFromGFF}}, some information
in the original GTF file will not be imported. This function aims to extract additionally information
from GTF file.

The function calls external Perl script, so you need to have Perl installed.
}
\value{
A vector in which 'primary_key' corresponds to the name and 'field' corresponds to the value.
}
\seealso{
\code{\link{available_gencode_fields}} lists all possible values for \code{primary_key} and \code{field}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
download.file("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
    destfile = "gencode.v19.annotation.gtf.gz")
extract_field_from_gencode("gencode.v19.annotation.gtf.gz")
extract_field_from_gencode("gencode.v19.annotation.gtf.gz", field = "gene_type")
}
NULL
}
