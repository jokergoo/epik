\name{available_gencode_fields}
\alias{available_gencode_fields}
\title{
Returns all supported fields in GTF file
}
\description{
Returns all supported fields in GTF file
}
\usage{
available_gencode_fields(file, level = "gene")
}
\arguments{

  \item{file}{the input GTF file}
  \item{level}{level of the annotation (e.g. gene, transcript, exon, the third column in GTF file)}

}
\details{
These fields are stored in the 9th column in the GTF file.
}
\value{
A vector of available fields.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
download.file("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
    destfile = "gencode.v19.annotation.gtf.gz")
available_gencode_fields("gencode.v19.annotation.gtf.gz", level = "gene")
available_gencode_fields("gencode.v19.annotation.gtf.gz", level = "transcript")
}
NULL
}
