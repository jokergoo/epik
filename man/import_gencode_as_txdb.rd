\name{import_gencode_as_txdb}
\alias{import_gencode_as_txdb}
\title{
Import gencode
}
\description{
Import gencode
}
\usage{
import_gencode_as_txdb(gtf, filter = NULL)
}
\arguments{

  \item{gtf}{path of the gtf file}
  \item{filter}{code which filter the gtf records}

}
\details{
For example, you can build a \code{\link{TxDb-class}} object only for protein coding genes by defining
\code{filter = gene_type == "protein_coding" & transcript_type == "protein_coding"}
}
\value{
a \code{\link{TxDb-class}} object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
