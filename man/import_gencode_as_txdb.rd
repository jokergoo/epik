\name{import_gencode_as_txdb}
\alias{import_gencode_as_txdb}
\title{
Import gencode GTF file as a TxDb object
}
\description{
Import gencode GTF file as a TxDb object
}
\usage{
import_gencode_as_txdb(gtf, filter = NULL)
}
\arguments{

  \item{gtf}{path of the GTF file}
  \item{filter}{code which filters the GTF records}

}
\details{
For example, you can build a \code{\link{TxDb-class}} object only for protein coding genes by defining

  \preformatted{
  import_gencode_as_txdb(GTF, gene_type == "protein_coding" & transcript_type == "protein_coding")  }

Here \code{gene_type} and \code{transcript_type} are attributes in the GTF file.

Please note, when building the \code{\link{TxDb-class}} object, the positions of genes are calculated by the union of all
its transcripts, while the positions in the GTF file are not used. This is important when only using
a subset of transcripts for a gene (e.g. only use protein coding transcripts) that the position of the 
gene may change. So, when number of transcripts for genes change, the corresponding \code{\link{TxDb-class}} object must be
re-generated accordingly.
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
