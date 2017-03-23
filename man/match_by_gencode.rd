\name{match_by_gencode}
\alias{match_by_gencode}
\title{
Filter one GTF annotation by another
}
\description{
Filter one GTF annotation by another
}
\usage{
match_by_gencode(gtf1, gtf2, filter = NULL)
}
\arguments{

  \item{gtf1}{path for GTF file 1}
  \item{gtf2}{path for GTF file 2}
  \item{filter}{code which additionally filters recodes in \code{gtf1}}

}
\details{
In some senarios, the analysis was done with an old version of Gencode and it is impossible to redo it
with a new version of Gencode (e.g. you dont have access to the original bam files). The only way is to remove those gene/transcript annotations which are not
consistent in the two versions. This \code{\link{match_by_gencode}} function only keeps genes/transcripts that have
same gene ids and positions in the two versions.

Similar as \code{\link{import_gencode_as_txdb}}, \code{filter} argument can be used to e.g. only mathch protein coding genes/transcripts.
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
