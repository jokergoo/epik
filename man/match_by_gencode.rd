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

  \item{gtf1}{path for gtf1}
  \item{gtf2}{path for gtf2}
  \item{filter}{code which filter the gtf1 records}

}
\details{
In some senarios, the analysis is done with an old version of Gencode and it is impossible to redo it
with a new version of Gencode. The only way is to remove those gene/transcript annotations which are not
consistent in the two version. This \code{\link{match_by_gencode}} function only keeps genes/transcripts that have
same gene ids and positions in the two versions.
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
