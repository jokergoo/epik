\name{get_peak_list}
\alias{get_peak_list}
\title{
Get a list of peak regions for a given histone mark
}
\description{
Get a list of peak regions for a given histone mark
}
\usage{
get_peak_list(mark, sample_id = chipseq_hooks$sample_id(mark), ...)
}
\arguments{

  \item{mark}{name of the histone mark}
  \item{sample_id}{a vector of sample IDs. If not defined, it is the total samples that are available for this histone mark.}
  \item{...}{more arguments pass to \code{chipseq_hooks$peak()}.}

}
\details{
It works after \code{\link{chipseq_hooks}} is set.
}
\value{
A list of \code{\link[GenomicRanges]{GRanges}} objects.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
