\name{get_peak_list}
\alias{get_peak_list}
\title{
Get a list of peak regions
}
\description{
Get a list of peak regions
}
\usage{
get_peak_list(mark, sample_id = chipseq_hooks$sample_id(mark))
}
\arguments{

  \item{mark}{mark type}
  \item{sample_id}{a vector of sample ids}

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
