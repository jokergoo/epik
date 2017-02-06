\name{get_chromHMM_list}
\alias{get_chromHMM_list}
\title{
Get a list of chromatin segmentation regions
}
\description{
Get a list of chromatin segmentation regions
}
\usage{
get_chromHMM_list(sample_id, ...)
}
\arguments{

  \item{sample_id}{a vector of sample IDs. If not defined, it is the total samples that are available for this histome mark.}
  \item{...}{more arguments pass to \code{chipseq_hooks$chromHMM()}.}

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
