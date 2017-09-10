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

  \item{sample_id}{a vector of sample IDs.}
  \item{...}{more arguments pass to \code{\link{chipseq_hooks}}$chromHMM().}

}
\details{
It works after \code{\link{chipseq_hooks}} is set.
}
\value{
A list of \code{\link[GenomicRanges]{GRanges}} objects.

If you e.g. set "chr" as the third argument when defining \code{\link{chipseq_hooks}}$peak(), "chr" can also be passed here through \code{...}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
