\name{genomic_corr_pintersect}
\alias{genomic_corr_pintersect}
\title{
Intersection between two sets of genomic regions
}
\description{
Intersection between two sets of genomic regions
}
\usage{
genomic_corr_pintersect(query, reference, ...)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{...}{pass to \code{\link{percentOverlaps}}}

}
\details{
For each region in \code{query}, it calculates the percent that is covered by \code{reference}.

The returned value is percent which is how much \code{query} is covered by \code{reference}.

Be careful with the \code{strand} in your GRanges object!
}
\value{
A single correlation value.
}
\seealso{
\code{\link{genomic_regions_correlation}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
genomic_corr_pintersect(gr1, gr2)
}
