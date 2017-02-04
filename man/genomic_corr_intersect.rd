\name{genomic_corr_intersect}
\alias{genomic_corr_intersect}
\title{
Intersections between two sets of genomic regions
}
\description{
Intersections between two sets of genomic regions
}
\usage{
genomic_corr_intersect(gr1, gr2, method = c("number", "percent", "length"), ...)
}
\arguments{

  \item{gr1}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{gr2}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{method}{how to calculate the intersection statistic, see "details"}
  \item{...}{pass to \code{\link[GenomicRanges]{countOverlaps}} or \code{\link{percentOverlaps}}}

}
\details{
There are three metrics for the intersection statistic:

\describe{
  \item{number}{It calculates number of regions in \code{gr1} that overlap with \code{gr2}. Please note this value is not equal to the number of intersections betweenn two sets of regions, because one region in \code{gr1} may overlap with more than one regions in \code{gr2}.}
  \item{percent}{It calculates for each region in \code{gr1}, how much it is covered by regions in \code{gr2}.}
  \item{length}{sum of length of the intersection of the two sets of regions.}
}

With methods of"number" and "percent", \code{genomic_corr_intersect(gr1, gr2)} is always not identical
to \code{genomic_corr_intersect(gr2, gr1)}.
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
genomic_corr_intersect(gr1, gr2, method = "number")
genomic_corr_intersect(gr1, gr2, method = "percent")
genomic_corr_intersect(gr1, gr2, method = "length")
}
