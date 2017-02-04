\name{genomic_corr_reldist}
\alias{genomic_corr_reldist}
\title{
Relative distance between two sets of genomic regions
}
\description{
Relative distance between two sets of genomic regions
}
\usage{
genomic_corr_reldist(gr1, gr2)
}
\arguments{

  \item{gr1}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{gr2}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}

}
\details{
For regions in \code{gr1} and \code{gr2}, they are all degenerated as single points
which are the middle points of regions. For each middle point in \code{gr1}, it looks 
for two nearest points in \code{gr2} on its left and right. The statistic is defined as the ratio of the distance
to the nearest neighbour point to the distance of two neighbour points. If \code{gr1} and \code{gr2} are not correlated at all,
It is expected that the ratio follows a uniform distribution. So final statisitic are the KS-statistic
between the real distribution of rations to the uniform distribution.
}
\references{
Favoriv A, et al. Exploring massive, genome scale datasets with the GenometriCorr package. PLoS Comput Biol. 2012 May; 8(5):e1002529
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
genomic_corr_reldist(gr1, gr2)
}
