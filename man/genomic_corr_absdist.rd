\name{genomic_corr_absdist}
\alias{genomic_corr_absdist}
\title{
Absolute distance between two sets of genomic regions
}
\description{
Absolute distance between two sets of genomic regions
}
\usage{
genomic_corr_absdist(query, reference, method = mean, ...)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{method}{function in which input is a vector of distance and output is a scalar}
  \item{...}{pass to \code{method}}

}
\details{
For regions in \code{query} and \code{reference}, they are all degenerated as single points
which are the middle points of corresponding regions. For each middle point in \code{query}, it looks 
for the nearest point in \code{reference}. Assuming the distance vector is \code{d}, the final statistic is \code{method(d)}.
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
genomic_corr_absdist(gr1, gr2)
}
