\name{genomic_corr_sintersect}
\alias{genomic_corr_sintersect}
\title{
Intersection between two sets of genomic regions
}
\description{
Intersection between two sets of genomic regions
}
\usage{
genomic_corr_sintersect(query, reference, background = NULL)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{background}{subset of regions that should be only looked into, a \code{\link[GenomicRanges]{GRanges}} object}

}
\details{
It calculates the total length of overlapped regions in \code{query}.

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
genomic_corr_sintersect(gr1, gr2)
}
