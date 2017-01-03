\name{genomic_corr_nintersect}
\alias{genomic_corr_nintersect}
\title{
Intersections between two sets of genomic regions
}
\description{
Intersections between two sets of genomic regions
}
\usage{
genomic_corr_nintersect(query, reference, ...)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{...}{pass to \code{\link[GenomicRanges]{countOverlaps}}}

}
\details{
It calculates number of regions in \code{query} that overlap with \code{reference}.

Please note this value is not equal to the number of intersections betweenn two sets of regions,
because one region in \code{query} may overlap with more than one
regions in \code{reference}

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
genomic_corr_nintersect(gr1, gr2)
}
