\name{genomic_corr_jaccard}
\alias{genomic_corr_jaccard}
\title{
Jaccard coefficient between two sets of genomic regions
}
\description{
Jaccard coefficient between two sets of genomic regions
}
\usage{
genomic_corr_jaccard(query, reference, background = NULL)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{background}{background regions that should be only looked in, a \code{\link[GenomicRanges]{GRanges}} object}

}
\details{
Jaccard coefficient is defined as the total length of intersection divided by total
length of union of two sets of genomic regions.

You can set the background when calculating Jaccard coefficient. For example,
if the interest is the Jaccard coefficient between CpG sites in \code{query} and in \code{reference}
\code{background} can be set with a \code{\link[GenomicRanges]{GRanges}} object which contains positions of CpG sites.

Be careful with the \code{strand} in your \code{\link[GenomicRanges]{GRanges}} object!
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
genomic_corr_jaccard(gr1, gr2)
}
