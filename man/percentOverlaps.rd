\name{percentOverlaps}
\alias{percentOverlaps}
\title{
Overlap genomic regions
}
\description{
Overlap genomic regions
}
\usage{
percentOverlaps(query, subject, ignore_strand = TRUE, ...)
}
\arguments{

  \item{query}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{subject}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{ignore_strand}{wether ignore strands}
  \item{...}{pass to \code{\link[GenomicRanges]{findOverlaps}}}

}
\details{
For every interval in \code{query}, it calculates the percent that is covered by \code{subject}.

Be careful with \code{strand} in your \code{\link[GenomicRanges]{GRanges}} object!
}
\value{
A numeric vector which has the same length as \code{query}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
percentOverlaps(gr1, gr2)
percentOverlaps(gr2, gr1)
}
