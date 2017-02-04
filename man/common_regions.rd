\name{common_regions}
\alias{common_regions}
\title{
Find common genomic regions across samples
}
\description{
Find common genomic regions across samples
}
\usage{
common_regions(gr_list, min_coverage = floor(length(gr_list)/4),
    gap = bp(1000), max_gap = Inf, min_width = 0)
}
\arguments{

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}}}
  \item{min_coverage}{minimal cross-sample coverage for the common regions}
  \item{gap}{gap to merge common regions, pass to \code{\link{reduce2}}}
  \item{max_gap}{maximum gap for merging common regions, pass to \code{\link{reduce2}}}
  \item{min_width}{minimal width for the common regions. It can be used to remove a lot of  very short regions.}

}
\details{
First a list of segments are filtered by the cross-sample coverage (\code{min_coverage}). Then close segments
are merged by \code{max_gap}.
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object contains coordinates of common regions. The columns in meta data
are percent of the common region which is covered by regions in every sample.
}
\seealso{
The returned variable can be sent to \code{\link{subgroup_specific_genomic_regions}} to find subgroup specific regions.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
gr_list = list(
    gr1 = GRanges(seqnames = "chr1", ranges = IRanges(1, 8)),
    gr2 = GRanges(seqnames = "chr1", ranges = IRanges(3, 9)),
    gr3 = GRanges(seqnames = "chr1", ranges = IRanges(2, 7)),
    gr4 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 6), c(4, 10))),
    gr5 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 9), c(3, 10)))
)
common_regions(gr_list, min_coverage = 4, gap = bp(1))
common_regions(gr_list, min_coverage = 4, gap = 0.5)
}
