\name{common_regions}
\alias{common_regions}
\title{
Find common genomic regions across several samples
}
\description{
Find common genomic regions across several samples
}
\usage{
common_regions(gr_list, min_width = 0, min_coverage = floor(length(gr_list)/4), gap = 0.5)
}
\arguments{

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}}}
  \item{min_width}{minimal width of the common regions}
  \item{min_coverage}{minimal cross-sample coverage for the common regions}
  \item{gap}{gap to merge common regions, pass to \code{\link{reduce2}}}

}
\details{
A common region is defined as a region which covers in at least k samples.
The output can be sent to \code{\link{subgroup_specific_genomic_regions}} to find subgroup
specific regions.
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
