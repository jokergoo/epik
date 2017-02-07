\name{common_regions}
\alias{common_regions}
\title{
Find common genomic regions across samples
}
\description{
Find common genomic regions across samples
}
\usage{
common_regions(gr_list, min_recurrency = floor(length(gr_list)/4),
    gap = bp(1000), max_gap = Inf, min_width = 0)
}
\arguments{

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}}}
  \item{min_recurrency}{minimal cross-sample recurrency for the common regions}
  \item{gap}{gap to merge common regions, pass to \code{\link{reduce2}}}
  \item{max_gap}{maximum gap for merging common regions, pass to \code{\link{reduce2}}}
  \item{min_width}{minimal width for the common regions. It can be used to remove a lot of  very short regions.}

}
\details{
A common region is defined as a region which is recurrent in at least k samples. The process of 
fiding common regions are as follows:

\itemize{
  \item merge regions in all samples into one object
  \item calculate coverage which is the recurrency, removed regions with recurrency less than the cutoff
  \item merge the segments and remove regions which are too short
}

Please note in each sample, regions should not be overlapped.
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
common_regions(gr_list, min_recurrency = 4, gap = bp(1))
common_regions(gr_list, min_recurrency = 4, gap = 0.5)
}
