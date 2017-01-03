\name{reduce2}
\alias{reduce2}
\title{
Merge genomic regions
}
\description{
Merge genomic regions
}
\usage{
reduce2(gr, max_gap = 1000, gap = bp(1000), add_n_column = TRUE)
}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{max_gap}{maximum gap to merge, measured in base pairs. Only work if \code{gap} is the ratio.}
  \item{gap}{a numeric value means to extend each region by \code{gap} times of its width before merging. If \code{gap} represents number of base pairs, use \code{\link{bp}}, \code{\link{kb}} or \code{\link{mb}} to wrap it. If \code{gap} represents absolute number of base pairs, the functionality is same as \code{\link[GenomicRanges]{reduce}}.}
  \item{add_n_column}{whether to add a column which represents number of regions merged, used internally.}

}
\details{
\code{\link[GenomicRanges]{reduce}} only merges regions with fixed gap width, but sometimes it is not reasonable to set gap
to a same width for all regions. Assuming we have a list of differentially methylated regions (DMRs) and we want to reduce
the number of DMRs by merging neighouring DMRs. DMRs distribute differently in different places in the genome, e.g. DMRs are dense
and short in CpG-rich regions (e.g. CpG islands) while long in CpG-sparse regions (e.g. gene bodies and intergenic regions),
thus the merging should be applied based to the width of every DMR itself. \code{\link{reduce2}} can merge regions by the width of every region itself.
This type of merging is dynamic because after each iteration of merging, width and number of regions change which results in changing of extension
of some regions. The whole merging will proceed iteratively unless there is no new merging or the gap between regions reaches \code{max_gap}.

If there are numeric meta columns, corresponding values will be summed up for the merged regions. There will
be a new column \code{.__n__.} added which represents number of regions that are merged.
}
\value{
a \code{\link[GenomicRanges]{GRanges}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 4, 8 ,12), 
    end = c(2, 5, 10, 20)))
reduce2(gr, gap = bp(2))
reduce2(gr, gap = 0.1)
}
