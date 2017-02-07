\name{reduce2}
\alias{reduce2}
\title{
Merge genomic regions
}
\description{
Merge genomic regions
}
\usage{
reduce2(gr, gap = 0.1, max_gap = Inf, .message = TRUE, .revmap = NULL, ...)
}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{gap}{a numeric value means to extend each region by \code{gap/2} times of its width before merging. If \code{gap} represents number of base pairs, use \code{\link{bp}}, \code{\link{kb}} or \code{\link{mb}} to wrap it. If \code{gap} represents absolute number of base pairs, the functionality is same as \code{\link[GenomicRanges]{reduce}} (\code{gap} is sent to \code{min.gapwidth}).}
  \item{max_gap}{maximum distance to merge, measured in base pairs. Only work if \code{gap} is a ratio value.}
  \item{.message}{internal use}
  \item{.revmap}{internal use}
  \item{...}{further arguments passed to \code{\link[GenomicRanges]{reduce}} (exclude \code{min.gapwidth} and \code{with.revmap})}

}
\details{
\code{\link[GenomicRanges]{reduce}} only merges regions with fixed gap width, but sometimes it is not reasonable to set gap
to a same width for all regions. Assuming we have a list of differentially methylated regions (DMRs) and we want to reduce
the number of DMRs by merging neighouring DMRs. DMRs distribute differently in different places in the genome, e.g. DMRs are dense
and short in CpG-rich regions (e.g. CpG islands) while long in CpG-poor regions (e.g. gene bodies and intergenic regions),
thus the merging should be applied based to the width of every DMR itself. \code{\link{reduce2}} can merge regions by the width of every region itself.
This type of merging is dynamic because after each iteration of merging, some regions are merged into a large region and 
it will has longer extension. The whole merging will proceed iteratively unless there is no new merging.

Note \code{with.revmap} is always set to \code{TRUE} when calling \code{\link[GenomicRanges]{reduce}}, thus there is always a \code{revmap}
meta column in the returned \code{\link[GenomicRanges]{GRanges}} object.
}
\value{
a \code{\link[GenomicRanges]{GRanges}} object with a \code{revmap} meta column.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 4, 8 ,16), 
    end = c(2, 5, 10, 30)), value = 1:4)
reduce2(gr, gap = bp(2))
reduce2(gr, gap = 0.6)
reduce2(gr, gap = 0.6, max_gap = 4)
}
