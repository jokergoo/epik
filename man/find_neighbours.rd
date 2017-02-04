\name{find_neighbours}
\alias{find_neighbours}
\title{
Find neighbour regions
}
\description{
Find neighbour regions
}
\usage{
find_neighbours(query, reference, upstream = 1000, downstream = 1000)
}
\arguments{

  \item{query}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{a \code{\link[GenomicRanges]{GRanges}} object }
  \item{upstream}{upstream that \code{query} is extended}
  \item{downstream}{downstream that \code{query} is extended}

}
\details{
With a certain extension of \code{query}, this funciton looks for \code{reference} which intersects the extended regions.
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object which contains regions in \code{query} for which the extended regisons are overlapped with \code{reference}.
There are three meta columns added:

\describe{
  \item{distance:}{distance from the \code{query} to corresponding \code{reference}}
  \item{query_index:}{index of regions in \code{query}}
  \item{reference_index:}{index of regions in \code{reference} }
}

Note one \code{reference} can correspond to multiple \code{query} regions.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
find_neighbours(gr1, gr2, upstream = 3, downstream = 3)
find_neighbours(gr1, gr2, upstream = 10, downstream = 10)
}
