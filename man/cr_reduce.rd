\name{cr_reduce}
\alias{cr_reduce}
\title{
Merge correlated regions
}
\description{
Merge correlated regions
}
\usage{
cr_reduce(cr, txdb, expr = NULL, gap = bp(1), mc.cores = 1)
}
\arguments{

  \item{cr}{correlated regions from \code{\link{correlated_regions}}. In most cases, it is correlated regions with significant correlations.}
  \item{txdb}{the transcriptome annotation which is same as the one used in \code{\link{correlated_regions}}}
  \item{expr}{the expression matrix which is same as the one used in \code{\link{correlated_regions}}. If it is set the correlation will be re-calculated for the merged regions.}
  \item{gap}{gap for the merging, pass to \code{\link{reduce2}}}
  \item{mc.cores}{cores for parallel computing. It is paralleled by gene}

}
\details{
As there are overlaps between two neighbouring correlated regions with the default settings, it is possible to merge them into
large regions. The mering is gene-wise, and all statistics (e.g. mean methylation, correlation) will be
re-calculated.
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
