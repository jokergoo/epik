\name{enrich_with_methylation}
\alias{enrich_with_methylation}
\title{
Normalize methylation to target
}
\description{
Normalize methylation to target
}
\usage{
enrich_with_methylation(target, sample_id, mode = rowMeans, ...)
}
\arguments{

  \item{target}{target regions}
  \item{sample_id}{a vector of sample ids}
  \item{mode}{how to summarize methylation among samples, by default it is the cross-sample mean methylation.  Since methylation is represented as matrix, here we use \code{row*}-family functions (e.g. \code{\link{rowMeans}}, \code{\link[matrixStats]{rowMedians}})}
  \item{...}{pass to \code{\link[EnrichedHeatmap]{normalizeToMatrix}}}

}
\value{
A matrix which can be directly visualized by \code{\link[EnrichedHeatmap]{EnrichedHeatmap}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
