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
  \item{mode}{how to summarize methylation among samples, by defualt is the cross-sample mean methylation}
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
