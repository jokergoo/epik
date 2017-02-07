\name{enrich_with_histone_mark}
\alias{enrich_with_histone_mark}
\title{
Normalize histone modification signals to target
}
\description{
Normalize histone modification signals to target
}
\usage{
enrich_with_histone_mark(target, mark, sample_id, mode = mean, return_arr = FALSE, ...)
}
\arguments{

  \item{target}{target regions}
  \item{mark}{histone mark name}
  \item{sample_id}{a vector of sample ids}
  \item{mode}{how to summarize histone modification signals among samples, by defualt is the cross-sample mean signal}
  \item{return_arr}{whether also return the three dimension array itself}
  \item{...}{pass to \code{\link[EnrichedHeatmap]{normalizeToMatrix}}}

}
\details{
For each sample, the signal is normalized as a matrix, which results an array in which the third
dimension corresponds to samples. The final normalized matrix which shows e.g. mean signal matrix is calcualted
by \code{apply(array, c(1, 2), mode)}.
}
\value{
If \code{return_arr} is set to \code{FALSE}, the funtion returns a matrix which can be directly sent to 
\code{\link[EnrichedHeatmap]{EnrichedHeatmap}}. If \code{return_arr} is \code{TRUE}, the returned value is a list in which
the first element is the original array that each slice in the third dimension is the normalize matrix
in each sample.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
