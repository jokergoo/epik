\name{correlated_regions_by_window}
\alias{correlated_regions_by_window}
\title{
Correlated regions in a specified region
}
\description{
Correlated regions in a specified region
}
\usage{
correlated_regions_by_window(site, meth, expr, chr, cov = NULL, cov_cutoff = 3, min_dp = 4,
    cor_method = "spearman", window_size = 5, window_step = window_size, subgroup = NULL, max_width = 10000)
}
\arguments{

  \item{site}{position of CpG sites in this region, should be sorted}
  \item{meth}{methylation matrix corresponding to \code{site}}
  \item{expr}{expression for the associated gene}
  \item{chr}{chromosome name, used to construct the \code{\link[GenomicRanges]{GRanges}} object}
  \item{cov}{CpG coverage matrix. CpG coverage is important when \code{meth} is the raw methylation which means CpG sites with extremely low coverage will be removed when calculating correlations}
  \item{cov_cutoff}{cutoff for the CpG coverage}
  \item{min_dp}{minimal number of non-NA values for calculating correlations. When \code{meth} is the raw methylation, values for which CpG coverage is too low will be replaced with \code{NA}, We only use non-NA values to calculate correlations.}
  \item{cor_method}{method for calcualting correlations}
  \item{window_size}{how many CpG sites in a window}
  \item{window_step}{step of the sliding window, measured in number of CpG sites}
  \item{subgroup}{subgroup information}
  \item{max_width}{maximum width of a window}

}
\details{
\code{cov} and \code{cov_cutoff} should be set when the methylation is unsmoothed, because
for the unsmoothed data, the methylation rate is not unreliable when the CpG coverage is low.
}
\value{
a \code{\link[GenomicRanges]{GRanges}} object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
