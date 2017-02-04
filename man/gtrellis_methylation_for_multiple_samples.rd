\name{gtrellis_methylation_for_multiple_samples}
\alias{gtrellis_methylation_for_multiple_samples}
\title{
Plot methylation for multiple samples as heatmaps
}
\description{
Plot methylation for multiple samples as heatmaps
}
\usage{
gtrellis_methylation_for_multiple_samples(sample_id, subgroup,
    chromosome = paste0("chr", 1:22), species = "hg19", nw = 1000,
    title = qq("genome-wide methylation for @{length(sample_id)} samples"), ...)
}
\arguments{

  \item{sample_id}{a vector of sample ids}
  \item{subtype}{annotation of samples (e.g. subtypes)}
  \item{chromosome}{a vector of chromosome names}
  \item{species}{species}
  \item{nw}{number of windows}
  \item{title}{title of the plot}
  \item{...}{pass to \code{\link[gtrellis]{gtrellis_layout}}}

}
\details{
The whole genome is segented by \code{nw} windows. Methylation in different subgroups are visualized as separated tracks.
}
\value{
No value is returned
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
