\name{plot_multiple_samples_methylation_on_genome}
\alias{plot_multiple_samples_methylation_on_genome}
\title{
Plot methylation for multiple samples
}
\description{
Plot methylation for multiple samples
}
\usage{
plot_multiple_samples_methylation_on_genome(sample_id, annotation,
    chromosome = paste0("chr", 1:22), species = "hg19", nw = 1000, ...)
}
\arguments{

  \item{sample_id}{a vector of sample ids}
  \item{annotation}{annotation of samples (e.g. subtypes)}
  \item{chromosome}{a vector of chromosome names}
  \item{species}{species}
  \item{nw}{number of windows}
  \item{...}{pass to \code{\link[gtrellis]{gtrellis_layout}}}

}
\details{
The whole genome is segented by \code{nw} windows. Methylation in different classes are visualized as separated tracks.
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
