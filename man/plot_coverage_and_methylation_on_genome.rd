\name{plot_coverage_and_methylation_on_genome}
\alias{plot_coverage_and_methylation_on_genome}
\title{
Plot coverage and methylation for a single sample
}
\description{
Plot coverage and methylation for a single sample
}
\usage{
plot_coverage_and_methylation_on_genome(sid, chromosome = paste0("chr", 1:22),
    species = "hg19", nw = 10000, ...)
}
\arguments{

  \item{sid}{a single sample id}
  \item{chromosome}{a vector of chromosome names}
  \item{species}{species}
  \item{nw}{number of windows}
  \item{...}{pass to \code{\link[gtrellis]{gtrellis_layout}}}

}
\details{
The whole genome is segented by \code{nw} windows and mean methylation and mean CpG coverage
are visualized as two tracks.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
