\name{gtrellis_coverage_and_methylation}
\alias{gtrellis_coverage_and_methylation}
\title{
Plot coverage and methylation for a single sample
}
\description{
Plot coverage and methylation for a single sample
}
\usage{
gtrellis_coverage_and_methylation(sid, chromosome = paste0("chr", 1:22),
    species = "hg19", nw = 10000, pch = 16, pt_gp = gpar(size = unit(1, "mm")), transparency = 0.8,
    title = qq("Distribution of CpG coverage and methylation for @{sid}"), ...)
}
\arguments{

  \item{sid}{a single sample id}
  \item{chromosome}{a vector of chromosome names}
  \item{species}{species}
  \item{nw}{number of windows to segment the genome}
  \item{pch}{point type}
  \item{pt_gp}{graphic parameters for points (\code{col} will be excluded)}
  \item{transparency}{transparency of points}
  \item{title}{title of the plot}
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
