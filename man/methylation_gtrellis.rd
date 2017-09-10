\name{methylation_gtrellis}
\alias{methylation_gtrellis}
\title{
Plot coverage and methylation for a single sample
}
\description{
Plot coverage and methylation for a single sample
}
\usage{
methylation_gtrellis(sid, chromosome = paste0("chr", 1:22),
    species = "hg19", nw = 10000, pch = 16, pt_gp = gpar(size = unit(1, "mm")), transparency = 0.8,
    title = qq("Distribution of CpG coverage and methylation for @{sid}"), ...)
}
\arguments{

  \item{sid}{a single sample ID}
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
\seealso{
\code{\link{methylation_gtrellis_multiple_samples}} visualizes methylation for multiple samples.
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
