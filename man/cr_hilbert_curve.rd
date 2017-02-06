\name{cr_hilbert_curve}
\alias{cr_hilbert_curve}
\title{
Visualize global correlation by Hilbert curve
}
\description{
Visualize global correlation by Hilbert curve
}
\usage{
cr_hilbert_curve(cr, species = "hg19", chromosome = paste0("chr", 1:22),
    merge_chr = TRUE, add_chr_name = TRUE, title = "cr", legend = lgd, ...)
}
\arguments{

  \item{cr}{correlated regions}
  \item{species}{species}
  \item{chromosome}{chromosomes}
  \item{merge_chr}{whether to merge all chromosomes in one plot}
  \item{add_chr_name}{whether add chromosome names to the plot}
  \item{title}{title of the plot}
  \item{legend}{legend}
  \item{...}{pass to \code{\link[HilbertCurve]{GenomicHilbertCurve}}}

}
\value{
no value is returned
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
