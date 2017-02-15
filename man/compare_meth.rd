\name{compare_meth}
\alias{compare_meth}
\title{
Compare raw and smoothed methylation
}
\description{
Compare raw and smoothed methylation
}
\usage{
compare_meth(gi, cr_smoothed, cr_raw, txdb = NULL, start = NULL, end = NULL)
}
\arguments{

  \item{gi}{a single gene id}
  \item{cr_smoothed}{correlated regions using smoothed methylation}
  \item{txdb}{transcriptome annotation if \code{start} and \code{end} are not set}
  \item{start}{start position of the region of interested (in the extended gene region)}
  \item{end}{end position of the region of interested (in the extended gene region)}

}
\details{
If \code{start} and \code{end} are not set, the whole extended gene will be plotted.

The aim of this function is see whether smoothing can improve the methylation dataset.

There will be six tracks:

\itemize{
  \item smoothed methylation
  \item correlation between methylation and gene expression
  \item raw methylation
  \item raw methylation for those CpG sites with coverage larger than 25th percential of all CpG coverage
  \item raw methylation for those CpG sites with coverage larger than 50th percential of all CpG coverage
  \item CpG coverage, the bottom, middle and top lines correspond to 25th, 50th and 75th percential from all samples
}
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
