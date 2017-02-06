\name{sig_cr_compare_cutoff}
\alias{sig_cr_compare_cutoff}
\title{
Compare cutoff for determing significant correlated regions
}
\description{
Compare cutoff for determing significant correlated regions
}
\usage{
sig_cr_compare_cutoff(cr, txdb, fdr_cutoff = c(0.1, 0.05, 0.01),
    meth_diff_cutoff = c(0, 0.1, 0.2, 0.3))
}
\arguments{

  \item{cr}{correlated regions, should be returned by \code{\link{cr_enrichedheatmap}}}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{fdr_cutoff}{a list of cutoffs to compare}
  \item{meth_diff_cutoff}{a list of cutoffs to compare}

}
\details{
It simply plot how correlated signals are enriched at extended gene regions for
negative correlated regions and positive correlated regions under different combination
of cutoffs.
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
