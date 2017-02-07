\name{sig_cr_compare_cutoff}
\alias{sig_cr_compare_cutoff}
\title{
Compare cutoff for determining significant correlated regions
}
\description{
Compare cutoff for determining significant correlated regions
}
\usage{
sig_cr_compare_cutoff(cr, txdb, fdr_cutoff = c(0.1, 0.05, 0.01),
    meth_diff_cutoff = c(0, 0.1, 0.2, 0.3))
}
\arguments{

  \item{cr}{correlated regions, should be returned by \code{\link{cr_enriched_heatmap}}}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{fdr_cutoff}{a list of cutoffs to compare}
  \item{meth_diff_cutoff}{a list of cutoffs to compare. If there are no subgroup information or only one subgroup, \code{meth_IQR} column is used for filtering. If there are more than one subgroups, \code{meth_diameter} column is used for filtering.}

}
\details{
It simply plots how correlated signals are enriched at extended gene regions for
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
