\name{get_sig_cr}
\alias{get_sig_cr}
\title{
Get significant CRs
}
\description{
Get significant CRs
}
\usage{
get_sig_cr(cr, fdr, meth_diff)
}
\arguments{

  \item{cr}{correlated regions}
  \item{fdr}{cutoff for FDRs}
  \item{meth_diff}{cutoff for methylation difference}

}
\details{
When there is no or only one subgroup, the filtering is \code{cr$corr_fdr <= fdr & cr$meth_IQR >= meth_diff},
while when there is more than one subgroups, the filtering is \code{cr$corr_fdr <= fdr & cr$meth_anova_fdr <= fdr & cr$meth_diameter >= meth_diff}.
}
\value{
Significant CRs. \code{\link{cr_reduce}} can be used to reduce number of regions.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
