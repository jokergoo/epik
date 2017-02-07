\name{sig_cr_enriched_heatmap}
\alias{sig_cr_enriched_heatmap}
\title{
Visualize significant correlated regions
}
\description{
Visualize significant correlated regions
}
\usage{
sig_cr_enriched_heatmap(cr, txdb, fdr_cutoff = 0.05, meth_diff_cutoff = 0.1)
}
\arguments{

  \item{cr}{correlated regions, should be returned by \code{\link{cr_enriched_heatmap}}.}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{fdr_cutoff}{cutoff for fdr}
  \item{meth_diff_cutoff}{cutoff for methylation difference. If there are no subgroup information or only one subgroup, \code{meth_IQR} column is used for filtering. If there are more than one subgroups, \code{meth_diameter} column is used for filtering.}

}
\details{
There are two heatmaps which corresponds to significant negative correlated regions and positive
correlated regions. Rows are same as the heatmaps produced by \code{\link{cr_enriched_heatmap}}.
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
