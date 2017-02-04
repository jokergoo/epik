\name{sig_cr_enrichedheatmap}
\alias{sig_cr_enrichedheatmap}
\title{
Visualize significant correlated regions
}
\description{
Visualize significant correlated regions
}
\usage{
sig_cr_enrichedheatmap(cr, txdb, fdr_cutoff = 0.05, meth_diff_cutoff = 0.1)
}
\arguments{

  \item{cr}{correlated regions, should be returned by \code{\link{cr_enrichedheatmap}}.}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{fdr_cutoff}{cutoff for fdr}
  \item{meth_diff_cutoff}{cutoff for methylation difference}

}
\details{
There are two heatmaps which corresponds to negative correlated regions and positive
correlated regions.
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
