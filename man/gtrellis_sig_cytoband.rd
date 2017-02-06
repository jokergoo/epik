\name{gtrellis_sig_cytoband}
\alias{gtrellis_sig_cytoband}
\title{
Visualize correlations in cytoband
}
\description{
Visualize correlations in cytoband
}
\usage{
gtrellis_sig_cytoband(cr, txdb, cytoband_list, color_head = TRUE)
}
\arguments{

  \item{cr}{correlated regions returned from \code{\link{cr_enrichedheatmap}}}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{cytoband_list}{a list of cytoband returned by \code{\link{gtrellis_cr_genes}}}
  \item{color_head}{internal use}

}
\details{
For each cytoband, there are several tracks:

\itemize{
  \item points showing correlations
  \item mean correlation in 50kb window
  \item genes
}
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
