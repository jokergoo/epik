\name{sig_cytoband_gtrellis}
\alias{sig_cytoband_gtrellis}
\title{
Visualize correlations in cytoband
}
\description{
Visualize correlations in cytoband
}
\usage{
sig_cytoband_gtrellis(cr, txdb, cytoband_list, color_head = TRUE)
}
\arguments{

  \item{cr}{correlated regions returned from \code{\link{cr_enriched_heatmap}}}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{cytoband_list}{a list of cytoband returned by \code{\link{gtrellis_cr_genes}}}
  \item{color_head}{internal use}

}
\details{
This function visualizes significant cytobands which have been found in \code{\link{gtrellis_cr_genes}}.

For each cytoband, there are several tracks:

\itemize{
  \item cytoband name. Green corresponds to significant cytobands in group 1 and red corresponds to group 4
  \item points showing correlations
  \item a one row heamtap showing mean correlation in 50kb window
  \item genes, green lines are genes in cluster 1, red lines are genes in cluster 4.
}
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
