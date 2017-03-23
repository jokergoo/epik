\name{cr_genes_gtrellis}
\alias{cr_genes_gtrellis}
\title{
Visualize CR genes in gtrellis layout
}
\description{
Visualize CR genes in gtrellis layout
}
\usage{
cr_genes_gtrellis(cr, txdb, expr)
}
\arguments{

  \item{cr}{correlated regions, should be returned by \code{\link{cr_enriched_heatmap}}}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{expr}{expression matrix which was used in \code{\link{correlated_regions}}}

}
\details{
CR genes in k-means group 1 and 4 (which correspond to negative correlated gene body
and positive correlated gene body) are visualized in gtrellis layout. Cytobands
which are significantly overlapped by CR genes are highlighted. In gtrellis layout, there
are following tracks:

\itemize{
  \item rainfall plot for genes in k-means cluster 1. The y-axis corresponds to the minimal distance to neighbouring genes. The color of points corresponds to the epxression value (blue is low expression and red is high expression) and size of points corresponds to gene length.
  \item a one row heatmap showing how much each cytoband is covered by CR genes in cluster 1.
  \item rainfall plot for genes in k-means cluster 4
  \item a one row heatmap showing how much each cytoband is covered by CR genes in cluster 4.
  \item cytoband
}
}
\value{
A list of two elements which shows how each cytoband is overlapped by CR genes
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
