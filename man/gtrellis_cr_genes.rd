\name{gtrellis_cr_genes}
\alias{gtrellis_cr_genes}
\title{
Visualize CR genes in gtrellis layout
}
\description{
Visualize CR genes in gtrellis layout
}
\usage{
gtrellis_cr_genes(cr, txdb, expr, species = "hg19")
}
\arguments{

  \item{cr}{correlated regions, should be returned by \code{\link{cr_enrichedheatmap}}}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{expr}{expression matrix which was used in \code{\link{correlated_regions}}}
  \item{species}{species}

}
\details{
CR genes in k-means group 1 and 4 are visualized in gtrellis layout. Cytobands
which are significantly overlapped by CR genes are highlighted.
}
\value{
A list of two elements which shows how each cytoband are overlapped by CR genes
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
