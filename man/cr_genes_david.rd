\name{cr_genes_david}
\alias{cr_genes_david}
\title{
DAVID analysis for CR genes
}
\description{
DAVID analysis for CR genes
}
\usage{
cr_genes_david(cr, david_user)
}
\arguments{

  \item{cr}{correlated regions returned by \code{\link{cr_enrichedheatmap}}}
  \item{david_user}{username for DAVID API (\url{https://david.ncifcrf.gov/content.jsp?file=WS.html)}}

}
\details{
Genes in k-means group 1 and 4 are sent to DAVID web server. The significant functions are visualized
as a heatmap.
}
\value{
a list of function enrichments
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
