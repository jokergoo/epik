\name{cr_genes_function_enrichment}
\alias{cr_genes_function_enrichment}
\title{
Functional enrichment for CR genes by DAVID
}
\description{
Functional enrichment for CR genes by DAVID
}
\usage{
cr_genes_function_enrichment(cr, david_user, count_cutoff = 50, fdr_cutoff = 0.01,
    pop_count_cutoff = 5000)
}
\arguments{

  \item{cr}{correlated regions returned by \code{\link{cr_enriched_heatmap}}}
  \item{david_user}{username for DAVID API (\url{https://david.ncifcrf.gov/content.jsp?file=WS.html} )}
  \item{count_cutoff}{minimal number of CR genes in a function term}
  \item{fdr_cutoff}{cutoff of fdr of the enrichment test}
  \item{pop_count_cutoff}{maximum number of population genes in a function term.}

}
\details{
Genes in k-means group 1 and 4 are sent to DAVID web server to do functional enrichment. 
The significant functions are visualized as a heatmap.

Only three Gene Ontology (biological process, molecular function and cellular component) categories are used.

There is also a heatmap which shows the significant enrichment.
}
\value{
A list of function enrichments.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
