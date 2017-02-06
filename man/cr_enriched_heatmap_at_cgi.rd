\name{cr_enriched_heatmap_at_cgi}
\alias{cr_enriched_heatmap_at_cgi}
\title{
Visualizing enrichment for epigenomic signals at TSS-CGIs
}
\description{
Visualizing enrichment for epigenomic signals at TSS-CGIs
}
\usage{
cr_enriched_heatmap_at_cgi(cr, txdb, expr, cgi,
    fdr_cutoff = 0.05, meth_diff_cutoff = 0.1, marks = NULL, type = "neg", extend = 5000,
    expr_annotation)
}
\arguments{

  \item{cr}{correalted regions}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{expr}{expression matrix which was used in \code{\link{correlated_regions}}}
  \item{cgi}{CpG island}
  \item{fdr_cutoff}{cutoff for fdr of correlation p-value and anove p-values}
  \item{meth_diff_cutoff}{cutoff for methylation difference}
  \item{marks}{names of histone marks}
  \item{type}{use negative correlated regions or positive correlated regions}
  \item{extend}{base pairs extended to upstream and downstream}
  \item{expr_annotation}{a \code{\link[ComplexHeatmap]{HeatmapAnnotation}} class object}

}
\details{
There are several heatmaps visualize various signals enriched at TSS-CGIs.

\itemize{
  \item heatmap for gene expression
  \item If \code{cr} is returned form \code{\link{cr_enrichedheatmap}}, there is a one column heatmap which shows the k-means cluters genes belong to
  \item heatmap for correlated regions
  \item a point plot 
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
