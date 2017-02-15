\name{cr_enriched_heatmap_at_tss}
\alias{cr_enriched_heatmap_at_tss}
\title{
Visualizing enrichment for epigenomic signals at TSS
}
\description{
Visualizing enrichment for epigenomic signals at TSS
}
\usage{
cr_enriched_heatmap_at_tss(cr, txdb, expr, cgi, fdr_cutoff = 0.05,
    meth_diff_cutoff = 0.1, marks = NULL, type = "neg", extend = c(5000, 10000),
    expr_annotation)
}
\arguments{

  \item{cr}{correalted regions}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{expr}{expression matrix which was used in \code{\link{correlated_regions}}}
  \item{cgi}{CpG island, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{fdr_cutoff}{cutoff for fdr}
  \item{meth_diff_cutoff}{cutoff for methylation difference. If there are no subgroup information or only one subgroup, \code{meth_IQR} column is used for filtering. If there are more than one subgroups, \code{meth_diameter} column is used for filtering.}
  \item{marks}{names of histone marks}
  \item{type}{visualize negative correlated regions or positive correlated regions}
  \item{extend}{base pairs extended to upstream and downstream}
  \item{expr_annotation}{a \code{\link[ComplexHeatmap]{HeatmapAnnotation}} class object.It is used for the expression heatmap}

}
\details{
There are several heatmaps visualize various signals enriched at gene TSS.

There are following heatmaps:

\itemize{
  \item heatmap for gene expression
  \item a point plot showing gene length
  \item If \code{cr} is returned form \code{\link{cr_enriched_heatmap}}, there is a one column heatmap which shows the k-means cluters genes belong to
  \item heatmap for CGI enrichment at TSS
  \item heatmap for significant correlated regions
  \item heatmap for correlation between methylation and gene expression
  \item heatmap for mean methylation
  \item heatmap for metnylation difference
  \item heatmap for correlation, mean signal and signal difference for histone marks
}

If there are more than 12 heatmaps, they will be put into two pages.

Heatmaps are split into two sub-clusters by k-means clustering on the mean methylation matrix.
If there are two subgroups in all samples, each subcluster are split by high expression/low expression
in subgroup 1. In each high expression/low expression, rows are split by the k-means clusters calculated
in \code{\link{cr_enriched_heatmap}}. Finally, rows are clustered by considering closeness of signals in the extended
TSS regions.
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
