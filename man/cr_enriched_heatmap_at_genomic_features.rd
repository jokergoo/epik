\name{cr_enriched_heatmap_at_genomic_features}
\alias{cr_enriched_heatmap_at_genomic_features}
\title{
Visualizing enrichment for epigenomic signals at TSS-CGIs
}
\description{
Visualizing enrichment for epigenomic signals at TSS-CGIs
}
\usage{
cr_enriched_heatmap_at_genomic_features(cr, txdb, expr, gf,
    fdr_cutoff = 0.05, meth_diff_cutoff = 0.1, marks = NULL, type = "neg", extend = 5000,
    min_reduce = 1, min_width = 1000, nearest_by = "tss", expr_ha)
}
\arguments{

  \item{cr}{correalted regions}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{expr}{expression matrix which was used in \code{\link{correlated_regions}}}
  \item{gf}{genomic features, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{fdr_cutoff}{cutoff for fdr}
  \item{meth_diff_cutoff}{cutoff for methylation difference. If there are no subgroup information or only one subgroup, \code{meth_IQR} column is used for filtering. If there are more than one subgroups, \code{meth_diameter} column is used for filtering.}
  \item{marks}{names of histone marks, should be supported in \code{\link{chipseq_hooks}}}
  \item{type}{visualize negative correlated regions or positive correlated regions}
  \item{extend}{base pairs extended to upstream and downstream}
  \item{min_reduce}{base pairs for merging neighbouring regions}
  \item{min_width}{minimal width of regions}
  \item{nearest_by}{"tss" or "gene", how to connect genomic features to genes}
  \item{expr_ha}{a \code{\link[ComplexHeatmap]{HeatmapAnnotation}} class object.It is used for the expression heatmap}

}
\details{
There are several heatmaps visualize various signals enriched at genomic features. After annotate to genes,
in the extended regions, each region can only have one gene.

There are following heatmaps:

\itemize{
  \item heatmap for gene expression
  \item If \code{cr} is returned form \code{\link{cr_enriched_heatmap}}, there is a one column heatmap which shows the k-means cluters genes belong to
  \item heatmap for significant correlated regions
  \item a point plot showing region length
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
gf regions.
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
