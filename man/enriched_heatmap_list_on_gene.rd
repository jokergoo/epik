\name{enriched_heatmap_list_on_gene}
\alias{enriched_heatmap_list_on_gene}
\title{
Enriched heatmaps to visualize how signals are enriched at genes
}
\description{
Enriched heatmaps to visualize how signals are enriched at genes
}
\usage{
enriched_heatmap_list_on_gene(cr, cgi, txdb, expr, hm_list = NULL, hm_name = NULL, on = "tss", by = "gene",
    hm_cor_p_cutoff = 0.05, show_expr = TRUE, ...)
}
\arguments{

  \item{cr}{filtered correlated regions from \code{\link{filter_correlated_regions}}}
  \item{cgi}{a \code{\link[GenomicRanges]{GRanges}} object which contains CpG islands}
  \item{txdb}{a \code{GenomicFeatures::GRanges} object.}
  \item{expr}{the expression matrix which is same as in \code{\link{correlated_regions}}}
  \item{hm_list}{a list of ChIP-Seq peaks}
  \item{hm_name}{names for the peaks}
  \item{on}{where the signals are enriched in, possible values are \code{tss} and \code{body}}
  \item{by}{if \code{on} is set to \code{tss}, whether it is tss of \code{gene} or \code{tx}}
  \item{hm_cor_p_cutoff}{cutoff for the correlation between peak intensity and gene expression}
  \item{show_expr}{whether show heatmap of gene expression}
  \item{...}{pass to \code{\link[EnrichedHeatmap]{draw,EnrichedHeatmapList-method}}}

}
\details{
Following signals are visualized around gene/tx tss or gene bodies:

\itemize{
  \item correlated regions (positive correlated regions and negative correlated regions can be subsetted by \code{corr} column in the object)
  \item correlation between peak intensity and gene expression
  \item overlapping with CpG islands
  \item length of genes
  \item relative gene expression
  \item expression level
  \item mean intensity of peaks in each subgroup
  \item mean methylation in each subgroup
}
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dfkz.de>
}
\examples{
# There is no example
NULL

}
