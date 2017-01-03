\name{enriched_heatmap_list_on_tss_cgi}
\alias{enriched_heatmap_list_on_tss_cgi}
\title{
Enriched heatmaps to visualize how signals are at enriched CpG islands
}
\description{
Enriched heatmaps to visualize how signals are at enriched CpG islands
}
\usage{
enriched_heatmap_list_on_tss_cgi(cr, cgi, txdb, expr, hm_list = NULL, hm_name = NULL, by = "gene", ...)
}
\arguments{

  \item{cr}{filtered correlated regions from \code{\link{filter_correlated_regions}}}
  \item{cgi}{a \code{\link[GenomicRanges]{GRanges}} object which contains CpG islands}
  \item{txdb}{a \code{GenomicFeatures::GRanges} object.}
  \item{expr}{the expression matrix which is same as in \code{\link{correlated_regions}}}
  \item{hm_list}{a list of ChIP-Seq peaks}
  \item{hm_name}{names for the peaks}
  \item{by}{by \code{gene} or \code{tx}}
  \item{...}{pass to \code{\link[EnrichedHeatmap]{draw,EnrichedHeatmapList-method}}}

}
\details{
Following signals are visualized around CpG islands which are close to gene/ts tss:

\itemize{
  \item correlated regions (positive correlated regions and negative correlated regions can be subsetted by \code{corr} column in the object)
  \item width of CpG islands
  \item strand of associated tss
  \item number of tss in the extended area
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
