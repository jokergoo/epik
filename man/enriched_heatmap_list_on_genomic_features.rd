\name{enriched_heatmap_list_on_genomic_features}
\alias{enriched_heatmap_list_on_genomic_features}
\title{
Enriched heatmaps to visualize how signals are enriched at a certain genomic feature
}
\description{
Enriched heatmaps to visualize how signals are enriched at a certain genomic feature
}
\usage{
enriched_heatmap_list_on_genomic_features(cr, gf, hm_list = NULL, hm_name = NULL, ...)
}
\arguments{

  \item{cr}{filtered correlated regions from \code{\link{filter_correlated_regions}}}
  \item{gf}{genomic features, e.g. TFBS or enhancers, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{hm_list}{a list of ChIP-Seq peaks}
  \item{hm_name}{names for the peaks}
  \item{...}{pass to \code{\link[EnrichedHeatmap]{draw,EnrichedHeatmapList-method}}}

}
\details{
Following signals are visualized around specified genomic features:

\itemize{
  \item correlated regions (positive correlated regions and negative correlated regions can be subsetted by \code{corr} column in the object)
  \item width of genomic features
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
