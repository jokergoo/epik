\name{heatmap_diff_methylation_in_genomic_features}
\alias{heatmap_diff_methylation_in_genomic_features}
\title{
Heatmap for differential methylation in genomic features
}
\description{
Heatmap for differential methylation in genomic features
}
\usage{
heatmap_diff_methylation_in_genomic_features(gr, subgroup,
    ha = HeatmapAnnotation(subgroup = subgroup, show_annotation_name = TRUE),
    genomic_features = NULL,
    meth_diff = 0, cutoff = 0.05, adj_method = "BH",
    cluster_columns = c("subgroup", "all", "none"), ...)
}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object returned from \code{\link{get_mean_methylation_in_genomic_features}}}
  \item{subgroup}{subgroup information}
  \item{ha}{column annotations, a \code{\link[ComplexHeatmap]{HeatmapAnnotation}} object}
  \item{genomic_features}{a single or a list of \code{\link[GenomicRanges]{GRanges}} ojects}
  \item{meth_diff}{minimal range between mean value in subgroups}
  \item{cutoff}{if classification information is provided, p-value for the oneway ANOVA test}
  \item{adj_method}{how to calculate adjusted p-values}
  \item{cluster_cols}{how to cluster columns}
  \item{...}{pass to \code{\link[ComplexHeatmap]{Heatmap}}}

}
\details{
Regions have differential methylation are only visualized.
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object which only contains regions with significant differential methylation.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
