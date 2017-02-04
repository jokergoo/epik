\name{global_methylation_distribution}
\alias{global_methylation_distribution}
\title{
Global methylation distribution
}
\description{
Global methylation distribution
}
\usage{
global_methylation_distribution(sample_id, subgroup,
    reorder_column = TRUE,
    ha = HeatmapAnnotation(subgroup = subgroup, show_annotation_name = TRUE),
    chromosome = paste0("chr", 1:22), by_chr = FALSE,
    # max_cov = 100,
    background = NULL, p = NULL, meth_range = c(0, 1))
}
\arguments{

  \item{sample_id}{a vector of sample ids}
  \item{subgroup}{subgroup information}
  \item{reorder_column}{if it is true, samples are first ordered by subgroups and in each subgroup, samples are ordered by median values}
  \item{ha}{additional annotation can be specified as a \code{\link[ComplexHeatmap]{HeatmapAnnotation}} object}
  \item{chromosome}{chromosomes}
  \item{by_chr}{whether make the plot by chromosome}
  \item{background}{background to look into. The value can be a single \code{\link[GenomicRanges]{GRanges}} object or a list of \code{\link[GenomicRanges]{GRanges}} objects.}
  \item{p}{probability to randomly sample CpG sites}
  \item{meth_range}{the range of methylation on the plot}

}
\details{
The distribution density is visualized as heatmaps.
}
\value{
If \code{by_chr} is set to \code{FALSE}, it returns a vector of column order.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
