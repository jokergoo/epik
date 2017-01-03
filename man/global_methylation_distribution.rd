\name{global_methylation_distribution}
\alias{global_methylation_distribution}
\title{
Global methylation distribution
}
\description{
Global methylation distribution
}
\usage{
global_methylation_distribution(sample_id, annotation,
    annotation_color = structure(seq_along(unique(annotation)), names = unique(annotation)),
    reorder_column = TRUE, ha = NULL, chromosome = paste0("chr", 1:22), by_chr = FALSE, max_cov = 100,
    background = NULL, p = 0.001)
}
\arguments{

  \item{sample_id}{a vector of sample ids}
  \item{annotation}{classification information}
  \item{annotation_color}{color for classifications}
  \item{reorder_column}{whether reorder the samples}
  \item{ha}{additional annotation can be specified as a \code{\link[ComplexHeatmap]{HeatmapAnnotation}} object}
  \item{chromosome}{chromosomes}
  \item{by_chr}{whether make the plot by chromosome}
  \item{max_cov}{maximum coverage (used to get rid of extremely high coverage which affects visualization of CpG coverage distribution)}
  \item{background}{background to look into. The value can be a single \code{\link[GenomicRanges]{GRanges}} object or a list of \code{\link[GenomicRanges]{GRanges}} objects.}
  \item{p}{probability to randomly sample CpG sites}

}
\details{
It visualize distribution of methylation values and CpG coverages through heatmaps.
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
