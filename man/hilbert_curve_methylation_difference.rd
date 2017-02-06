\name{hilbert_curve_methylation_difference}
\alias{hilbert_curve_methylation_difference}
\title{
Visualize methylation by Hilbert curve
}
\description{
Visualize methylation by Hilbert curve
}
\usage{
hilbert_curve_methylation_difference(subgroup, comparison, chromosome = paste0("chr", 1:22),
    species = "hg19", type = c("global_mean", "subgroup_mean", "difference"))
}
\arguments{

  \item{subgroup}{subgroup information which corresponds to sample IDs stored in \code{\link{methylation_hooks}}. The value can be a vector with same length as sample IDs or a named vector that names are sample IDs.}
  \item{comparison}{if there are more than one subgroups, the comparison of two subgroups which shows the methylation difference The value is a vector of length of two and the difference is calculated as subgroup[1] - subgroup[2]}
  \item{chromosome}{a vector of chromosome}
  \item{species}{species}
  \item{type}{Three types of visualization supported, see "details" section}

}
\details{
Genome is segmented by 1kb window and mean methylation for each 1kb window is calculated, later visualized
by Hilbert curve.

There are three types of visualization methods:

\describe{
  \item{global_mean}{the mean methylation averaged from all samples}
  \item{subgroup_mean}{the mean methylation averaged in every subgroup}
  \item{difference}{the difference of methylation in two subgroups}
}
}
\value{
a \code{\link[GenomicRanges]{GRanges}} object which contains mean methylation for the 1kb segmentation and other statistics.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
