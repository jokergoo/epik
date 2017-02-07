\name{hilbert_curve_chipseq_difference}
\alias{hilbert_curve_chipseq_difference}
\title{
Visualize ChIP-Seq signal difference by Hilbert curve
}
\description{
Visualize ChIP-Seq signal difference by Hilbert curve
}
\usage{
hilbert_curve_chipseq_difference(mark, subgroup, comparison, chromosome = paste0("chr", 1:22),
    species = "hg19", type = c("global_mean", "subgroup_mean", "abs_difference", "rel_difference"),
    density_column = "density")
}
\arguments{

  \item{mark}{name of the histone mark, should also be supported in \code{\link{chipseq_hooks}}}
  \item{subgroup}{subgroup information which corresponds to sample IDs stored in \code{\link{methylation_hooks}}. The value should be a named vector that names are sample IDs.}
  \item{comparison}{if there are more than one subgroups, the comparison of two subgroups which shows the methylation difference The value is a vector of length of two and the difference is calculated as \code{subgroup[1] - subgroup[2]}}
  \item{chromosome}{a vector fo chromosome names}
  \item{species}{species}
  \item{type}{four types of visualization supported, see "details" section}
  \item{density_column}{the column name of the signal defined in \code{\link{chipseq_hooks}}$peak}

}
\details{
Genome is segmented by 1kb window and mean signal for each 1kb window is calculated, later visualized
by Hilbert curve.

There are four types of visualization methods:

\describe{
  \item{global_mean}{the mean signal averaged from all samples}
  \item{subgroup_mean}{the mean signal averaged in every subgroup}
  \item{abs_difference}{the absolute difference of signal in two subgroups}
  \item{rel_difference}{the relative difference in two subgroups. The value is calculated as absolute difference divided by mean singal.}
}
}
\value{
a \code{\link[GenomicRanges]{GRanges}} object which contains mean signal for the 1kb segments and other statistics.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
