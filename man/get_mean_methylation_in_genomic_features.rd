\name{get_mean_methylation_in_genomic_features}
\alias{get_mean_methylation_in_genomic_features}
\title{
Calculate mean methylation in a list of genomic features
}
\description{
Calculate mean methylation in a list of genomic features
}
\usage{
get_mean_methylation_in_genomic_features(sample_id, genomic_features, chromosome = paste0("chr", 1:22))
}
\arguments{

  \item{sample_id}{a vector of sample IDs}
  \item{genomic_features}{a list or a single \code{\link[GenomicRanges]{GRanges}} objects}
  \item{chromosome}{a vector of chromosome names}

}
\value{
A list of or a single \code{\link[GenomicRanges]{GRanges}} objects (according to \code{genomic_features} you specified) in which mean methylation matrix and number of CpG in each region
are attached. The variable can be sent to \code{\link{heatmap_diff_methylation_in_genomic_features}} to visualize.

Note it should be kept in mind that it doesn't make any sense to calculate mean methylation in long regions where
there are hetergenuous methylation patterns.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
