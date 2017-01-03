\name{get_mean_methylation_in_genomic_features}
\alias{get_mean_methylation_in_genomic_features}
\title{
Calculate mean methylation value in a list of genomic features
}
\description{
Calculate mean methylation value in a list of genomic features
}
\usage{
get_mean_methylation_in_genomic_features(sample_id, gf_list, average = TRUE, p = 0.001,
    chromosome = paste0("chr", 1:22),
    filter_fun = function(s) length(s) >= 10 && any(diff(s) < 50))
}
\arguments{

  \item{sample_id}{a vector of sample IDs}
  \item{gf_list}{a list of genomic features in \code{\link[GenomicRanges]{GRanges}} class}
  \item{average}{whether to calcualte average methylation in a interval? if not, the function will randomly sample CpG sites from the input regions.}
  \item{p}{if average is \code{FALSE}, the probability to randomly pick CpG sites}
  \item{chromosome}{a vector of chromosome names}
  \item{filter_fun}{filtering function on CpG sites in each intersection (e.g. exclude regions which contains too few CpGs). The object sent to this function is a vector of positions for CpGs that locate in the current region.}

}
\value{
A list of \code{\link[GenomicRanges]{GRanges}} objects in which mean methylation matrix are attached.
The variable can be sent to \code{\link{heatmap_diff_methylation_in_genomic_features}} to visualize.

= author
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
