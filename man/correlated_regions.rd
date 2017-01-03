\name{correlated_regions}
\alias{correlated_regions}
\title{
Correlation between methylation and expression
}
\description{
Correlation between methylation and expression
}
\usage{
correlated_regions(sample_id, expr, txdb, chr, extend = 50000,
    cov_filter = function(x) sum(x > 0, na.rm = TRUE) > length(x)/2,
    cor_method = "spearman", factor = NULL, window_size = 5, max_width = 10000,
    raw_meth = FALSE, cov_cutoff = 3, min_dp = 4, col = NULL)
}
\arguments{

  \item{sample_id}{a vector of sample id}
  \item{expr}{expression matrix in which columns correspond to sample ids}
  \item{txdb}{a \code{GenomicFeatures::GRanges} object. Gene names should be same type as row names in \code{expr}}
  \item{chr}{a single chromosome}
  \item{extend}{extension of gene model, both upstream and downstream}
  \item{cov_filter}{if \code{coverage} hook is set in \code{\link{methylation_hooks}}, this option can be set to filter out CpG sites with low coverage across samples. the value for this option is a function for which the argument is a vector of coverage values for current CpG in all samples.}
  \item{cor_method}{method to calculate correlation}
  \item{factor}{classes of samples}
  \item{window_size}{number of CpGs in a window}
  \item{max_width}{maximum width of a window}
  \item{raw_meth}{whether use raw methylation value (values from \code{raw} hook set in \code{\link{methylation_hooks}})}
  \item{cov_cutoff}{cutoff for coverage}
  \item{min_dp}{minimal non-NA values for calculating correlations}
  \item{col}{color for classes}

}
\details{
The detection for correlated regions is gene-centric. For every gene, the process are as follows:

\itemize{
  \item extend to both upstream and downstream
  \item from the most upstream, use a sliding window which contains \code{windows_size} CpG sites
  \item filter each window by CpG coverage (by \code{cov_filter} and \code{cov_cutoff})
  \item calculate correlation between methylation and gene expression for this window
}

Following meth columns are attached to the \code{\link[GenomicRanges]{GRanges}} objects:

\describe{
  \item{n}{number of CpG sites}
  \item{mean_meth_*}{mean methylation in each window in every sample.}
  \item{corr}{correlation}
  \item{corr_p}{p-value for the correlation test}
  \item{meth_IQR}{IQR of mean methylation if \code{factor} is not set}
  \item{meth_anova}{p-value from oneway ANOVA test if \code{factor} is set}
  \item{meth_diameter}{range between maximum mean and minimal mean in all subgroups if \code{factor} is set}
  \item{gene_id}{gene id}
  \item{gene_tss_dist}{distance to tss of genes}
  \item{tx_tss_dist}{if genes have multiple transcripts, this is the distance to the nearest transcript}
  \item{nearest_txx_tss}{transcript id of the nearest transcript}
}

This function keeps all the information for all CpG windows. Users can get \code{\link{filter_correlated_regions}} to get correlated regions
with significant correlations and use \code{\link{reduce_cr}} to merge neighbouring windows.

Since information for all CpG windows are kept, the size of the object is always very huge, thus, it is reasonable
to analyze each chromosome separately and save each object as a single file. Some downstream functions expect a formatted
path of the cr file.
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object which contains associated statistics for every CpG windows.
}
\seealso{
\code{\link{filter_correlated_regions}}, \code{\link{reduce_cr}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
