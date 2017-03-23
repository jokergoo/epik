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
    cor_method = "spearman", subgroup = NULL, window_size = 5, window_step = window_size,
    max_width = 10000, raw_meth = FALSE, cov_cutoff = 3, min_dp = 4, col = NULL, species = "hg19")
}
\arguments{

  \item{sample_id}{a vector of sample IDs}
  \item{expr}{expression matrix}
  \item{txdb}{a \code{\link[GenomicFeatures]{TxDb-class}} object.}
  \item{chr}{a single chromosome name}
  \item{extend}{extension of gene model, both upstream and downstream}
  \item{cov_filter}{if \code{coverage} hook is set in \code{\link{methylation_hooks}}, this option can be set to filter out CpG sites with low coverage across samples. the value for this option is a function for which the argument is a vector of coverage values for current CpG in all samples. The default setting means the CpG should have coverage in more than half of samples.}
  \item{cor_method}{method for calcualting correlations, pass to \code{\link[stats]{cor}}.}
  \item{subgroup}{subgroup information. If provided, ANOVA test and group mean are applied on each correlated region.}
  \item{window_size}{how many CpG sites in a window}
  \item{window_step}{step of the sliding window, measured in number of CpG sites}
  \item{max_width}{maximum width of a window}
  \item{raw_meth}{whether use raw methylation value (values from \code{raw} hook set in \code{\link{methylation_hooks}})}
  \item{cov_cutoff}{cutoff for CpG coverage when using raw methylation rate, used for raw methylation. Note when the CpG coverage is too low, the raw methylation rate is not reliable. Raw methylation rate for those CpGs with coverage less this this cutoff is set to \code{NA} will be further filtered by \code{min_gp}.}
  \item{min_dp}{minimal number of non-NA values for calculating correlations. When \code{meth} is the raw methylation, values for which CpG coverage is too low will be replaced with \code{NA}, We only use non-NA values to calculate correlations. If the number of data points for calculating correlation is less than \code{min_dp}, the CpG window is just excluded.}
  \item{col}{color for subgroups. This setting will be saved in the returned object and will be used in downstream analysis. If not set, random colors are assigned.}
  \item{species}{species. This setting will be saved and used in downstream analysis}

}
\details{
A correlated region is defined as a region where methylation is correlated with the expression of associated gene.
The detection for correlated regions is gene-centric. For every gene, the processes are as follows:

\itemize{
  \item extend to both upstream and downstream by \code{extend};
  \item filter CpG sites by CpG coverage (by \code{cov_filter});
  \item from the most upstream, use a sliding window which contains \code{windows_size} CpG sites, moving step of \code{window_step} CpG sites;
  \item calculate correlation between methylation and gene expression for this window;
  \item calculate other statistics.
}

Following meta columns are attached to the \code{\link[GenomicRanges]{GRanges}} objects:

\describe{
  \item{ncpg}{number of CpG sites}
  \item{mean_meth_*}{mean methylation in each window in every sample.}
  \item{corr}{correlation between methylation and expression}
  \item{corr_p}{p-value for the correlation test}
  \item{meth_IQR}{IQR of mean methylation if \code{subgroup} is not set}
  \item{meth_anova}{p-value from oneway ANOVA test if \code{subgroup} is set}
  \item{meth_diameter}{range between maximum mean and minimal mean in all subgroups if \code{subgroup} is set}
  \item{meth_diff}{when there are two subgroups, the mean methylation in subgroup 1 substracting mean methylation in subgroup 2.}
  \item{gene_id}{gene id}
  \item{gene_tss_dist}{distance to tss of genes}
  \item{tx_tss_dist}{if genes have multiple transcripts, this is the distance to the nearest transcript}
  \item{nearest_txx_tss}{transcript id of the nearest transcript}
}

This function keeps all the information for all CpG windows. Uses can use \code{\link{cr_add_fdr_column}} to add fdr columns to the object,
filter significant correlated regions by p-value, fdr and meth_diff columns, or use \code{\link{cr_reduce}} to reduce the significant regions.
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object which contains correlations and associated statistics for every CpG windows.

The settings for finding correlated regions are stored as the meta data of the \code{\link[GenomicRanges]{GRanges}} object.
}
\seealso{
Internally, the calculation is done by \code{\link{correlated_regions_by_window}}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
