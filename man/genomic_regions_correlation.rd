\name{genomic_regions_correlation}
\alias{genomic_regions_correlation}
\title{
Correlation between two sets of genomic regions
}
\description{
Correlation between two sets of genomic regions
}
\usage{
genomic_regions_correlation(gr_list_1, gr_list_2, background = NULL,
    chromosome = paste0("chr", 1:22), species = "hg19",
    nperm = 0, mc.cores = 1, stat_fun = genomic_corr_jaccard, ...,
    bedtools_binary = Sys.which("bedtools"), tmpdir = tempdir())
}
\arguments{

  \item{gr_list_1}{a list of \code{\link[GenomicRanges]{GRanges}} objects, should be a named list, e.g. low methylated regions in different samples.}
  \item{gr_list_2}{a list of \code{\link[GenomicRanges]{GRanges}} objects, should be a named list, e.g. a list of different genomic features.}
  \item{background}{a \code{\link[GenomicRanges]{GRanges}} object. The correlation is only looked in background regions.}
  \item{chromosome}{a vector of chromosome names}
  \item{species}{species, used for random shuffling genomic regions}
  \item{nperm}{number of random shufflings. If it is set to 0, no random shuffling will be performed.}
  \item{mc.cores}{number of cores for parallel calculation}
  \item{stat_fun}{method to calculate correlations. There are some pre-defined functions: \code{\link{genomic_corr_reldist}}, \code{\link{genomic_corr_absdist}} measure how two sets of genomic regions are close; \code{\link{genomic_corr_jaccard}}, \code{\link{genomic_corr_nintersect}}, \code{\link{genomic_corr_pintersect}}, \code{\link{genomic_corr_sintersect}} measures how two sets of genomic regions are overlapped. The self-defined function should accept at least two arguments which are two GRanges object. The third argument is \code{...} which is passed from the main function. The function should only return a numeric value.}
  \item{...}{pass to \code{stat_fun}}
  \item{bedtools_binary}{random shuffling is perfomed by \code{bedtools}. If \code{bedtools} is not in \code{PATH}, the path of \code{bedtools} can be set here.}
  \item{tmpdir}{dir for temporary files}

}
\details{
The correlation between two sets of genomic regions basically means how much the first type of genomic regions
are overlapped or close to the second type of genomic regions.

The significance of the correlation is calculated by random shuffling the regions. 
In random shuffling, regions in \code{gr_list_1} will be shuffled. So if you want to shuffle \code{gr_list_2},
just switch the first two arguments.

Pleast note random shuffling is done by bedtools, so bedtools should be installed and exists in \code{PATH}
and should support \code{-i -g -incl} options.
}
\value{
A list containing following elements:

\describe{
  \item{stat}{statistic value}
  \item{fold_change}{stat/E(stat), stat divided by expected value which is generated from random shuffling}
  \item{p.value}{p-value for over correlated. So, 1 - p.value is the significance for being less correlated}
  \item{stat_random_mean}{mean value of stat in random shuffling}
  \item{stat_random_sd}{standard deviation in random shuffling}
}

If \code{perm} is set to 0 or 1, \code{fold_change}, \code{p.value}, \code{stat_random_mean} and \code{stat_random_sd} are all \code{NULL}.
}
\seealso{
\code{\link{genomic_corr_reldist}}, \code{\link{genomic_corr_jaccard}}, \code{\link{genomic_corr_absdist}}, \code{\link{genomic_corr_nintersect}}, 
\code{\link{genomic_corr_pintersect}}, \code{\link{genomic_corr_sintersect}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
genomic_regions_correlation(gr1, gr2, nperm = 0)
genomic_regions_correlation(list(gr1 = gr1), list(gr2 = gr2), nperm = 0)
}
