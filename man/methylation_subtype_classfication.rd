\name{methylation_subtype_classfication}
\alias{methylation_subtype_classfication}
\title{
Classify subtypes by methylation data
}
\description{
Classify subtypes by methylation data
}
\usage{
methylation_subtype_classfication(gr, n_class, pct_cutoff = 1, corr_cutoff = 0.5,
    k, ha = NULL, cgi = NULL, shore = NULL)
}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} which contains mean methylation, should be generated by \code{\link{get_mean_methylation_in_genomic_features}}}
  \item{n_class}{number of classes expected}
  \item{pct_cutoff}{percent of most variable rows}
  \item{corr_cutoff}{cutoff for absolute correlation}
  \item{k}{number of correlated windows}
  \item{ha}{additional annotation}
  \item{cgi}{a \code{\link[GenomicRanges]{GRanges}} object which contains CpG islands}
  \item{shore}{a \code{\link[GenomicRanges]{GRanges}} object which contains CpG shores}

}
\details{
For the subtype classification which is based on clustering, if there are clear subtypes, 
it is expected that there must be a group of rows that show high correlation to each other. 
Based on this correlation feature, we select rows that under cutoff of \code{corr_cutoff}, 
each row should correlate to at least other \code{k} rows. On the second hand, since difference between 
subtypes are not in an identical position, we first separate samples into two groups based on consensus clustering, 
then, for the subgroup which contains more samples, we separate them again into two subgroups. 
We apply it repeatedly until there are \code{n_class} subtypes. On every step of clustering, 
we select rows based on the correlation criterion and the final rows are union of rows in all iterations.

CpG islands and shores will be added as row annotations to the heatmap.
}
\value{
A vector with predicted classification of samples
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}