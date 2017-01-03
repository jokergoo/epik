\name{filter_correlated_regions}
\alias{filter_correlated_regions}
\title{
Get correlated regions with significant correlations
}
\description{
Get correlated regions with significant correlations
}
\usage{
filter_correlated_regions(chromosome = paste0("chr", 1:22), template,
    cutoff = 0.05, adj_method = "BH", meth_diameter_cutoff = 0.25, meth_IQR_cutoff = 0.25,
    anova_cutoff = 0.05)
}
\arguments{

  \item{chromosome}{a vector of chromosome names}
  \item{template}{template path to find cr files}
  \item{cutoff}{cutoff of adjusted correlation p-values}
  \item{adj_method}{method for calculating adjusted p-values}
  \item{meth_diameter_cutoff}{cutoff for methylation diameters}
  \item{meth_IQR_cutoff}{cutoff for IQR, if there is no subtype information, IQR is used to remove less variable methylation}
  \item{anova_cutoff}{cutoff for adjust ANOVA p-values}

}
\details{
As explained in \code{\link{correlated_regions}}, original cr is huge and is always saved as a separated file for single chromosome.
Here \code{template} defined how to get the cr files. E.g. if \code{template}is defined as "path_to/@{chr}_cr.rds", the funciton
will replace "@{chr}" to every chromosome and read the data.

Two additional columns are attached:

\describe{
  \item{corr_fdr}{FDR for correlation p-values}
  \item{meth_anova_fdr}{FDR for anova test, will be added only if \code{factor} is set in \code{\link{correlated_regions}}.}
}
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
