\name{cr_overlap_to_genomic_features}
\alias{cr_overlap_to_genomic_features}
\title{
Enrichment of cr to other genomic features
}
\description{
Enrichment of cr to other genomic features
}
\usage{
cr_overlap_to_genomic_features(cr, gf_list, species = NULL, chromosome = paste0("chr", 1:22))
}
\arguments{

  \item{cr}{filtered correlated regions from \code{\link{filter_correlated_regions}}}
  \item{gf_list}{a list of \code{\link[GenomicRanges]{GRanges}} objects}
  \item{species}{species}
  \item{chromosome}{a vector of chromosomes}

}
\details{
There will be two plots generates:

\itemize{
  \item Fold change for the intersection of genomic features to correlated regions compared to whole genome.
  \item Intersection between genomic features and correlated regions
}
}
\value{
A list of two matrix:

\itemize{
  \item a matrix with how much each genomic feature that is overlapped by correlated regions
  \item a matrix with total number of base pairs of each genomic feature that is overlapped by correlated regions.
}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
