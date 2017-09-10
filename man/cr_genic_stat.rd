\name{cr_genic_stat}
\alias{cr_genic_stat}
\title{
Plot general statistics for the annotations to genes
}
\description{
Plot general statistics for the annotations to genes
}
\usage{
cr_genic_stat(cr_reduced, txdb)
}
\arguments{

  \item{cr_reduced}{significant correlated regions which have been reduced by \code{\link{cr_reduce}}}
  \item{txdb}{the transcriptome annotation which is same as the one used in \code{\link{correlated_regions}}}

}
\details{
There are five plots:

\itemize{
  \item mean methylation difference for negative/positive CRs in promoters/gene body/intergenic regions
  \item length of negative/positive CRs in promoters/gene body/intergenic regions
  \item number of genes for which promoters/gene body/intergenic regions are affected by negative/positive CRs
  \item sum of length of negative/positive CRs in promoters/gene body/intergenic regions
  \item Number of CpG per 1kb window for negative/positive CRs in promoters/gene body/intergenic regions
}

Here the promoter is defined as (-1000, 2000) of TSS.
}
\value{
No value is returned
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
