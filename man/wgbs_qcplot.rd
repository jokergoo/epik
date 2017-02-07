\name{wgbs_qcplot}
\alias{wgbs_qcplot}
\title{
Basic qc plot for distribution of methylation and CpG coverage
}
\description{
Basic qc plot for distribution of methylation and CpG coverage
}
\usage{
wgbs_qcplot(sample_id, chromosome = paste0("chr", 1:22), background = NULL)
}
\arguments{

  \item{sample_id}{a vector of sample ids. You can generate plots for a list of samples in a same time while is faster than make it one by one.}
  \item{chromosome}{a vector of chromosome names}
  \item{background}{background regions where the CpG sites will only be looked into}

}
\details{
For each sample id, it will produce five plots:

\itemize{
  \item mean/median CpG coverage per chromosome, red area corresponds to the 25th and 75th percential.
  \item histogram of CpG coverage
  \item mean/median methylation per chromosome, red area corresponds to the 25th and 75th percential.
  \item histogram of methylation
  \item mean/median Methylation at each CpG coverage , red area corresponds to the 25th and 75th percential at each CpG coverage.
}
}
\value{
A list of corresponding statistics
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
