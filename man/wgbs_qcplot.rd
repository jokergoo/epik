\name{wgbs_qcplot}
\alias{wgbs_qcplot}
\title{
Basic qc plot for bisulfite sequencing data
}
\description{
Basic qc plot for bisulfite sequencing data
}
\usage{
wgbs_qcplot(sample_id, chromosome = paste0("chr", 1:22))
}
\arguments{

  \item{sample_id}{a vector of sample ids}
  \item{chromosome}{a vector of chromosome names}

}
\details{
For each sample id, it will produce five plots:

\itemize{
  \item mean/median CpG coverage per chromosome
  \item histogram of CpG coverage
  \item methylation per chromosome 
  \item histogram of methylation
  \item mean Methylation for each CpG coverage 
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
