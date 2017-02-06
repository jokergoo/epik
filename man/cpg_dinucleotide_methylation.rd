\name{cpg_dinucleotide_methylation}
\alias{cpg_dinucleotide_methylation}
\title{
Merge methylation for CpG dinucleotide
}
\description{
Merge methylation for CpG dinucleotide
}
\usage{
cpg_dinucleotide_methylation(pos, meth, cov)
}
\arguments{

  \item{pos}{positions of CpG sites, must be sorted}
  \item{meth}{methylation matrix}
  \item{cov}{CpG coverage matrix}

}
\details{
normally methylation for the two Cs in a CpG dinucleotide is very similar. This function
helps to reduce the redundency of the methylation dataset.

For two Cs in a CpG dinucleotide, the merged methylation value is calculated by weighting
the CpG coverage of the two Cs. Also the merged coverage is also calculated by weighting
the coverage itself.
}
\value{
a list containg a \code{pos} vector, a \code{meth} matrix and a \code{cov} matrix
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
