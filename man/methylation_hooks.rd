\name{methylation_hooks}
\alias{methylation_hooks}
\title{
Read methylation dataset
}
\description{
Read methylation dataset
}
\usage{
methylation_hooks(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE)
}
\arguments{

  \item{...}{please ignore, see 'details' section.}
  \item{RESET}{please ignore}
  \item{READ.ONLY}{please ignore}
  \item{LOCAL}{please ignore}

}
\details{
Methylation dataset from whole genome bisulfite sequencing is always huge and it does not
make sense to read them all into the memory. Normally, the methylation dataset is stored
by chromosome and this hook function can be set to read methylation data in a per-chromosome
manner.

Generally, for methylation dataset, there are methylation rate, CpG coverage and genomic positions
for CpG sites. Sometimes there is also smoothed methylation rate. All these datasets can be set
by defining a proper \code{methylation_hooks$get_by_chr}. The value for \code{methylation_hooks$get_by_chr}
is a function with only one argument which is the chromosome name. This function defined how to
read methylation dataset for a single chromosome. The function must return a list which contains
following mandatory elements:

\describe{
  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object which contains genomic positions for CpG sites.}
  \item{meth}{a matrix which contains methylation rate. This will be the main methylation dataset the epik package uses, so it should be smoothed methylation rate if the CpG coverage is not high. Note, this matrix must have column names which is sample names and will be used to match other datasets (e.g. RNASeq)}
  \item{cov}{a matrix which contains CpG coverage.}
}

It can also contain some optional elements and they are not needed for the core analysis:

\describe{
  \item{raw}{a matrix which contains unsmoothed methylation rate (or the original methylation rate calculatd as the fraction of methylated CpG in a CpG site)}
}

Note each row in above datasets should correspond to the same CpG site.

After \code{methylation_hooks$get_by_chr} is properly set, the "current chromosome" for the methylation dataset
can be set by \code{methylation_hooks$set_chr(chr)} where \code{chr} is the chromosome name you want to go.
After validating the dataset, following variables can be used directly:

\itemize{
  \item \code{methylation_hooks$gr}
  \item \code{methylation_hooks$meth}
  \item \code{methylation_hooks$sample_id}
  \item \code{methylation_hooks$cov}
  \item \code{methylation_hooks$raw} if available
}
}
\value{
Hook functions
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
