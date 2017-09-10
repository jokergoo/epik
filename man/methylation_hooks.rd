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
  \item{RESET}{remove all hooks}
  \item{READ.ONLY}{please ignore}
  \item{LOCAL}{please ignore}

}
\details{
Methylation dataset from whole genome bisulfite sequencing is always huge and it does not
make sense to read them all into the memory. Normally, the methylation dataset is stored
by chromosome and this hook function can be set to read methylation data in a per-chromosome
manner. In the package, there are many functions use it internally to read methylation datasets.

Generally, for methylation dataset, there are methylation rate (ranging from 0 to 1), CpG coverage and genomic positions
for CpG sites. Sometimes there is also smoothed methylation rate. All these datasets can be set
by defining a proper \code{methylation_hooks$get_by_chr}. The value for \code{methylation_hooks$get_by_chr}
is a function with only one argument which is the chromosome name. This function defines how to
read methylation dataset for a single chromosome. The function must return a list which contains
following mandatory elements:

\describe{
  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object which contains genomic positions for CpG sites. Positions should be sorted.}
  \item{meth}{a matrix which contains methylation rate. This will be the main methylation dataset the epik package uses, so it should be smoothed methylation rate if the CpG coverage is not high. Note, this matrix must have column names which is sample names and will be used to match other datasets (e.g. RNASeq)}
  \item{cov}{a matrix which contains CpG coverage.}
}

It can also contain some optional elements and they are not needed for the core analysis:

\describe{
  \item{raw}{a matrix which contains unsmoothed methylation rate (or the original methylation rate calculatd as the fraction of methylated CpG in a CpG site)}
}

Note each row in above datasets should correspond to the same CpG site.

In following example code, assume the methylation data has been processed by bsseq package and saved as
\code{path/bsseq_$chr.rds}, then the definition of \code{methylation_hooks$get_by_chr} is:

  \preformatted{
  methylation_hooks$get_by_chr = function(chr) \{
      obj = readRDS(paste0("path/bsseq_", chr, ".rds"))
      lt = list(gr   = granges(obj),
                raw  = getMeth(obj, type = "raw"),
                cov  = getCoverage(obj, type = "Cov"),
                meth = getMeth(obj, type = "smooth")
      return(lt)
  \}  }

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

\code{methylation_hooks$set_chr(chr)} tries to reload the data only when the current chromosome changes.
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
