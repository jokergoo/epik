\name{chipseq_hooks}
\alias{chipseq_hooks}
\title{
Read ChIP-Seq dataset
}
\description{
Read ChIP-Seq dataset
}
\usage{
chipseq_hooks(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE)
}
\arguments{

  \item{...}{please ignore, see 'details' section.}
  \item{RESET}{please ignore}
  \item{READ.ONLY}{please ignore}
  \item{LOCAL}{please ignore}

}
\details{
Unlike methylation dataset which is always stored as matrix, ChIP-Seq dataset is stored
as a list of peak regions that each one corresponds to peaks in one sample. In many cases, 
there are ChIP-Seq datasets for multiple histone marks that each mark does not include all
samples sequenced in e.g. whole genome bisulfite sequencing or RNA-Seq, thus, to import
such type of flexible data format, users need to define following hook functions:

\describe{
  \item{sample_id}{This self-defined function returns a list of sample IDs given the name of a histome mark.}
  \item{peak}{This function should return a \code{\link[GenomicRanges]{GRanges}} object which are peaks for a given histome mark in a given sample. The \code{\link[GenomicRanges]{GRanges}} object should better have a meta column  which is the intensity of the histome modification signals.}
  \item{chromHMM}{This hook is optional. If chromatin segmentation by chromHMM is avaialble, this hoook can be defined as a function which accepts sample ID as argument and returns a \code{\link[GenomicRanges]{GRanges}} object. The \code{\link[GenomicRanges]{GRanges}} object should have a meta column named "states" which is the chromatin states inferred by chromHMM.}
}

The \code{chipseq_hooks$peak} must have two arguments \code{mark} and \code{sid} which are the name of the histone mark
and the sample id, there can be more arguments such as chromosomes.

The \code{chipseq_hooks$chromHMM} must have one argument \code{sid} which is the sample id, also there can be more arguments such as chromosomes.
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
