\name{chipseq_hooks}
\alias{chipseq_hooks}
\title{
Hook functions to extract ChIP-Seq peak data
}
\description{
Hook functions to extract ChIP-Seq peak data
}
\usage{
chipseq_hooks(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE)
}
\arguments{

  \item{...}{Arguments for the parameters, see "details" section}
  \item{RESET}{reset to default values}
  \item{READ.ONLY}{whether only return read-only options}
  \item{LOCAL}{switch local mode}

}
\details{
This hook defines how to get sample ids for a specific marks and how to get peak regions
by given mark type and sample id:

\describe{
  \item{sample_id}{how to extract sample ids. The argument is the name of the mark and it returns a vector of sample ids.}
  \item{peak}{how to get peak regions. The argument is the name of the mark and a single sample id. The function returns a \code{\link[GenomicRanges]{GRanges}} object.}
  \item{chromHMM}{how to get chromHMM data. The argument is a single sample id and the function returns a \code{\link[GenomicRanges]{GRanges}} object. Note the first column in meta columns should be the chromatin states.}
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
