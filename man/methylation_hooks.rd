\name{methylation_hooks}
\alias{methylation_hooks}
\title{
Hook functions to extract methylation data
}
\description{
Hook functions to extract methylation data
}
\usage{
methylation_hooks(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE)
}
\arguments{

  \item{...}{Arguments for the parameters, see "details" section}
  \item{RESET}{reset to default values}
  \item{READ.ONLY}{whether only return read-only options}
  \item{LOCAL}{switch local mode}

}
\details{
Methylation from whole genome bisulfite sequencing is always huge and it does not
make sense to read them all into the memory (imaging there are 20M CpG sites on single strand in human genome). 
This hook sets how to read the methylation
data and how to return methylation data (e.g. CpG coverage, methylation rate...).

All downstream functions which analyze methylation data needs this hook to be already set.

There are following hooks:

\describe{
  \item{get_data}{how to get the object which contains methylation data. The function accepts a single chromosome name and  returns an object which is used as the first argument in other hook functions}
  \item{meth}{how to extract methylation rate. The function should have three arguments: the object returned from \code{set()}, index of rows and index of columns. Normally, the first argument (\code{obj}) can be ignored when calling this hook. Note the methylation matrix should have column names. The methylation rate should between 0 and 1.}
  \item{raw}{how to extract raw methylation value, same setting as \code{meth}. This hook is optional.}
  \item{site}{the function should return a vector of positions of CpG sites}
  \item{coverage}{how to extract CpG coverage, same setting as \code{meth}.}
  \item{GRanges}{howt to extract CpG sites as a \code{\link[GenomicRanges]{GRanges}} object.}
}

Following two hooks can be used if above hooks are set:

\describe{
  \item{set}{select chromosome as the current chromosome}
  \item{sample_id}{a vector of sample ids which contains methylation data}
}

Note: positions of CpG sites in a chromosome should be sorted.
}
\value{
Hook functions
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
NULL
\dontrun{
# following are examples of setting `methylation_hooks`
methylation_hooks$set = function(chr) {

    if(!is.null(methylation_hooks$obj)) {
        if(attr(methylation_hooks$obj, "chr") == chr) {
            qqcat("[@{chr}] @{chr} is already set.\n")
            return(invisible(NULL))
        }
    }
    obj = readRDS(GetoptLong::qq("path_to/@{chr}_methylation_data.rds"))
    attr(obj, "chr") = chr
    methylation_hooks$obj = obj

    return(invisible(NULL))
}

methylation_hooks$meth = function(obj = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

   if(is.null(row_index) && is.null(col_index)) {
        obj$meth[, , drop = FALSE]
    } else if(is.null(row_index)) {
        obj$meth[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        obj$meth[row_index, , drop = FALSE]
    } else {
        obj$meth[row_index, col_index, drop = FALSE]
    }

}

methylation_hooks$raw = function(obj = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

   if(is.null(row_index) && is.null(col_index)) {
        obj$meth[, , drop = FALSE]
    } else if(is.null(row_index)) {
        obj$meth[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        obj$meth[row_index, , drop = FALSE]
    } else {
        obj$meth[row_index, col_index, drop = FALSE]
    }
}

methylation_hooks$site = function(obj = methylation_hooks$obj, index = NULL) {
    if(is.null(index))
        start(obj$gr)
    else start(obj$gr[index])
}

methylation_hooks$GRanges = function(obj = methylation_hooks$obj) {
    obj$gr
}

methylation_hooks$coverage = function(obj = methylation_hooks$obj,
    row_index = NULL, col_index = NULL) {

    if(is.null(row_index) && is.null(col_index)) {
        obj$cov[, , drop = FALSE]
    } else if(is.null(row_index)) {
        obj$cov[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        obj$cov[row_index, , drop = FALSE]
    } else {
        obj$cov[row_index, col_index, drop = FALSE]
    }
}
 
}

}
