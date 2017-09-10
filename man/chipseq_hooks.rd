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
  \item{RESET}{remove all hooks}
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
  \item{sample_id}{This self-defined function returns a list of sample IDs given the name of a histone mark.}
  \item{peak}{This function should return a \code{\link[GenomicRanges]{GRanges}} object which are peaks for a given histone mark in a given sample. The \code{\link[GenomicRanges]{GRanges}} object should better have a meta column named "density" which is the density of the histone modification signals. (**Note when you want to take the histone modification signals as quatitative analysis, please make sure they are properly normalized between samples**)}
  \item{chromHMM}{This hook is optional. If chromatin segmentation by chromHMM is avaialble, this hook can be defined as a function which accepts sample ID as argument and returns a \code{\link[GenomicRanges]{GRanges}} object. The \code{\link[GenomicRanges]{GRanges}} object should have a meta column named "states" which is the chromatin states inferred by chromHMM.}
}

The \code{chipseq_hooks$peak()} must have two arguments \code{mark} and \code{sid} which are the name of the histone mark
and the sample id. There can also be more arguments such as chromosomes.

As an example, let's assume the peak files are stored in a format of \code{path/$sample_id/$mark.bed}, then we can define
hooks functions as:

  \preformatted{
  # here `qq` is from GetoptLong package which allows simple variable interpolation
  chipseq_hooks$sample_id = function(mark) \{
      peak_files = scan(pipe(qq("ls path/*/@\{mark\}.bed")), what = "character")
      sample_id = gsub("^path/(.*?)/.*$", "\\1", peak_files)
      return(sample_id)
  \}

  # here ... is important that the epik package will pass more arguments to it
  chipseq_hooks$peak = function(mark, sid, ...) \{
      peak_file = qq("path/@\{sid\}/@\{mark\}.bed")
      df = read.table(peak_file, sep = "\t", stringsAsFactors = FALSE)
      GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), density = df[[5]])
  \}  }

Normally \code{chipseq_hooks$peak()} are not directly used, it is usually used by \code{\link{get_peak_list}} to read peaks from all samples as a list.
You can also add more arguments when defining \code{chipseq_hooks$peak()} that these arguments can be passed from \code{\link{get_peak_list}} as well. 
For example, you can add chromosome name as the third argument that you do not need to read the full dataset at a time:

  \preformatted{
  # to make it simple, let's assume it only allows one single chromosome
  chipseq_hooks$peak = function(mark, sid, chr) \{
      peak_file = qq("path/@\{sid\}/@\{mark\}.bed")
      df = read.table(pipe(qq("awk '$1==\"@\{chr\}\"' @\{peak_file\}")), sep = "\t", stringsAsFactors = FALSE)
      GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), density = df[[5]])
  \}  }

then you can call \code{\link{get_peak_list}} as:

  \preformatted{
  get_peak_list(mark, chr = "chr1")  }

The \code{chipseq_hooks$chromHMM()} must have one argument \code{sid} which is the sample id, also there can be more arguments such as chromosomes.
The usage for the additional argumetns are same as \code{chipseq_hooks$peak()}.
}
\value{
Hook functions
}
\seealso{
\code{\link{get_peak_list}}, \code{\link{get_chromHMM_list}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
