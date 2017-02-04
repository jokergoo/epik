\name{make_transition_matrix_from_chromHMM}
\alias{make_transition_matrix_from_chromHMM}
\title{
Generate transition matrix from chromHMM results
}
\description{
Generate transition matrix from chromHMM results
}
\usage{
make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2, window = NULL,
    min_1 = floor(length(gr_list_1)/2), min_2 = floor(length(gr_list_2)/2), methylation_diff = 0,
    chromosome = paste0("chr", 1:22))
}
\arguments{

  \item{gr_list_1}{a list of \code{\link[GenomicRanges]{GRanges}} objects which contain chromatin states in group 1. The first column in meta columns should be the states. Be careful when importing bed files to  \code{\link[GenomicRanges]{GRanges}} objects (start positions in bed files are 0-based while 1-based in \code{GRanges} objects.}
  \item{gr_list_2}{a list of \code{\link[GenomicRanges]{GRanges}} objects which contains chromatin states in group 2.}
  \item{window}{window size which was used to do chromHMM states segmentation If it is not specified, the greatest common divisor of the width of all regions is used.}
  \item{min_1}{If there are multiple samples in the group, it is possible that a segment has more than one states asigned to it. If the recurrency of each state is relatively low, it means there is no one dominant state for this segment and it should  be removed. This argument controls the minimal value for the recurrency of states in a given segment.}
  \item{min_2}{same as \code{min_1}, but for samples in group 2.}
  \item{methylation_diff}{If methylation dataset is provided, the segments for which the methylation difference between two groups is less than this value are removed.}
  \item{chromosome}{subset of chromosomes}

}
\details{
The whole genome is segmentated by size of \code{window} and states with highest occurence among samples are assigned to segments.

To make the function run successfully, number of segments (after binned by \code{window}) in all samples 
should be all the same and there should not be gaps between segments
}
\value{
A transition matrix in which values represent total width of segments that transite from one state to the other in the two groups. Rows correspond
to group 1 and columns correspond to group 2.

If methylation dataset is provided, the mean methylation for each state in each group is attached, which will be used to calculate
mean methylation difference in \code{\link{chromatin_states_transition_chord_diagram}}.
}
\seealso{
The matrix can be sent to \code{\link{chromatin_states_transition_chord_diagram}} to visualize.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
set.seed(123)
gr_list_1 = lapply(1:5, function(i) {
	pos = sort(c(0, sample(1:9999, 99), 10000))*200
	GRanges(seqnames = "chr1", ranges = IRanges(pos[-101] + 1, pos[-1]), 
	    states = paste0("state_", sample(1:9, 100, replace = TRUE)))
})
gr_list_2 = lapply(1:5, function(i) {
	pos = sort(c(0, sample(1:9999, 99), 10000))*200
	GRanges(seqnames = "chr1", ranges = IRanges(pos[-101] + 1, pos[-1]), 
	    states = paste0("state_", sample(1:9, 100, replace = TRUE)))
})
mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2)
}
