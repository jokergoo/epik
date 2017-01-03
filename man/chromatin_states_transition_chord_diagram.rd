\name{chromatin_states_transition_chord_diagram}
\alias{chromatin_states_transition_chord_diagram}
\title{
Chord diagram for chromatin states transistion
}
\description{
Chord diagram for chromatin states transistion
}
\usage{
chromatin_states_transition_chord_diagram(mat, max_mat = mat,
    remove_unchanged_transition = TRUE, state_col = NULL, legend_position = NULL, ...)
}
\arguments{

  \item{mat}{the transition matrix. It should be a square matrix in which row names and column names are the same. If it is not, the function will try to re-format it.}
  \item{max_mat}{if there are several transition matrix to be compared, set it to the matrix with maximum absolute and it will make scales of all matrix the same and comparable.}
  \item{remove_unchanged_transition}{whether to remove transitions that states are not changed (set the values in diagonal to 0)}
  \item{state_col}{color for states. It should be a vector of which names correspond to states.}
  \item{legend_position}{positions of legends. Possible values are "bottomleft", "bottomright", "topright" and "topleft". If the value is specified as vector with length larger than two, the legend will be split into several parts. Set the value to \code{NULL} to suppress legends.}
  \item{...}{pass to \code{\link[circlize]{chordDiagram}}}

}
\details{
Rows of \code{mat} locate at the bottom of the circle by default.

The chord diagram visualizes how much chromatin states change. In the diagram, width of each link represents the total
width of regions in a certain chromatin state in group 1 that transite to other chromatin state in group 2. The width of 
each grid represents total width of regions in a certain chromatin in group 1 that transite to all states in group 2.

Chord diagram is implemented in base graphic system, which means, you can add titles or other graphics by base graphic 
functions (e.g. \code{\link[graphics]{title}}, \code{\link[graphics]{text}}, ...)

If you want to adjust order of states in the chord diagram, directly change row and column order of the matrix.
}
\value{
No value is returned.
}
\seealso{
\code{\link{make_transition_matrix_from_chromHMM}} which generates transition matrix from chromHMM results.
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
chromatin_states_transition_chord_diagram(mat, legend_position = "bottomleft")
}
