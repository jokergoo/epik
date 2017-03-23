\name{chromatin_states_transition_chord_diagram}
\alias{chromatin_states_transition_chord_diagram}
\title{
Chord diagram for visualizing chromatin states transitions
}
\description{
Chord diagram for visualizing chromatin states transitions
}
\usage{
chromatin_states_transition_chord_diagram(mat, group_names = NULL, max_mat = mat,
    remove_unchanged_transition = TRUE, state_col = NULL, legend_position = NULL, ...)
}
\arguments{

  \item{mat}{the transition matrix. It should be a square matrix in which row names and column names should be all the same. If it is not, the function will try to re-format it. If the matrix is from \code{\link{make_transition_matrix_from_chromHMM}} with methylation dataset, there will be additional tracks showing methylation difference between two groups.}
  \item{group_names}{name for the two groups under comparison. You also add it afterwards by using \code{\link[graphics]{text}}}
  \item{max_mat}{if there are several transition matrix to be compared, set it to the matrix with maximum absolute and it will make scales of all matrix the same and comparable.}
  \item{remove_unchanged_transition}{whether to remove transitions that states are not changed (set the values in diagonal to 0)}
  \item{state_col}{color for states. It should be a vector of which names correspond to states.}
  \item{legend_position}{positions of legends. Possible values are "bottomleft", "bottomright", "topright" and "topleft". If the value is specified as vector with length larger than two, the legend will be split into several parts. Set the value to \code{NULL} to suppress legends.}
  \item{...}{pass to \code{\link[circlize]{chordDiagram}}}

}
\details{
Rows of \code{mat} locate at the bottom of the circle by default. You can transpose the matrix to move rows to the top of the circle.

The chord diagram visualizes how much chromatin states change. In the diagram, width of each link represents the total
width of segments in a certain chromatin state in group 1 that transite to other chromatin state in group 2. The width of 
each grid represents total width of segments in a certain chromatin in group 1 that transite to all states in group 2.

If methylation dataset is provided when making the transistion matrix by using \code{\link{make_transition_matrix_from_chromHMM}},
there will be extra tracks on the outside of the circlie to represenst the mean methylation difference in two groups.

Chord diagram is implemented in base graphic system, which means, you can add titles or other graphics by base graphic 
functions (e.g. \code{\link[graphics]{title}}, \code{\link[graphics]{text}}, ...)

If you want to adjust order of states in the chord diagram, directly change row and column order of the matrix. If the matrix
is from \code{\link{make_transition_matrix_from_chromHMM}}, use \code{\link{state_names}} to retrieved or modify chromatin state names.
}
\value{
No value is returned.
}
\seealso{
\code{\link{make_transition_matrix_from_chromHMM}} which generates transition matrix directly from chromHMM results.
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
