\name{mat_dist}
\alias{mat_dist}
\title{
Visualize distribution of a matrix or a list
}
\description{
Visualize distribution of a matrix or a list
}
\usage{
mat_dist(x, subgroup = NULL, reorder_column = TRUE, od = if(is.matrix(x)) seq_len(ncol(x)) else seq_along(x),
    ha = if(is.null(subgroup)) NULL else HeatmapAnnotation(subgroup = subgroup, show_annotation_name = TRUE),
    type = c("densityHeatmap", "MDS"), title = title, ...)
}
\arguments{

  \item{x}{a matrix or a list. If it is a matrix, distribution in columns are visualized}
  \item{subgroup}{subgroup information}
  \item{reorder_column}{if it is true, samples are first ordered by subgroups and in each subgroup, samples are ordered by median values}
  \item{od}{order of columns. If \code{reorder_column} is set to \code{TRUE}, this argument is ignored.}
  \item{ha}{Annotations that are specified as a \code{\link[ComplexHeatmap]{HeatmapAnnotation-class}} object. The annotations will be put on top of the heatmap.}
  \item{type}{three types of plots are supported, see details}
  \item{title}{title for the plot}
  \item{...}{pass to \code{\link[ComplexHeatmap]{densityHeatmap}}}

}
\details{
Three types of plots for visualizing distributions are supported:

\describe{
  \item{densityHeatmap:}{density of distribution is visualized as heatmaps, use \code{\link[ComplexHeatmap]{densityHeatmap}}}
  \item{lineplot:}{distribution is visualized as normal line plot, use \code{\link[graphics]{matplot}}}
  \item{MDS:}{multiple dimension scaling by calling \code{\link[stats]{cmdscale}}}
}
}
\value{
Order of columns in the heatmap
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
