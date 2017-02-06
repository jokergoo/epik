\name{cr_enrichedheatmap}
\alias{cr_enrichedheatmap}
\title{
Visualize landscape of genome-wide correlations
}
\description{
Visualize landscape of genome-wide correlations
}
\usage{
cr_enrichedheatmap(cr, txdb, expr, expr_annotation)
}
\arguments{

  \item{cr}{correlated regions}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{expr}{expressio matrix which was used in \code{\link{correlated_regions}}}
  \item{expr_annotation}{a \code{\link[ComplexHeatmap]{HeatmapAnnotation}} objects}

}
\details{
The landscape of genome-wide correlations is visualized by a list of heatmaps.
Each row corresponds to a single gene:

\itemize{
  \item an enriched heatmap in which correlation is normalized at gene bodies
  \item a point plot showing the length of genes
  \item a heatmap of gene expression
  \item an heatmap showing the mean methylation in the extended gene regions.
  \item an heatmap showing the methylation difference in the extended gene regions.
}

K-means clustering with four groups is applied on the correlation normalized matrix
and the four row subclusters are ordered by mean correlation.

There are also general statistic plots generated.
}
\value{
An updated \code{cr} that includes the partitioning.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
