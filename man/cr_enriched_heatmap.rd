\name{cr_enriched_heatmap}
\alias{cr_enriched_heatmap}
\title{
Visualize landscape of genome-wide correlations
}
\description{
Visualize landscape of genome-wide correlations
}
\usage{
cr_enriched_heatmap(cr, txdb, expr, expr_ha)
}
\arguments{

  \item{cr}{correlated regions}
  \item{txdb}{transcriptome annotation which was used in \code{\link{correlated_regions}}}
  \item{expr}{expression matrix which was used in \code{\link{correlated_regions}}}
  \item{expr_ha}{a \code{\link[ComplexHeatmap]{HeatmapAnnotation}} object for the expression heatmap}

}
\details{
The landscape of genome-wide correlations is visualized by a list of heatmaps.
Each row corresponds to a single gene:

\itemize{
  \item an enriched heatmap in which correlation signals are normalized at gene bodies
  \item a point plot showing the length of genes
  \item a heatmap of gene expression
  \item a heatmap showing the mean methylation in the extended gene regions.
  \item a heatmap showing the methylation difference in the extended gene regions.
}

K-means clustering with four groups is applied on the correlation matrix which has been normalized.
The four row subclusters are ordered by mean correlation. So basically, the four groups correspond to
negative gene body correlation, weak negative gene body correlation, weak positive gene body correlation
and positive gene body correlation. For each subcluster, rows are clustered by the mean methylation matrix.

There is another plot which shows quantitative statistics in the four groups:

\itemize{
  \item mean gene body correlation
  \item gene length
  \item expression difference
  \item methylation difference
}
}
\value{
An updated \code{cr} that includes the partitioning, this information is important for many downstream analysis.
Users should update it by \code{cr = cr_enriched_heatmap(cr, ...)}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
