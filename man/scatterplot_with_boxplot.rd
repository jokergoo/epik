\name{scatterplot_with_boxplot}
\alias{scatterplot_with_boxplot}
\title{
Scatterplot with boxplots on both sides
}
\description{
Scatterplot with boxplots on both sides
}
\usage{
scatterplot_with_boxplot(x, y, annotation = rep("unknown", length(x)),
    annotation_color = structure(seq_along(levels(annotation)), names = levels(annotation)),
    main = NULL, xlab = NULL, ylab = NULL, xlim = range(x), ylim = range(y), text_list = NULL)
}
\arguments{

  \item{x}{values on x-axis}
  \item{y}{values on y-axis}
  \item{annotation}{annotations which show groups of data points}
  \item{annotation_color}{colors for annotation}
  \item{main}{title for the plot}
  \item{xlab}{labels on x-axis}
  \item{ylab}{labels on y-axis}
  \item{xlim}{range on x-axis}
  \item{ylim}{range on y-axis}
  \item{text_list}{additional text which is a named vector or list (if the text is mixed with character and numbers)}

}
\details{
On the left and bottom, there are boxplots and on the top right, there is the scatter plot.
}
\value{
No value is returned
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
