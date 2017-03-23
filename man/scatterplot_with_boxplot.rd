\name{scatterplot_with_boxplot}
\alias{scatterplot_with_boxplot}
\title{
Scatterplot with boxplots on both sides
}
\description{
Scatterplot with boxplots on both sides
}
\usage{
scatterplot_with_boxplot(x, y, subgroup = rep("unknown", length(x)),
    subgroup_col = structure(seq_along(levels(subgroup)), names = levels(subgroup)),
    main = NULL, xlab = NULL, ylab = NULL, xlim = range(x), ylim = range(y), text_list = NULL)
}
\arguments{

  \item{x}{values on x-axis}
  \item{y}{values on y-axis}
  \item{subgroup}{groups of data points}
  \item{subgroup_col}{colors for groups}
  \item{main}{title for the plot}
  \item{xlab}{labels on x-axis}
  \item{ylab}{labels on y-axis}
  \item{xlim}{range on x-axis}
  \item{ylim}{range on y-axis}
  \item{text_list}{additional text which is a named vector or list (if the text is mixed with character and numbers)}

}
\details{
Boxplots are on the left and bottom, the scatter plot is on the top right.
}
\value{
No value is returned
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
scatterplot_with_boxplot(rnorm(40), rnorm(40), subgroup = sample(letters[1:2], 40, replace = TRUE))
}
