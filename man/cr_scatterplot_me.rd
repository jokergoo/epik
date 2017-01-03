\name{cr_scatterplot_me}
\alias{cr_scatterplot_me}
\title{
Scatter plot between methylation and expression in a correlated region
}
\description{
Scatter plot between methylation and expression in a correlated region
}
\usage{
cr_scatterplot_me(cr, expr, gi = NULL, text_column = NULL,
    xlab = "Methylation", ylab = "Expression")
}
\arguments{

  \item{cr}{correlated regions from \code{\link{correlated_regions}} or \code{\link{filter_correlated_regions}}}
  \item{expr}{the expression matrix which is same as in \code{\link{correlated_regions}}}
  \item{gi}{gene id}
  \item{text_column}{which column in \code{cr} should be put as annotation text in the plot}
  \item{xlab}{xlab in the plot}
  \item{ylab}{ylab in the plot}

}
\details{
Scatterplot for all CRs corresponding to the gene will be made.
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
