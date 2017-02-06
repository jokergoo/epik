\name{cr_scatterplot}
\alias{cr_scatterplot}
\title{
Scatter plot between methylation and expression in one correlated region
}
\description{
Scatter plot between methylation and expression in one correlated region
}
\usage{
cr_scatterplot(cr, expr, gi = NULL, text_column,
    xlab = "Methylation", ylab = "Expression")
}
\arguments{

  \item{cr}{correlated regions}
  \item{expr}{the expression matrix which was used in \code{\link{correlated_regions}}}
  \item{gi}{gene id}
  \item{text_column}{the column name in \code{cr} which will be plotted as text annotations.}
  \item{xlab}{xlab in the plot}
  \item{ylab}{ylab in the plot}

}
\details{
Scatterplot for all CRs corresponding to the gene will be made. If you want to make
a subset of CRs, directly subset \code{cr}.
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
