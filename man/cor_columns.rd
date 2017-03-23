\name{cor_columns}
\alias{cor_columns}
\title{
Number of columns which are highly correlated to other columns
}
\description{
Number of columns which are highly correlated to other columns
}
\usage{
cor_columns(x, abs_cutoff = 0.5, size = 1000, mc = 1, ...)
}
\arguments{

  \item{x}{a matrix, correlation is calculated by columns}
  \item{abs_cutoff}{cutoff of absolute correlation. It can be a numeric vector with more than one cutoffs.}
  \item{size}{size of blocks}
  \item{mc}{number of cores}
  \item{...}{pass to \code{\link[stats]{cor}}}

}
\details{
For each column, it looks for number of other columns which correlate with absolute correlation coefficient larger tham \code{abs_cutoff}.
The calculation involves pair-wise correlation of all columns in the matrix.
When number of columns is huge in the matrix, it is out of ability of R to store such long vector. This function
solves this problem by splitting the columns into k blocks and looks at each block sequentially or in parallel.

The code is partially adapted from \url{https://rmazing.wordpress.com/2013/02/22/bigcor-large-correlation-matrices-in-r/}
}
\value{
A matrix that represents how many other columns correlate to current column under the correlation cutoff.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
mat = matrix(rnorm(20000 * 10), ncol = 20000, nrow = 20)
cor_columns(mat, abs_cutoff = c(0.5, 0.6, 0.7))
}
NULL
}
