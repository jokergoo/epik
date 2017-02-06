\name{binary_search}
\alias{binary_search}
\title{
Find intervals by binary search
}
\description{
Find intervals by binary search
}
\usage{
binary_search(breaks, search, left_index = TRUE)
}
\arguments{

  \item{breaks}{a non-decreasing integer vector}
  \item{search}{an integer vector}
  \item{left_index}{whether to use the index of left break or right break}

}
\value{
A vector of index.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
site = c(2, 5, 9, 10, 15, 20)
binary_search(site, c(1, 5, 12, 30), FALSE)
binary_search(site, c(1, 5, 12, 30), TRUE)
}
