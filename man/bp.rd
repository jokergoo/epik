\name{bp}
\alias{bp}
\title{
Mark the numbers to be number of base pairs
}
\description{
Mark the numbers to be number of base pairs
}
\usage{
bp(x)
}
\arguments{

  \item{x}{a numeric vector. It will be convert to integers by \code{\link[base]{as.integer}}.}

}
\details{
It adds a new \code{bp} class to the vector.
}
\value{
A same vector as input
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
bp(10)
bp(10.1)
}
