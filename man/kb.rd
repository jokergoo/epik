\name{kb}
\alias{kb}
\title{
Mark that the numbers represent number of kilo bases
}
\description{
Mark that the numbers represent number of kilo bases
}
\usage{
kb(x)
}
\arguments{

  \item{x}{a numeric vector.}

}
\details{
The input values are multiplied by 1000 and send to \code{\link{bp}}.
}
\value{
A numeric vector measured in bp
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
kb(10)
kb(10.01)
}
