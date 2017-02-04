\name{kb}
\alias{kb}
\title{
Add unit "kb" to the number
}
\description{
Add unit "kb" to the number
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
