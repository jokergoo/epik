\name{mb}
\alias{mb}
\title{
Add unit "mb" to the number
}
\description{
Add unit "mb" to the number
}
\usage{
mb(x)
}
\arguments{

  \item{x}{a numeric vector.}

}
\details{
The input values are multiplied by 1000000 and send to \code{\link{bp}}.
}
\value{
A numeric vector measured in bp
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
mb(10)
mb(10.01)
}
