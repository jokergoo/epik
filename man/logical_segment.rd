\name{logical_segment}
\alias{logical_segment}
\title{
Segmentation by a logical vector
}
\description{
Segmentation by a logical vector
}
\usage{
logical_segment(l)
}
\arguments{

  \item{l}{a logical vector}

}
\details{
The logical vector will be segmented according to their values.
It returns intervals for continuous \code{\link{TRUE}} values.
}
\value{
A data frame in which the first column is the index of start sites in the original vector and
the second column is the index of end sites.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
l = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
logical_segment(l)
}
