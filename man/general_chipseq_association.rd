\name{general_chipseq_association}
\alias{general_chipseq_association}
\title{
General association between histome marks
}
\description{
General association between histome marks
}
\usage{
general_chipseq_association(gr_list, q = 0.9)
}
\arguments{

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}} objects which show signal difference in two groups. Each \code{\link[GenomicRanges]{GRanges}} object must have a \code{diff} column. The object is usually from \code{\link{hilbert_curve_chipseq_difference}}.}
  \item{q}{quantile of difference}

}
\details{
For each pair of histome marks, the Jaccard coefficient for the regions which show higher difference
than \code{q} is calcualted and visualized.
}
\value{
no value is returned
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
