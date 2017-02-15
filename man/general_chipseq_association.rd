\name{general_chipseq_association}
\alias{general_chipseq_association}
\title{
General association between histone modifications
}
\description{
General association between histone modifications
}
\usage{
general_chipseq_association(gr_list, q = 0.9)
}
\arguments{

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}} objects that each one corresponds to one histone modification Each \code{\link[GenomicRanges]{GRanges}} object must have a \code{diff} column which is the signal difference in two groups. The object is usually from \code{\link{hilbert_curve_chipseq_difference}}.}
  \item{q}{quantile of difference where Venn diagrams are added}

}
\details{
For every pair of histone marks, the Jaccard coefficient for the regions which show higher difference
than quantile \code{q} is calculated and visualized. There are three plots:

For each pair of histone marks,

\itemize{
  \item one mark having higher signals in group 1 and the other having higher signals in group 2, or vise versa
  \item both marks having higher signals in group 1
  \item both marks having higher signals in group 2
}
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
