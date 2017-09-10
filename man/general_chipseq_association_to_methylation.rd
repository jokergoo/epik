\name{general_chipseq_association_to_methylation}
\alias{general_chipseq_association_to_methylation}
\title{
General association between histone modifications and methylations
}
\description{
General association between histone modifications and methylations
}
\usage{
general_chipseq_association_to_methylation(gr_list, gr_meth)
}
\arguments{

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}} objects that each one corresponds to one histone modification Each \code{\link[GenomicRanges]{GRanges}} object must have a \code{diff} column which is the signal difference in two groups. The object is usually from \code{\link{hilbert_curve_chipseq_difference}}.}
  \item{gr_meth}{a \code{\link[GenomicRanges]{GRanges}} object which shows methylation difference in two groups. The object is usually from \code{\link{hilbert_curve_methylation_difference}}}

}
\details{
Each histone mark corresponds to one panel. In each panel, there are two plots:

\itemize{
  \item distribution of methylation difference in regions where histone mark signals are larger than corresponding quantile.
  \item proportion of regions which show higher histone modification signal in group 1 and in group2.
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
