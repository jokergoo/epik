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

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}} objects which show signal difference in two groups. Each \code{\link[GenomicRanges]{GRanges}} object must have a \code{diff} column. The object is usually from \code{\link{hilbert_curve_chipseq_difference}}.}
  \item{gr_meth}{a \code{\link[GenomicRanges]{GRanges}} object which shows methylation difference in two groups. The object is usually from \code{\link{hilbert_curve_methylation_difference}}}

}
\details{
For each histone mark, the distribution of methylation difference in regions which show
high histome modification signal difference is visualized.
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
