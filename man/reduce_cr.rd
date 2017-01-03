\name{reduce_cr}
\alias{reduce_cr}
\title{
Merge neighbouring cr regions
}
\description{
Merge neighbouring cr regions
}
\usage{
reduce_cr(cr, expr, txdb, max_gap = 1000, gap = 1.0, mc.cores = 1)
}
\arguments{

  \item{cr}{filtered correlated regions from \code{\link{filter_correlated_regions}}}
  \item{expr}{the expression matrix which is same as in \code{\link{correlated_regions}}}
  \item{txdb}{a \code{GenomicFeatures::GRanges} object.}
  \item{max_gap}{maximum gap for merging.}
  \item{gap}{gap for merging, a numeric value represents the ratio of width of itself and use \code{\link{bp}}, \code{\link{kb}} or \code{\link{mb}} to represent the number is absoltue base pairs. Pass to \code{\link{reduce2}}.}
  \item{mc.cores}{number of cores}

}
\details{
Since cr with positive correlation and negative correlation has different distribution patterns
in the genome, i.e. generally, pos_cr are long and CpG density in it is low while neg_cr is short,
clustered and has high CpG density, thus, pos cr and neg cr are reduced separatedly.

Even only look at e.g. pos cr, the pattern for the distribution in the genome is still different, 
thus, we recommend to merge crs by width itself while not by an absolute value for which you can set
\code{gap} by a numeric value.

Since original regions are merged and columns related to calculation of correlation will be dropped.

The merging is applied by gene, so it is still possible that regions associated with gene A overlap to gene B.
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
