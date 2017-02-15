\name{cr_add_subtype_specificity}
\alias{cr_add_subtype_specificity}
\title{
Add subtype specificity columns in cr
}
\description{
Add subtype specificity columns in cr
}
\usage{
cr_add_subtype_specificity(cr, cutoff = 0.05, suffix = "_ss")
}
\arguments{

  \item{cr}{correlated regions}
  \item{cutoff}{cutoff for p-values of ANOVA test}
  \item{suffix}{suffix of column names}

}
\details{
If \code{subgroup} is set in \code{\link{correlated_regions}}, this function can assign subtype specificity to each subtype.

We use following digits to represent subtype specificity: 1 is defined as the methylation is higher than all other subtypes and the difference is significant.
-1 is defined as the methylation is lower than all other subtypes and the difference is significant.
All the others are defined as 0.
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
