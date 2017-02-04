\name{cr_add_fdr_column}
\alias{cr_add_fdr_column}
\title{
Calcualte FDRs for CRs
}
\description{
Calcualte FDRs for CRs
}
\usage{
cr_add_fdr_column(cr, fdr_method = "BH")
}
\arguments{

  \item{cr}{correlated regions }
  \item{fdr_method}{method to calculate FDR}

}
\details{
Since correlated region detection is per-chromosome, after merging CRs from all chromosomes, FDR
can be calcualted based on \code{corr_p} and/or \code{meth_anno} column.

Please note, FDRs are calculated for negative CRs and positive CRs separatedly.
}
\value{
CR with two/one columns (\code{corr_fdr}, \code{meth_anova_fdr})
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
