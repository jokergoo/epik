\name{cr_gviz}
\alias{cr_gviz}
\title{
Customized Gviz plot for a single gene
}
\description{
Customized Gviz plot for a single gene
}
\usage{
cr_gviz(sig_cr, gi, expr, txdb, species = "hg19", gf_list = NULL,
    hm_list = NULL, title = gi)
}
\arguments{

  \item{cr}{correlated regions which show significant correlations}
  \item{gi}{gene id}
  \item{expr}{the expression matrix which was used in \code{\link{correlation_regions}}}
  \item{txdb}{the transcriptome annotation which was used in \code{\link{correlation_regions}}}
  \item{species}{species}
  \item{gf_list}{a list of \code{\link[GenomicRanges]{GRanges}} objects which contains additional genomic annotations}
  \item{hm_list}{a list of \code{\link[GenomicRanges]{GRanges}} objects which contains histome modification peaks. The value is a two-layer list. The first layer is histome modifications and the second layer is the peaks in each sample which has current histome modification data. Name of the first layer is histome mark name and the name of the second layer is sample ID.}
  \item{title}{title of the plot}

}
\details{
There are following Gviz tracks:

\itemize{
  \item gene models. Multiple transcripts will also be plotted.
  \item correlation between methylation and expression
  \item heatmap for methylation
  \item significant correlated regions
  \item CpG density
  \item annotation to other genomic features, if provided
  \item histome modification data, if provided
}
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
