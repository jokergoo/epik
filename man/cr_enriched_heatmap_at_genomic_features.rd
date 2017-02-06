\name{cr_enriched_heatmap_at_genomic_features}
\alias{cr_enriched_heatmap_at_genomic_features}
\title{
Visualizing enrichment for epigenomic signals at TSS-CGIs
}
\description{
Visualizing enrichment for epigenomic signals at TSS-CGIs
}
\usage{
cr_enriched_heatmap_at_genomic_features(cr, txdb, expr, gf,
    fdr_cutoff = 0.05, meth_diff_cutoff = 0.1, marks = NULL, type = "neg", extend = 5000,
    min_reduce = 1, min_width = 1000, nearest_by = "tss", expr_annotation)
}
\arguments{

  \item{cr}{-cr}
  \item{txdb}{-txdb}
  \item{expr}{-expr}
  \item{gf}{-gf}
  \item{fdr_cutoff}{-fdr_cutoff}
  \item{meth_diff_cutoff}{-meth_diff_cutoff}
  \item{marks}{-marks}
  \item{type}{-type}
  \item{extend}{-extend}
  \item{min_reduce}{-min_reduce}
  \item{min_width}{-min_width}
  \item{nearest_by}{-nearest_by}
  \item{expr_annotation}{}

}
\details{

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
