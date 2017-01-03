\name{cr_enriched_at_tss}
\alias{cr_enriched_at_tss}
\title{
Enrichment of cr around tss
}
\description{
Enrichment of cr around tss
}
\usage{
cr_enriched_at_tss(cr, txdb)
}
\arguments{

  \item{cr}{filtered correlated regions from \code{\link{filter_correlated_regions}}}
  \item{txdb}{a \code{GenomicFeatures::GRanges} object.}

}
\details{
It visualizes how neg cr and pos cr are distributed within c(-10kb, 20kb) of gene tss.
Also the distribution of transcript will also be added. Values on y-axis represent
how many crs/transcripts cover each position relative to the tss.
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
