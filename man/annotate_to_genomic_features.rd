\name{annotate_to_genomic_features}
\alias{annotate_to_genomic_features}
\title{
Annotate to genomic features
}
\description{
Annotate to genomic features
}
\usage{
annotate_to_genomic_features(gr, genomic_features,
    name = NULL, type = c("percent", "number"), prefix = "overlap_to_", ...)
}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{genomic_features}{a single \code{\link[GenomicRanges]{GRanges}} object or a list of \code{\link[GenomicRanges]{GRanges}} objects}
  \item{name}{names for the genomic features if there is no name in \code{genomic_features}}
  \item{type}{How to calculate the values for the annotation. \code{number} means numbers of genomic features that each region in \code{gr} overlap; \code{percent} means the  percent of each region in \code{gr} that is overlapped by genomic features}
  \item{prefix}{prefix for names of the annotation columns}
  \item{...}{pass to \code{\link[GenomicRanges]{countOverlaps}} or \code{\link{percentOverlaps}}}

}
\details{
It adds new columns in \code{gr} which tell you how \code{gr} is overlaped by \code{genomic_features}.

Note for the annotation, strand information is ignored.
}
\value{
A \code{\link[GenomicRanges]{GRanges}} with additional columns of annotations.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
df1 = generateRandomBed(nr = 1000)
df2 = generateRandomBed(nr = 1000)
df3 = generateRandomBed(nr = 1000)
gr1 = GRanges(seqnames = df1[[1]], ranges = IRanges(df1[[2]], df1[[3]]))
gr2 = GRanges(seqnames = df2[[1]], ranges = IRanges(df2[[2]], df2[[3]]))
gr3 = GRanges(seqnames = df3[[1]], ranges = IRanges(df3[[2]], df3[[3]]))
annotate_to_genomic_features(gr1, list(gr2 = gr2, gr3 = gr3))
annotate_to_genomic_features(gr1, list(gr2 = gr2, gr3 = gr3), type = "number", prefix = "#")
}
