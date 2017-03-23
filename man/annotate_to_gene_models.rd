\name{annotate_to_gene_models}
\alias{annotate_to_gene_models}
\title{
Annotate genomic regions to gene models
}
\description{
Annotate genomic regions to gene models
}
\usage{
annotate_to_gene_models(gr, txdb, gene_model =c("gene", "tx"),
    species = "hg19", promoters_upstream = 2000, promoters_downstream = 200,
    annotation_type = c("percent", "number"),
    annotation_prefix = "overlap_to_")
}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{txdb}{a \code{\link[GenomicFeatures]{TxDb-class}} object.}
  \item{gene_model}{type of gene model, by transcript or by genes}
  \item{species}{We need this information to find out proper intergenic regions}
  \item{promoters_upstream}{length of upstream promoter from TSS, pass to \code{\link[GenomicRanges]{promoters}}}
  \item{promoters_downstream}{length of downstream promoter from TSS, pass o \code{\link[GenomicRanges]{promoters}}}
  \item{annotation_type}{pass to \code{\link{annotate_to_genomic_features}}}
  \item{annotation_prefix}{pass to \code{\link{annotate_to_genomic_features}}}

}
\value{
Following columns are attached to \code{gr}:

\describe{
  \item{nearest_*_tss}{the nearest tss (depending on \code{gene_model})}
  \item{dist_to_*_tss}{distance to the closest tss (depending on \code{gene_model})}
  \item{nearest_*}{the closest gene model (depending on \code{gene_model})}
  \item{dist_to_*}{distance to the closest gene model (depending on \code{gene_model})}
  \item{'prefix_to'_exon}{percent of the region which is covered by exons or number of exons overlapped to the region}
  \item{'prefix_to'_intron}{percent of the region which is covered by introns or number of introns overlapped to the region}
  \item{'prefix_to'_promoter}{percent of the region which is covered by promoters or number of promoters overlapped to the region}
  \item{'prefix_to'_intergenic}{percent of the region which is covered by intergenic regions or number of intergenic regions overlapped to the region}
  \item{'prefix_to'_fiveUTR}{percent of the region which is covered by 5'UTRs or number of 5'UTRs overlapped to the region}
  \item{'prefix_to'_threeUTR}{percent of the region which is covered by 3'UTRs or number of 3'UTRs overlapped to the region}
}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
