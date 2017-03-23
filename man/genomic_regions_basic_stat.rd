\name{genomic_regions_basic_stat}
\alias{genomic_regions_basic_stat}
\title{
Visualize basic statistics on genomic regions
}
\description{
Visualize basic statistics on genomic regions
}
\usage{
genomic_regions_basic_stat(gr_list, subgroup = NULL, subgroup_color = NULL,
    title = paste0("Basic statistics for genomic regions"), species = "hg19",
    type = c("proportion", "number", "median_width"),
    chromosome = paste0("chr", c(1:22, "X", "Y")), by_chr = FALSE)
}
\arguments{

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}}.}
  \item{subgroup}{a vector which contains groups of samples. If it has names which correspond to \code{gr_list}, the order of this vector is automatically adjusted.}
  \item{subgroup_color}{colors corresponding to subgroups}
  \item{title}{title of the plot}
  \item{species}{species, necessary if \code{type} is set to \code{proportion}.}
  \item{type}{type of statistics}
  \item{chromosome}{subset of chromosomes}
  \item{by_chr}{take all chromosomes as a whole or calculate statistics for every chromosome}

}
\details{
The function makes barplot to visualize different statistics in all samples.

For \code{type} settings:

\describe{
  \item{proportion}{proportion of total length of regions in the whole genome.}
  \item{number}{number of regions. Sometimes only looking at the number of regions gives biased estimation of amount of regions if the width of regions are very viarable.}
  \item{median_width}{median width of regions}
}
}
\value{
A data frame which contains statistics for each chromosome in each sample.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
set.seed(123)
gr_list = lapply(1:10, function(i) {
	df = generateRandomBed(1000)[1:sample(100:1000, 1), ]
	GRanges(df[[1]], ranges = IRanges(df[[2]], df[[3]]))
})
names(gr_list) = paste0("sample_", 1:10)
genomic_regions_basic_stat(gr_list)
genomic_regions_basic_stat(gr_list, subgroup = rep(letters[1:2], each = 5), 
    subgroup_color = c("a" = "red", "b" = "blue"))
genomic_regions_basic_stat(gr_list, subgroup = rep(letters[1:2], each = 5), 
    subgroup_color = c("a" = "red", "b" = "blue"), type = "number")
genomic_regions_basic_stat(gr_list, subgroup = rep(letters[1:2], each = 5), 
    subgroup_color = c("a" = "red", "b" = "blue"), type = "median_width")
genomic_regions_basic_stat(gr_list, subgroup = rep(letters[1:2], each = 5), 
    subgroup_color = c("a" = "red", "b" = "blue"), by_chr = TRUE)
}
