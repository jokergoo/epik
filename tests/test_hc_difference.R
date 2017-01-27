source("/home/guz/project/development/epik/tests/test_head.R")
library(memoise)
source("/home/guz/project/development/epik/R/hilbert_curve_difference.R")
source("/home/guz/project/development/epik/R/methylation_genomic_features.R")
source("/home/guz/project/development/epik/R/common_utils.R")
source("/home/guz/project/development/epik/R/genomic_region_correlation.R")

library(HilbertCurve)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(GenomicFeatures)
library(circlize)

gr_meth = hilbert_curve_methylation_difference(subgroup = subgroup, comparison = c("group1", "group2"),
	chromosome = c("chr21", "chr22"), type = "none")


MARKS = c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3")

gr_list = lapply(MARKS, function(mk) hilbert_curve_chipseq_difference(mk, subgroup = subgroup, comparison = c("group1", "group2"),
	chromosome = c("chr21", "chr22"), type = "none"))
names(gr_list) = MARKS

gr_meth$diff = gr_meth$mean_group1 - gr_meth$mean_group2
gr_list = lapply(gr_list, function(gr) {
	gr$diff = gr$mean_group1 - gr$mean_group2
	gr
})

general_chipseq_association_to_methylation(gr_list, gr_meth)

pdf(width = 12, height = 10)
general_chipseq_association(gr_list)
general_chipseq_association(gr_list, q = seq(0.1, 0.9, by = 0.1))
dev.off()


# test average_in_window
gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 11, 21, 31, 41, 51), c(10, 20, 30, 40, 50, 60)))
gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(6, 16, 26), c(10, 20, 30)))
v = 1:3
average_in_window(gr1, gr2, v)
