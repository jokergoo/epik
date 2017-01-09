source("test_head.R")

gr_meth = hilbert_curve_methylation_difference(subgroup = subgroup, comparison = c("group1", "group2"),
	chromosome = c("chr21", "chr22"), type = "none")


MARKS = c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3")

gr_list = lapply(MARKS, function(mk) hilbert_curve_chipseq_difference(mk, subgroup = subgroup, comparison = c("group1", "group2"),
	chromosome = c("chr21", "chr22"), type = "none"))
names(gr_list) = MARKS

gr_meth$diff = gr_meth$mean_group1 - gr_meth$mean_group2
gr_list = lapply(gr_list, function(gr) {
	mean = (gr$mean_group1 + gr$mean_group2)/2
	s0 = quantile(mean[mean > 1e-6], 0.05)
	# gr$diff = (gr$mean_group1 - gr$mean_group2)/(mean + s0)
	gr$diff = gr$mean_group1 - gr$mean_group2
	gr
})

general_chipseq_association_to_methylation(gr_list, gr_meth)

pdf(width = 12, height = 8)
general_chipseq_association(gr_list)
general_chipseq_association(gr_list, q = seq(0.1, 0.9, by = 0.1))
dev.off()
