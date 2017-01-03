source("test_head.R")

gr = hilbert_curve_methylation_difference(subgroup = subgroup, comparison = c("group1", "group2"),
	chromosome = c("chr21", "chr22"))

gr = hilbert_curve_chipseq_difference("H3K4me1", subgroup = subgroup, comparison = c("group1", "group2"),
	chromosome = c("chr21", "chr22"))
