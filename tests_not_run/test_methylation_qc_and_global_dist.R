library(circlize)
library(gtrellis)
source("/home/guz/project/development/epik/tests/test_head.R")
source("/home/guz/project/development/epik/R/methylation_qc_and_distribution.R")
library(EnrichedHeatmap)
library(ComplexHeatmap)

methylation_hooks$set_chr("chr21")
sample_id = methylation_hooks$sample_id

wgbs_qcplot(sample_id[1])
wgbs_qcplot(sample_id[1], background = CGI)
wgbs_qcplot(sample_id[1], background = CGI_SHORE)


methylation_gtrellis(sample_id[1], nrow = 3, compact = TRUE)

methylation_gtrellis_multiple_samples(sample_id, subgroup = subgroup, nrow = 3, compact = TRUE)

ha = HeatmapAnnotation(subgroup = subgroup, col = list(subgroup = c("group1" = "red", "group2" = "blue")))
methylation_global_distribution(sample_id, subgroup = subgroup, ha = ha)
methylation_global_distribution(sample_id, subgroup = subgroup, ha = ha, by_chr = TRUE)
methylation_global_distribution(sample_id, subgroup = subgroup, ha = ha, background = CGI)
methylation_global_distribution(sample_id, subgroup = subgroup, ha = ha, background = CGI, meth_range = c(0, 0.1), by_chr = TRUE)
bg = lapply(seq_along(sample_id), function(i) {
	df = generateRandomBed(nr = 1000)
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
})
methylation_global_distribution(sample_id, subgroup = subgroup, ha = ha, background = bg, chromosome = c("chr21", "chr22"))
