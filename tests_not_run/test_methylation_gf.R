source("/home/guz/project/development/epik/tests/test_head.R")
source("/home/guz/project/development/epik/R/methylation_genomic_features.R")

methylation_hooks$set_chr("chr21")
sample_id = methylation_hooks$sample_id

library(GenomicFeatures)
txdb = loadDb(qq("@{PROJECT_DIR}/data/gen10.long.sqlite"))
genomic_features = list(
	gene = genes(txdb),
	exon = exons(txdb)
)
strand(genomic_features[[1]]) = "*"
strand(genomic_features[[2]]) = "*"

gr_list = get_mean_methylation_in_genomic_features(sample_id, genomic_features, chromosome = c("chr21", "chr22"))

identical(ranges(gr_list[[1]]), ranges(genomic_features[[1]][seqnames(genomic_features[[1]]) %in% c("chr21", "chr22")]))

gr_list2 = get_mean_methylation_in_genomic_features(sample_id, genomic_features$exon, chromosome = c("chr21", "chr22"))
heatmap_diff_methylation_in_genomic_features(gr_list[[1]], subgroup = subgroup, genomic_features = genomic_features)
