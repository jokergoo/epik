source("test_head.R")

methylation_hooks$set_chr("chr21")
sample_id = methylation_hooks$sample_id

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/"

library(GenomicFeatures)
txdb = loadDb(qq("@{BASE_DIR}/data/gen10.long.sqlite"))
genomic_features = list(
	gene = genes(txdb),
	exon = exons(txdb)
)
strand(genomic_features[[1]]) = "*"
strand(genomic_features[[2]]) = "*"

gr = get_mean_methylation_in_genomic_features(sample_id, genomic_features, chromosome = c("chr21", "chr22"))
gr2 = get_mean_methylation_in_genomic_features(sample_id, genomic_features$exon, chromosome = c("chr21", "chr22"))
heatmap_diff_methylation_in_genomic_features(gr[[1]], subgroup = subgroup, genomic_features = genomic_features)
