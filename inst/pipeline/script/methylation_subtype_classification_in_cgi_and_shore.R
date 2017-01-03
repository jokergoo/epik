# this script depends on `differential_methylation_in_cgi_and_shore.R`

suppressPackageStartupMessages(library(GetoptLong))
GetoptLong("config=s", "configuration R script",
           "n_class=i", "number of classes expected")

library(epic)
load_config(config)

library(ComplexHeatmap)
library(EnrichedHeatmap)

ha = HeatmapAnnotation(df = SAMPLE, col = COLOR)


if(!all(file.exists(c(qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi.rds"), qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi_shore.rds"), qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_neither_cgi_nor_shore.rds"))))) {
	sample_id = rownames(SAMPLE)
	n_sample = length(sample_id)


	chromInfo = getChromInfoFromUCSC(GENOME)
	chromInfo = chromInfo[chromInfo$chrom %in% CHROMOSOME, ]
	chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))

	cat("split cgi by 1kb window and calculate mean methylation in it.\n")
	cgi_1kb_window = makeWindows(GENOMIC_FEATURE_LIST$cgi, w = 1000, short.keep = TRUE)
	cgi_1kb_window = cgi_1kb_window[width(cgi_1kb_window) > 500]
	gr_cgi = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(cgi_1kb_window = cgi_1kb_window))[[1]]

	shore = GENOMIC_FEATURE_LIST$cgi_shore
	shore_1kb_window = makeWindows(shore, w = 1000, short.keep = TRUE)
	shore_1kb_window = shore_1kb_window[width(shore_1kb_window) > 500]
	gr_shore = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(shore_1kb_window = shore_1kb_window))[[1]]

	cat("split genome by 1kb window and calculate mean methylation in it.\n")
	complement_1kb_window = makeWindows(setdiff(chromGr, union(GENOMIC_FEATURE_LIST$cgi, GENOMIC_FEATURE_LIST$cgi_shore)), w = 1000)
	complement_1kb_window = complement_1kb_window[width(complement_1kb_window) > 500]
	gr_complement = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(complement_1kb_window = complement_1kb_window))[[1]]

	saveRDS(gr_cgi, file = qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi.rds"))
	saveRDS(gr_shore, file = qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi_shore.rds"))
	saveRDS(gr_complement, file = qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_neither_cgi_nor_shore.rds"))

}



cat("subtype classification based on cgi\n")
gr_cgi = readRDS(qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi.rds"))
pdf(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi.pdf"), width = 14, height = 14)
methylation_subtype_classfication(gr_cgi, n_class = n_class, pct_cutoff = 0.2, corr_cutoff = 0.5, k = 50, ha = ha)
dev.off()


cat("subtype classification based on cgi shores\n")
gr_shore = readRDS(qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi_shore.rds"))
pdf(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi_shore.pdf"), width = 14, height = 14)
methylation_subtype_classfication(gr_shore, n_class = n_class, pct_cutoff = 0.2, corr_cutoff = 0.6, k = 500, ha = ha)
dev.of()


cat("subtype classification based on other parts in genome\n")
gr_complement = readRDS(qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_neither_cgi_nor_shore.rds"))
pdf(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_neither_cgi_nor_shores.pdf"), width = 14, height = 14)
methylation_subtype_classfication(gr_complement, n_class = n_class, pct_cutoff = 0.02, corr_cutoff = 0.8, k = 1000, ha = ha)
dev.of()
