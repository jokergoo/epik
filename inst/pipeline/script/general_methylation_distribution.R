
suppressPackageStartupMessages(library(GetoptLong))
GetoptLong("config=s", "configuration R script")

library(epic)
load_config(config)

library(ComplexHeatmap)

sample_id = rownames(SAMPLE)
n_sample = length(sample_id)

ha = HeatmapAnnotation(df = SAMPLE, col = COLOR)

cat("general methylation distribution whole genome-wide\n")
#pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution.pdf"), width = 0.16*n_sample + 2, height = 10)
pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution.pdf"), width = 10, height = 10)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, ha = ha)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, by_chr = TRUE, ha = ha)
dev.off()

cat("general methylation distribution in cpg islands\n")
#pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi.pdf"), width = 0.16*n_sample + 2, height = 10)
pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi.pdf"), width = 10, height = 10)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, background = GENOMIC_FEATURE_LIST$cgi, p = 0.01, ha = ha, 
	meth_range = c(0, 0.1))
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, by_chr = TRUE, background = GENOMIC_FEATURE_LIST$cgi, p = 0.01, ha = ha, 
	meth_range = c(0, 0.1))
dev.off()


cat("general methylation distribution in cgi shores\n")
#pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi_shores.pdf"), width = 0.16*n_sample + 2, height = 10)
pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi_shores.pdf"), width = 10, height = 10)
shore = GENOMIC_FEATURE_LIST$cgi_shore
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, background = shore, p = 0.01, ha = ha)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, by_chr = TRUE, background = shore, p = 0.01, ha = ha)
dev.off()

cat("general methylation distribution in neither cgi nor shores\n")
#pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution_neither_cgi_nor_shores.pdf"), width = 0.16*n_sample + 2, height = 10)
pdf(qq("@{OUTPUT_DIR}/general_methylation_distribution_neither_cgi_nor_shores.pdf"), width = 10, height = 10)
chromInfo = getChromInfoFromUCSC(GENOME)
chromInfo = chromInfo[chromInfo$chrom %in% CHROMOSOME, ]
chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
complement = setdiff(chromGr, union(GENOMIC_FEATURE_LIST$cgi, GENOMIC_FEATURE_LIST$cgi_shore))
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, background = complement, ha = ha)
global_methylation_distribution(sample_id = sample_id, annotation = SAMPLE$class, 
	annotation_color = COLOR$class, chromosome = CHROMOSOME, by_chr = TRUE, background = complement, ha = ha)
dev.off()

pdf(qq("@{OUTPUT_DIR}/basic_WGBS_qc.pdf"), width = 12, height = 8)
wgbs_qcplot(sample_id, chromosome = CHROMOSOME)
dev.off()

pdf(qq("@{OUTPUT_DIR}/basic_WGBS_qc_in_cgi.pdf"), width = 12, height = 8)
wgbs_qcplot(sample_id, chromosome = CHROMOSOME, background = GENOMIC_FEATURE_LIST$cgi)
dev.off()

pdf(qq("@{OUTPUT_DIR}/basic_WGBS_qc_in_cgi_shore.pdf"), width = 12, height = 8)
wgbs_qcplot(sample_id, chromosome = CHROMOSOME, background = GENOMIC_FEATURE_LIST$cgi_shore)
dev.off()

pdf(qq("@{OUTPUT_DIR}/basic_WGBS_qc_neither_cgi_nor_shores.pdf"), width = 12, height = 8)
wgbs_qcplot(sample_id, chromosome = CHROMOSOME, background = complement)
dev.off()
