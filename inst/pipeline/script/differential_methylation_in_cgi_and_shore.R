
suppressPackageStartupMessages(library(GetoptLong))
GetoptLong("config=s", "configuration R script")

library(epic)
load_config(config)
library(EnrichedHeatmap)

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

library(ComplexHeatmap)
ha = HeatmapAnnotation(df = SAMPLE, col = COLOR)

cat("making heatmaps for differential methylation.\n")
pdf(qq("@{OUTPUT_DIR}/heatmap_diff_methylation_1kb_window.pdf"), width = 14, height = 14)
gr_diff_genome = heatmap_diff_methylation_in_genomic_features(gr_complement, annotation = SAMPLE$class, 
    annotation_color = COLOR$class, title = "genome 1kb window (exclude cgi and shore)", ha = ha)
gr_diff_cgi = heatmap_diff_methylation_in_genomic_features(gr_cgi, annotation = SAMPLE$class, 
    annotation_color = COLOR$class, title = "cgi 1kb window", ha = ha, cutoff = 1, min_mean_range=0)
gr_diff_shore = heatmap_diff_methylation_in_genomic_features(gr_shore, annotation = SAMPLE$class, 
    annotation_color = COLOR$class, title = "shore 1kb window", ha = ha)
dev.off()


cat("correlate differentially methylated cgi/shores to other genomic features.\n")
MR_list = list(gr_diff_genome = gr_diff_genome,
               gr_diff_cgi = gr_diff_cgi,
               gr_diff_shore = gr_diff_shore)
res = genomic_regions_correlation(MR_list, GENOMIC_FEATURE_LIST[names(GENOMIC_FEATURE_LIST) != "cgi"], chromosome = CHROMOSOME, nperm = 0)
pdf(qq("@{OUTPUT_DIR}/genome_diff_1kb_window_correlation.pdf"), width = 6, height = 6)
ht = Heatmap(res$stat, name = "jaccard", column_title = qq("jaccard coefficient"), cluster_columns = FALSE)
draw(ht)
dev.off()

