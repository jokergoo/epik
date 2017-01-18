
library(GetoptLong)

GetoptLong("rds=s", "path of rds",
	"output_prefix=s", "prefix")

library(GenomicRanges)
cr = readRDS(rds)

if(all(cr$mean_corr > 0)) {
	color = "red"
} else if(all(cr$mean_corr < 0)) {
	color = "darkgreen"
} else {
	color = "blue"
}

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))


gm = genes(TXDB)
gm = gm[seqnames(gm) %in% CHROMOSOME]
g = gm[gm$gene_id %in% unique(cr$gene_id)]
tss = promoters(g, upstream = 1, downstream = 0)
tes = GRanges(seqnames = seqnames(g), ranges = IRanges(
	start = ifelse(strand(g) == "+", end(g), start(g)),
	end = ifelse(strand(g) == "+", end(g), start(g))),
	strand = strand(g))
names(tes) = names(g)


mat_body_list = list()
mat_tss_list = list()
mat_tes_list = list()
mat_upstream_list = list()
mat_downstream_list = list()
for(chr in CHROMOSOME) {
	qqcat("normalize in @{chr}\n")
	mat_body_list[[chr]] = normalizeToMatrix(cr[seqnames(cr) == chr], g[seqnames(g) == chr], 
		mapping_column = "gene_id", extend = 0, k = 200, trim = 0, target_ratio = 1, empty_value = 0)
	mat_tss_list[[chr]] = normalizeToMatrix(cr[seqnames(cr) == chr], tss[seqnames(g) == chr], 
		mapping_column = "gene_id", extend = 5000, w = 50, trim = 0, empty_value = 0)
	mat_tes_list[[chr]] = normalizeToMatrix(cr[seqnames(cr) == chr], tes[seqnames(g) == chr], 
		mapping_column = "gene_id", extend = 5000, w = 50, trim = 0, empty_value = 0)
	mat_upstream_list[[chr]] = normalizeToMatrix(cr[seqnames(cr) == chr], tss[seqnames(g) == chr], 
		mapping_column = "gene_id", extend = c(50000, 0), w = 100, trim = 0, empty_value = 0)
	mat_downstream_list[[chr]] = normalizeToMatrix(cr[seqnames(cr) == chr], tes[seqnames(g) == chr], 
		mapping_column = "gene_id", extend = c(0, 50000), w = 100, trim = 0, empty_value = 0)
}
mat_body = do.call("rbind", mat_body_list)
for(at in setdiff(names(attributes(mat_body_list[[1]])), c("dim", "dimnames"))) {
	attr(mat_body, at) = attr(mat_body_list[[1]], at)
}
colnames(mat_body) = colnames(mat_body_list[[1]])

mat_tss = do.call("rbind", mat_tss_list)
for(at in setdiff(names(attributes(mat_tss_list[[1]])), c("dim", "dimnames"))) {
	attr(mat_tss, at) = attr(mat_tss_list[[1]], at)
}
colnames(mat_tss) = colnames(mat_tss_list[[1]])

mat_tes = do.call("rbind", mat_tes_list)
for(at in setdiff(names(attributes(mat_tes_list[[1]])), c("dim", "dimnames"))) {
	attr(mat_tes, at) = attr(mat_tes_list[[1]], at)
}
colnames(mat_tes) = colnames(mat_tes_list[[1]])

mat_upstream = do.call("rbind", mat_upstream_list)
for(at in setdiff(names(attributes(mat_upstream_list[[1]])), c("dim", "dimnames"))) {
	attr(mat_upstream, at) = attr(mat_upstream_list[[1]], at)
}
colnames(mat_upstream) = colnames(mat_upstream_list[[1]])

mat_downstream = do.call("rbind", mat_downstream_list)
for(at in setdiff(names(attributes(mat_downstream_list[[1]])), c("dim", "dimnames"))) {
	attr(mat_downstream, at) = attr(mat_downstream_list[[1]], at)
}
colnames(mat_downstream) = colnames(mat_downstream_list[[1]])


## filter out rows that have too many NAs in the body
w = width(g[rownames(mat_body)])
rainfall = rainfallTransform(as.data.frame(g[rownames(mat_body)]))

ylim_max = max(c(colMeans(mat_upstream, na.rm = TRUE),
	             colMeans(mat_tss, na.rm = TRUE),
	             colMeans(mat_body, na.rm = TRUE),
	             colMeans(mat_tes, na.rm = TRUE),
	             colMeans(mat_downstream, na.rm = TRUE)))

ht_list = EnrichedHeatmap(mat_upstream, name = "upstream", col = colorRamp2(c(0, 1), c("white", color)),
	use_raster = TRUE, raster_quality = 2, cluster_rows = TRUE, show_row_dend = FALSE,
	top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = color), ylim = c(0, ylim_max))), 
	top_annotation_height = unit(2, "cm"),
	show_heatmap_legend = TRUE, column_title = "upstream (50kb)", width = 1, axis_name = c("-50kb", "TSS"))
ht_list = ht_list + EnrichedHeatmap(mat_tss, name = "tss", col = colorRamp2(c(0, 1), c("white", color)),
	use_raster = TRUE, raster_quality = 2, cluster_rows = TRUE, show_row_dend = FALSE,
	top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = color), ylim = c(0, ylim_max))), 
	top_annotation_height = unit(2, "cm"),
	show_heatmap_legend = FALSE, column_title = "tss +-5kb", width = 1, axis_name = c("-5kb", "TSS", "5kb"))
ht_list = ht_list + EnrichedHeatmap(mat_body, name = "body", col = colorRamp2(c(0, 1), c("white", color)),
	use_raster = TRUE, raster_quality = 2, cluster_rows = TRUE, show_row_dend = FALSE,
	top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = color), ylim = c(0, ylim_max))), 
	top_annotation_height = unit(2, "cm"),
	show_heatmap_legend = FALSE, column_title = "gene body", width = 1, axis_name = c("TSS", "TES"))
ht_list = ht_list + EnrichedHeatmap(mat_tes, name = "tes", col = colorRamp2(c(0, 1), c("white", color)),
	use_raster = TRUE, raster_quality = 2, cluster_rows = TRUE, show_row_dend = FALSE,
	top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = color), ylim = c(0, ylim_max))), 
	top_annotation_height = unit(2, "cm"),
	show_heatmap_legend = FALSE, column_title = "tes +-5kb", width = 1, axis_name = c("-5kb", "TES", "5kb"))
ht_list = ht_list + EnrichedHeatmap(mat_downstream, name = "downstream", col = colorRamp2(c(0, 1), c("white", color)),
	use_raster = TRUE, raster_quality = 2, cluster_rows = TRUE, show_row_dend = FALSE,
	top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = color), ylim = c(0, ylim_max))), 
	top_annotation_height = unit(2, "cm"),
	show_heatmap_legend = FALSE, column_title = "downstream 50kb", width = 1, axis_name = c("TES", "50kb"))

ht_list = ht_list + Heatmap(w, name = "gene_length", col = colorRamp2(quantile(w, c(0, 0.95)), c("black", "white")), show_row_names = FALSE,
	width = unit(5, "mm")) + 
Heatmap(rainfall[[4]], name = "dist", col = colorRamp2(quantile(rainfall[[4]], c(0, 0.95)), c("black", "white")), show_row_names = FALSE,
	width = unit(5, "mm"))

row_order = hclust(dist(cbind(mat_upstream, mat_tss, mat_body, mat_tes, mat_downstream)))$order

pdf(qq("@{OUTPUT_DIR}/enrichedheatmap_all_@{output_prefix}.pdf"), width = 24, height = 12)
draw(ht_list, main_heatmap = "body", column_title = qq("@{nrow(mat_body)} genes"), cluster_rows = FALSE, row_order = row_order)
dev.off()


# all_files = dir("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/rds", 
# 	pattern = "all_(neg|pos)_cr_w6s3(.|.+\\d.)rds")

# name = gsub('.rds$', '', all_files)

# for(i in seq_along(all_files)) {
# 	cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/correlated_regions_general_heatmap.R --rds /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/rds/@{all_files[i]} --output_prefix @{name[i]}")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=10:00:00,mem=5G -N roadmap_heatmap_@{name[i]}' '@{cmd}'")
# 	system(cmd)
# }

