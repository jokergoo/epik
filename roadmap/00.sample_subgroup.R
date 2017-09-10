library(GetoptLong)

rerun = FALSE
GetoptLong("rerun!", "whether ignore cache rdata files")

library(epik)
library(epik.cmd)

PROJECT_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/epik_roadmap/"
initialize_project_directory(PROJECT_DIR)

## roadmap_configure_primary.R is same as roadmap_configure.R except the two samples (E053, E058)
## are not removed because in this script we want to show there are inconsistency in RNASeq
## classification and WGBS classficiation for these two samples.
source(qq("@{PROJECT_DIR}/scripts/configure/roadmap_configure_primary.R"))

sample_id = rownames(SAMPLE)
n_sample = length(sample_id)

library(ConsensusClusterPlus)
library(gtrellis)

####################################################
# cluster by expressoin

## top 2k genes with highest MAD values
od = order(rowMads(EXPR), decreasing = TRUE)[1:2000]
set.seed(123)
res_expr = ConsensusClusterPlus(EXPR[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_expr = res_expr[[2]]$consensusClass[sample_id])
COLOR$km_expr = structure(res_expr[[2]]$clrs[[3]], names = unique(res_expr[[2]]$consensusClass))

pdf(qq("@{PROJECT_DIR}/image/expression_group_classification.pdf"), width = 8, height = 10)

ht = Heatmap(EXPR[od, ], name = "expr", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_expr[[2]]$consensusTree,
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_expr = SAMPLE[sample_id, ]$km_expr, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_expr = COLOR$km_expr), show_annotation_name = TRUE),
	column_title = qq("@{length(od)} rows, based on MAD of expression"))
draw(ht)

# rainfall plot for the top 2k genes to show they distribute wide spread in the genome,
# thus can represent the whole genome
g = genes(TXDB)[rownames(EXPR[od, ])]
rainfall = rainfallTransform(as.data.frame(g))
den = genomicDensity(as.data.frame(g), window.size = 5e6)
gtrellis_layout(n_track = 2, nrow = 6, compact = TRUE, category = CHROMOSOME,
    track_axis = TRUE, track_height = c(3, 1),
    track_ylim = c(0, ceiling(max(log10(rainfall[[4]]))), 0, max(den[[4]])),
    track_ylab = c("log10(inter_dist)", "density"),
    add_name_track = TRUE, add_ideogram_track = TRUE, title = "rainfall plot for top 2000 genes with highest MAD of expression")
add_points_track(rainfall, log10(rainfall[[4]]), pch = 16, size = unit(1, "mm"), gp = gpar(col = "#FF000080"))
add_lines_track(den, den[[4]], area = TRUE, gp = gpar(fill = "pink"))
dev.off()

## classification based on WGBS in CGI, CGI-shore and complement regions
if(!all(file.exists(c(qq("@{PROJECT_DIR}/rds/mean_meth_1kb_cgi.rds"), 
					  qq("@{PROJECT_DIR}/rds/mean_meth_1kb_cgi_shore.rds"), 
					  qq("@{PROJECT_DIR}/rds/mean_meth_1kb_neither_cgi_nor_shore.rds")))) && !rerun) {
	
	chromInfo = getChromInfoFromUCSC(GENOME)
	chromInfo = chromInfo[chromInfo$chrom %in% CHROMOSOME, ]
	chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))

	# split by 1kb window and calculate the mean methylation
	cat("split cgi by 1kb window and calculate mean methylation in it.\n")
	cgi_1kb_window = makeWindows(CGI, w = 1000, short.keep = TRUE)
	cgi_1kb_window = cgi_1kb_window[width(cgi_1kb_window) > 500]
	gr_cgi = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(cgi_1kb_window = cgi_1kb_window), filter_fun = function(s) TRUE)[[1]]

	shore = CGI_SHORE
	shore_1kb_window = makeWindows(shore, w = 1000, short.keep = TRUE)
	shore_1kb_window = shore_1kb_window[width(shore_1kb_window) > 500]
	gr_shore = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(shore_1kb_window = shore_1kb_window), filter_fun = function(s) TRUE)[[1]]

	cat("split genome by 1kb window and calculate mean methylation in it.\n")
	complement_1kb_window = makeWindows(setdiff(chromGr, union(CGI, CGI_SHORE)), w = 1000)
	complement_1kb_window = complement_1kb_window[width(complement_1kb_window) > 500]
	gr_complement = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(complement_1kb_window = complement_1kb_window), filter_fun = function(s) TRUE)[[1]]

	saveRDS(gr_cgi, file = qq("@{PROJECT_DIR}/rds/mean_meth_1kb_cgi.rds"))
	saveRDS(gr_shore, file = qq("@{PROJECT_DIR}/rds/mean_meth_1kb_cgi_shore.rds"))
	saveRDS(gr_complement, file = qq("@{PROJECT_DIR}/rds/mean_meth_1kb_neither_cgi_nor_shore.rds"))
}


cat("subtype classification based on cgi\n")
gr_cgi = readRDS(qq("@{PROJECT_DIR}/rds/mean_meth_1kb_cgi.rds"))
mat = as.data.frame(mcols(gr_cgi))
mat = as.matrix(mat[, -ncol(mat)])

od = order(rowSds(mat), decreasing = TRUE)[1:4000]
set.seed(123)
res_cgi = ConsensusClusterPlus(mat[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_cgi = res_cgi[[2]]$consensusClass[sample_id])
COLOR$km_cgi = structure(res_cgi[[2]]$clrs[[3]], names = unique(res_cgi[[2]]$consensusClass))

pdf(qq("@{PROJECT_DIR}/image/methylation_classification_wgbs_cgi.pdf"), width = 8, height = 10)

ht = Heatmap(mat[od, ], name = "methylation", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_cgi[[2]]$consensusTree,
	col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_cgi = SAMPLE[sample_id, ]$km_cgi, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_cgi = COLOR$km_cgi), show_annotation_name = TRUE),
	column_title = qq("@{length(od)} rows, based on SD, CGI"))
draw(ht)

gr = as.data.frame(gr_cgi[od, 1:3])
rainfall = rainfallTransform(gr)
den = genomicDensity(gr, window.size = 5e6)

gtrellis_layout(n_track = 2, nrow = 6, compact = TRUE, category = CHROMOSOME,
    track_axis = TRUE,  track_height = c(3, 1),
    track_ylim = c(0, ceiling(max(log10(rainfall[[4]]))), 0, max(den[[4]])),
    track_ylab = c("log10(inter_dist)", "density"),
    add_name_track = TRUE, add_ideogram_track = TRUE, title = "rainfall plot for top 4000 1kb windows in CGI with highest variance")
add_points_track(rainfall, log10(rainfall[[4]]), pch = 16, size = unit(1, "mm"), gp = gpar(col = "#FF000080"))
add_lines_track(den, den[[4]], area = TRUE, gp = gpar(fill = "pink"))
dev.off()


cat("subtype classification based on cgi shores\n")
gr_shore = readRDS(qq("@{PROJECT_DIR}/rds/mean_meth_1kb_cgi_shore.rds"))
mat = as.data.frame(mcols(gr_shore))
mat = as.matrix(mat[, -ncol(mat)])

od = order(rowSds(mat), decreasing = TRUE)[1:8000]
set.seed(123)
res_shore = ConsensusClusterPlus(mat[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_shore = res_shore[[2]]$consensusClass[sample_id])
COLOR$km_shore = structure(res_shore[[2]]$clrs[[3]], names = unique(res_shore[[2]]$consensusClass))

pdf(qq("@{PROJECT_DIR}/image/methylation_classification_wgbs_shore.pdf"), width = 8, height = 10)

ht = Heatmap(mat[od, ], name = "methylation", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_shore[[2]]$consensusTree,
	col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_shore = SAMPLE[sample_id, ]$km_shore, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_shore = COLOR$km_shore), show_annotation_name = TRUE),
	column_title = qq("@{length(od)} rows, based on SD, CGI shore"))
draw(ht)

gr = as.data.frame(gr_shore[od, 1:3])
rainfall = rainfallTransform(gr)
den = genomicDensity(gr, window.size = 5e6)
gtrellis_layout(n_track = 2, nrow = 6, compact = TRUE, category = CHROMOSOME,
    track_axis = TRUE,  track_height = c(3, 1),
    track_ylim = c(0, ceiling(max(log10(rainfall[[4]]))), 0, max(den[[4]])),
    track_ylab = c("log10(inter_dist)", "density"),
    add_name_track = TRUE, add_ideogram_track = TRUE, title = "rainfall plot for top 8000 1kb windows in CGI shore with highest variance")
add_points_track(rainfall, log10(rainfall[[4]]), pch = 16, size = unit(1, "mm"), gp = gpar(col = "#FF000080"))
add_lines_track(den, den[[4]], area = TRUE, gp = gpar(fill = "pink"))
dev.off()


cat("subtype classification based on other parts in genome\n")
gr_complement = readRDS(qq("@{PROJECT_DIR}/rds/mean_meth_1kb_neither_cgi_nor_shore.rds"))
mat = as.data.frame(mcols(gr_complement))
mat = as.matrix(mat[, -ncol(mat)])

od = order(rowSds(mat), decreasing = TRUE)[1:8000]
set.seed(123)
res_complement = ConsensusClusterPlus(mat[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_complement = res_complement[[2]]$consensusClass[sample_id])
COLOR$km_complement = structure(res_complement[[2]]$clrs[[3]], names = unique(res_complement[[2]]$consensusClass))

pdf(qq("@{PROJECT_DIR}/image/methylation_classification_wgbs_complement.pdf"), width = 8, height = 10)

ht = Heatmap(mat[od, ], name = "methylation", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_complement[[2]]$consensusTree,
	col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_complement = SAMPLE[sample_id, ]$km_complement, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_complement = COLOR$km_complement), show_annotation_name = TRUE),
	column_title = qq("@{length(od)} rows, based on SD, complement of CGI and shore"))
draw(ht)

gr = as.data.frame(gr_complement[od, 1:3])
rainfall = rainfallTransform(gr)
den = genomicDensity(gr, window.size = 5e6)
gtrellis_layout(n_track = 2, nrow = 6, compact = TRUE, category = CHROMOSOME,
    track_axis = TRUE, track_height = c(3, 1), 
    track_ylim = c(0, ceiling(max(log10(rainfall[[4]]))), 0, max(den[[4]])),
    track_ylab = c("log10(inter_dist)", "density"),
    add_name_track = TRUE, add_ideogram_track = TRUE, title = "rainfall plot for top 8000 1kb windows in complement with highest variance")
add_points_track(rainfall, log10(rainfall[[4]]), pch = 16, size = unit(1, "mm"), gp = gpar(col = "#FF000080"))
add_lines_track(den, den[[4]], area = TRUE, gp = gpar(fill = "pink"))
dev.off()

saveRDS(SAMPLE, file = qq("@{PROJECT_DIR}/rds/group_classification.rds"))
saveRDS(COLOR, file = qq("@{PROJECT_DIR}/rds/group_classification_color.rds"))

pdf(qq("@{PROJECT_DIR}/image/group_classification.pdf"), width = 4)
ht_list = Heatmap(SAMPLE$group, name = "group", show_row_names = FALSE, col = COLOR$group) +
          Heatmap(SAMPLE$sample_type, name = "sample_type", show_row_names = FALSE, col = COLOR$sample_type) +
          Heatmap(SAMPLE[, c("km_expr", "km_cgi", "km_shore", "km_complement")], name = "km", show_row_names = TRUE,
	cluster_columns = FALSE, cluster_rows = FALSE, col = COLOR$km_expr)
draw(ht_list)
dev.off()

# cmd = qq("Rscript /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/00.sample_subgroup.R")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=30:00:00,mem=10G -N sample_subgroup.R' '@{cmd}'")
# system(cmd)

