
suppressPackageStartupMessages(library(GetoptLong))

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

sample_id = rownames(SAMPLE)

library(ConsensusClusterPlus)


####################################################
# cluster by expressoin
pdf(qq("@{OUTPUT_DIR}/plots/expression_group_classification.pdf"), width = 8, height = 10)

x = epic:::cor_columns(t(EXPR), abs_cutoff = c(0.5, 0.6, 0.7, 0.8))
od = x[, "0.7"] > 500

set.seed(123)
res_expr = ConsensusClusterPlus(EXPR[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_expr = res_expr[[2]]$consensusClass[sample_id])
COLOR$km_expr = structure(res_expr[[2]]$clrs[[3]], names = unique(res_expr[[2]]$consensusClass))

ht = Heatmap(EXPR[od, ], name = "expr", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_expr[[2]]$consensusTree,
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_expr = SAMPLE[sample_id, ]$km_expr, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_expr = COLOR$km_expr)),
	column_title = qq("@{sum(od)} rows, based on correlation"))
draw(ht)

od = order(rowSds(EXPR), decreasing = TRUE)[1:sum(od)]
set.seed(123)
res_expr = ConsensusClusterPlus(EXPR[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_expr_var = res_expr[[2]]$consensusClass[sample_id])
COLOR$km_expr_var = structure(res_expr[[2]]$clrs[[3]], names = unique(res_expr[[2]]$consensusClass))

ht = Heatmap(EXPR[od, ], name = "expr", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_expr[[2]]$consensusTree,
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_expr_var = SAMPLE[sample_id, ]$km_expr_var, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_expr_var = COLOR$km_expr_var)),
	column_title = qq("@{length(od)} rows, based on variance"))
draw(ht)
dev.off()



# pdf(qq("@{OUTPUT_DIR}/plots/global_methylation_distribution.pdf"))
# global_methylation_distribution(sample_id, annotation = SAMPLE[, "km"], annotation_color = COLOR$km,
# 	ha = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type,
# 		col = list(group = COLOR$group, sample_type = COLOR$sample_type)))
# dev.off()


ha = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type,
		col = list(group = COLOR$group, sample_type = COLOR$sample_type))


sample_id = rownames(SAMPLE)
n_sample = length(sample_id)

if(!all(file.exists(c(qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi.rds"), qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi_shore.rds"), qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_neither_cgi_nor_shore.rds"))))) {
	
	chromInfo = getChromInfoFromUCSC(GENOME)
	chromInfo = chromInfo[chromInfo$chrom %in% CHROMOSOME, ]
	chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))

	cat("split cgi by 1kb window and calculate mean methylation in it.\n")
	cgi_1kb_window = makeWindows(CGI, w = 1000, short.keep = TRUE)
	cgi_1kb_window = cgi_1kb_window[width(cgi_1kb_window) > 500]
	gr_cgi = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(cgi_1kb_window = cgi_1kb_window))[[1]]

	shore = CGI_SHORE
	shore_1kb_window = makeWindows(shore, w = 1000, short.keep = TRUE)
	shore_1kb_window = shore_1kb_window[width(shore_1kb_window) > 500]
	gr_shore = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(shore_1kb_window = shore_1kb_window))[[1]]

	cat("split genome by 1kb window and calculate mean methylation in it.\n")
	complement_1kb_window = makeWindows(setdiff(chromGr, union(CGI, CGI_SHORE)), w = 1000)
	complement_1kb_window = complement_1kb_window[width(complement_1kb_window) > 500]
	gr_complement = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(complement_1kb_window = complement_1kb_window))[[1]]

	saveRDS(gr_cgi, file = qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi.rds"))
	saveRDS(gr_shore, file = qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi_shore.rds"))
	saveRDS(gr_complement, file = qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_neither_cgi_nor_shore.rds"))

}


cat("subtype classification based on cgi\n")
gr_cgi = readRDS(qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi.rds"))
mat = as.data.frame(mcols(gr_cgi))


mat = as.matrix(mat[, -ncol(mat)])
x = epic:::cor_columns(t(mat), abs_cutoff = c(0.5, 0.6, 0.7, 0.8))
od = x[, "0.8"] > 2000

pdf(qq("@{OUTPUT_DIR}/plots/methylation_classification_wgbs_cgi.pdf"), width = 14, height = 14)

set.seed(123)
res_cgi = ConsensusClusterPlus(mat[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_cgi = res_cgi[[2]]$consensusClass[sample_id])
COLOR$km_cgi = structure(res_cgi[[2]]$clrs[[3]], names = unique(res_cgi[[2]]$consensusClass))

ht = Heatmap(mat[od, ], name = "methylation", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_cgi[[2]]$consensusTree,
	col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_cgi = SAMPLE[sample_id, ]$km_cgi, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_cgi = COLOR$km_cgi)),
	column_title = qq("@{sum(od)} rows, based on correlation"))
draw(ht)

od = order(rowSds(mat), decreasing = TRUE)[1:sum(od)]
set.seed(123)
res_cgi = ConsensusClusterPlus(mat[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_cgi_var = res_cgi[[2]]$consensusClass[sample_id])
COLOR$km_cgi_var = structure(res_cgi[[2]]$clrs[[3]], names = unique(res_cgi[[2]]$consensusClass))

ht = Heatmap(mat[od, ], name = "methylation", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_cgi[[2]]$consensusTree,
	col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_cgi_var = SAMPLE[sample_id, ]$km_cgi_var, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_cgi_var = COLOR$km_cgi_var)))
draw(ht)
dev.off()


cat("subtype classification based on cgi shores\n")
gr_shore = readRDS(qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_cgi_shore.rds"))
mat = as.data.frame(mcols(gr_shore))
mat = as.matrix(mat[, -ncol(mat)])
x = epic:::cor_columns(t(mat), abs_cutoff = c(0.5, 0.6, 0.7, 0.8))
od = x[, "0.8"] > 4000

pdf(qq("@{OUTPUT_DIR}/plots/methylation_classification_wgbs_shore.pdf"), width = 14, height = 14)

set.seed(123)
res_shore = ConsensusClusterPlus(mat[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_shore = res_shore[[2]]$consensusClass[sample_id])
COLOR$km_shore = structure(res_shore[[2]]$clrs[[3]], names = unique(res_shore[[2]]$consensusClass))

ht = Heatmap(mat[od, ], name = "methylation", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_shore[[2]]$consensusTree,
	col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_shore = SAMPLE[sample_id, ]$km_shore, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_shore = COLOR$km_shore)),
	column_title = qq("@{sum(od)} rows, based on correlation"))
draw(ht)

od = order(rowSds(mat), decreasing = TRUE)[1:sum(od)]
set.seed(123)
res_shore = ConsensusClusterPlus(mat[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_cgi_var = res_shore[[2]]$consensusClass[sample_id])
COLOR$km_shore_var = structure(res_shore[[2]]$clrs[[3]], names = unique(res_shore[[2]]$consensusClass))

ht = Heatmap(mat[od, ], name = "methylation", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_shore[[2]]$consensusTree,
	col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_shore_var = SAMPLE[sample_id, ]$km_shore_var, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_shore_var = COLOR$km_shore_var)))
draw(ht)
dev.off()


cat("subtype classification based on other parts in genome\n")
gr_complement = readRDS(qq("@{OUTPUT_DIR}/rds/mean_meth_1kb_neither_cgi_nor_shore.rds"))
mat = as.data.frame(mcols(gr_complement))
mat = as.matrix(mat[, -ncol(mat)])
mat = mat[sample(c(FALSE, TRUE), nrow(mat), replace = TRUE, p = c(0.9, 0.1)), ]
x = epic:::cor_columns(t(mat), abs_cutoff = c(0.5, 0.6, 0.7, 0.8))
od = x[, "0.8"] > 4000

pdf(qq("@{OUTPUT_DIR}/plots/methylation_classification_wgbs_complement.pdf"), width = 14, height = 14)

set.seed(123)
res_complement = ConsensusClusterPlus(mat[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_complement = res_complement[[2]]$consensusClass[sample_id])
COLOR$km_complement = structure(res_complement[[2]]$clrs[[3]], names = unique(res_complement[[2]]$consensusClass))

ht = Heatmap(mat[od, ], name = "methylation", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_complement[[2]]$consensusTree,
	col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_complement = SAMPLE[sample_id, ]$km_complement, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_complement = COLOR$km_complement)),
	column_title = qq("@{sum(od)} rows, based on correlation"))
draw(ht)

mat = as.data.frame(mcols(gr_complement))
mat = as.matrix(mat[, -ncol(mat)])
od = order(rowSds(mat), decreasing = TRUE)[1:sum(od)]
set.seed(123)
res_complement = ConsensusClusterPlus(mat[od, ], maxK = 6, reps = 1000, clusterAlg = "km", distance = "euclidean")

SAMPLE = cbind(SAMPLE, km_complement_var = res_complement[[2]]$consensusClass[sample_id])
COLOR$km_complement_var = structure(res_complement[[2]]$clrs[[3]], names = unique(res_complement[[2]]$consensusClass))

ht = Heatmap(mat[od, ], name = "methylation", show_row_names = FALSE, show_column_names = TRUE, cluster_columns = res_complement[[2]]$consensusTree,
	col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, km_complement_var = SAMPLE[sample_id, ]$km_complement_var, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, km_complement_var = COLOR$km_complement_var)))
draw(ht)
dev.off()


saveRDS(SAMPLE, file = qq("@{OUTPUT_DIR}/rds/group_classification.rds"))



# cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/sample_groups.R")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=30:00:00,mem=10G -N sample_groups' '@{cmd}'")
# system(cmd)

