
suppressPackageStartupMessages(library(GetoptLong))

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))


neg_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3.rds")))
pos_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3.rds")))

cr = c(neg_cr, pos_cr)

all_cr_gi = unique(cr$gene_id)
sample_id = attr(neg_cr, "sample_id")
extend = attr(neg_cr, "extend")

gene = genes(TXDB)
gene = gene[all_cr_gi]

target_ratio = 1/3
window = 200

rdata_file = qq("@{OUTPUT_DIR}/rds/mat_all_cr_enriched_to_gene_extend_50000_target_0.33.rds")
if(!file.exists(rdata_file)) {
	mat_list = list()
	for(chr in CHROMOSOME) {
		qqcat("normalize cr to gene for @{chr}.\n")
		mat_list[[chr]] = normalizeToMatrix(cr[seqnames(cr) == chr], gene[seqnames(gene) == chr], 
			mapping_column = "gene_id", value_column = "corr",
			target_ratio = target_ratio, extend = extend, w = window, empty_value = 0)
	}
	mat = do.call("rbind", mat_list)
	mat = copyAttr(mat_list[[1]], mat)
} else {
	mat = readRDS(rdata_file)
}

l = apply(mat, 1, function(x) sum(x == 0)/length(x)) < 0.25
qqcat("@{sum(!l)} genes are filtered out because more than 25% of the bins have no methylation.\n")
mat = mat[l, ]

gene = gene[rownames(mat)]

gene_length = width(gene[rownames(mat)])

sample_id_subgroup1 = intersect(sample_id, rownames(SAMPLE[SAMPLE$subgroup == "subgroup1", ]))
sample_id_subgroup2 = intersect(sample_id, rownames(SAMPLE[SAMPLE$subgroup == "subgroup2", ]))


expr = EXPR[rownames(mat), ]
expr_mean = rowMeans(expr[, sample_id_subgroup1]) - 
			rowMeans(expr[, sample_id_subgroup2])
expr_split = ifelse(expr_mean > 0, "high", "low")
expr_split = factor(expr_split, levels = c("high", "low"))


rdata_file = qq("@{OUTPUT_DIR}/rds/meth_mat_all_cr_enriched_to_gene_extend_50000_target_0.33.RData")
if(!file.exists(rdata_file)) {
	meth_mat = enrich_with_methylation(gene, sample_id, target_ratio = target_ratio, extend = extend, w = window)
	failed_rows = attr(meth_mat, "failed_rows")
	qqcat("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets.\n")
	meth_mat[failed_rows, ] = 0.5

	meth_mat_1 = enrich_with_methylation(gene, sample_id_subgroup1, target_ratio = target_ratio, extend = extend, w = window)
	failed_rows = attr(meth_mat_1, "failed_rows")
	qqcat("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets.\n")
	meth_mat_1[failed_rows, ] = 0.5

	meth_mat_2 = enrich_with_methylation(gene, sample_id_subgroup2, target_ratio = target_ratio, extend = extend, w = window)
	failed_rows = attr(meth_mat_2, "failed_rows")
	qqcat("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets.\n")
	meth_mat_2[failed_rows, ] = 0.5

	meth_mat_diff = meth_mat_1 - meth_mat_2

	save(meth_mat, meth_mat_1, meth_mat_2, meth_mat_diff, file = rdata_file)
} else {
	load(rdata_file)
}

target_index = attr(mat, "target_index")

rdata_file = qq("@{OUTPUT_DIR}/rds/mat_all_cr_enriched_to_gene_extend_50000_target_0.33_row_order_km3_and_km4.rds")
if(!file.exists(rdata_file)) {
	set.seed(123)
	km3 = kmeans(mat[, target_index], centers = 3)$cluster
	od = order(tapply(rowMeans(mat[, target_index]), km3, mean))
	km3 = od[km3]

	km4 = kmeans(mat[, target_index], centers = 4)$cluster
	od = order(tapply(rowMeans(mat[, target_index]), km4, mean))
	km4 = od[km4]

	saveRDS(list(km3, km4), )
} else {
	res = readRDS(rdata_file)
	km3 = res[[1]]
	km4 = res[[2]]
}

combined_split = paste(km4, expr_split, sep = ",")

row_order = NULL
for(i1 in 1:4) {
	for(i2 in c("high", "low")) {
		label = paste(i1, i2, sep = ",")
		dend1 = as.dendrogram(hclust(dist(meth_mat[combined_split == label, ])))
		dend1 = reorder(dend1, rowMeans(meth_mat[combined_split == label, ]))
		row_od1 = order.dendrogram(dend1)

		row_order = c(row_order, which(combined_split == label)[row_od1])
	}
}

cluster_color = brewer.pal(4, "Set1")[c(3, 4, 5, 1)]
names(cluster_color) = 1:4

corr_col_fun = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))

gene_length2 = gene_length
q9 = quantile(gene_length, 0.9)
gene_length2[gene_length2 > q9] = q9

generate_diff_color_fun = function(x) {
	q = quantile(x, c(0.05, 0.95))
	max_q = max(abs(q))
	colorRamp2(c(-max_q, 0, max_q), c("#3794bf", "#FFFFFF", "#df8640"))
}

pdf(qq("@{OUTPUT_DIR}.plots/enrichedheatmap_all.pdf"), width = 16, height = 10)
lines_lgd = Legend(at = c("1", "2", "3", "4"), title = "Lines", legend_gp = gpar(col = cluster_color), type = "lines")

axis_name = c("-50kb", "TSS", "TES", "50kb")

ht_list = Heatmap(km4, name = "cluster", show_row_names = FALSE, col = cluster_color, width = unit(5, "mm")) +
Heatmap(expr_split, name = "expr_direction", show_row_names = FALSE, col = c("high" = "red", "low" = "darkgreen"),
	width = unit(5, "mm")) +
EnrichedHeatmap(mat, name = "correlation", col = corr_col_fun, 
	top_annotation = HeatmapAnnotation(lines1 = anno_enriched(yaxis_side = "left", gp = gpar(col = cluster_color))),
	use_raster = TRUE, raster_quality = 2, column_title = "Correlation", axis_name = axis_name) +
rowAnnotation(gene_length = anno_points(gene_length2, pch = ".", gp = gpar(col = "#00000080")), width = unit(2, "cm")) +
EnrichedHeatmap(expr, name = "expression", 
	top_annotation = HeatmapAnnotation(df = SAMPLE, col = COLOR, show_annotation_name = TRUE, annotation_name_side = left,
		annotation_name_gp = gpar(fontsize = 10)), use_raster = TRUE, raster_quality = 2, column_title = "Expression",
	axis_name = axis_name) +
EnrichedHeatmap(meth_mat, name = "methylation", col = colorRamp2(c(0, 0.5, 1),c("blue", "white", "red")),
	top_annotation = HeatmapAnnotation(lines2 = anno_enriched(gp = gpar(col = cluster_color))), 
	column_title = "Methylation", axis_name = axis_name, combined_name_fun = NULL) +
EnrichedHeatmap(meth_mat_diff, name = "meth_diff", col = generate_diff_color_fun(meth_mat_diff),
	top_annotation = HeatmapAnnotation(lines3 = anno_enriched(gp = gpar(col = cluster_color))),
	column_title = "Methylation difference", axis_name = axis_name)

draw(ht_list, main_heatmap = "methylation", row_order = row_order, split = km4, annotation_legend_list = list(lgd))
decorate_annotation("gene_length", slice = 4, {
	grid.text("gene length", 0.5, -unit(5, "mm"), gp = gpar(fontsize = 8))
})
for(i in 1:4) {
	decorate_heatmap_body("cluster", slice = i, {
		grid.text(i)
	})
	decorate_heatmap_body("expression", slice = i, {
		grid.rect(gp = gpar(fill = "transparent"))
	})
}
dev.off()



pdf(qq("@{OUTPUT_DIR}.plots/enrichedheatmap_all.pdf"), width = 16, height = 10)
lt = tapply(seq_len(nrow(mat)), km4, function(ind) mat[ind, target_index])
group_mean = sapply(lt, mean)
boxplot(lt, col = corr_col_fun(group_mean), outline = FALSE)

lt1 = split(rowMeans(expr[, sample_id_subgroup1]), km4)
lt2 = split(rowMeans(expr[, sample_id_subgroup2]), km4)
lt = c(lt1, lt2)
lt = lt[c(1, 5, 2, 6, 3, 7, 4, 8)]
boxplot(lt, col = corr_col_fun(rep(group_mean, each = 2)), outline = FALSE, lty = 1:2)

lt_gene_length = split(gene_length, km4)
boxplot(lt_gene_length, col = corr_col_fun(group_mean), outline = FALSE)

lt = split(rowMeans(expr[, sample_id_subgroup1]), combined_split)
boxplot(lt, col = corr_col_fun(rep(group_mean, each = 2)), outline = FALSE, lty = 1:2)

lt = split(rowMeans(expr[, sample_id_subgroup2]), combined_split)
boxplot(lt, col = corr_col_fun(rep(group_mean, each = 2)), outline = FALSE, lty = 1:2)

lt_diff = split(rowMeans(meth_mat_diff[, target_index]), km4)
lt = sapply(1:4, function(i) sum(lt_diff[[i]]*lt_gene_length[[i]])/sum(lt_gene_length[[i]]))
barplot(lt, col = corr_col_fun(group_mean), ylim = c(0, max(lt)*1.1))

dev.off()


################################################
## gtrellis plot for group1 and group4