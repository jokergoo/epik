
library(methods)
library(GetoptLong)

cutoff = 0.05
meandiff = 0.1
rerun = FALSE
GetoptLong("cutoff=f", "0.05",
	       "meandiff=s", "0",
	       "rerun!", "rerun")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

neg_cr = readRDS(qq("@{OUTPUT_DIR}/rds/all_neg_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds"))
pos_cr = readRDS(qq("@{OUTPUT_DIR}/rds/all_pos_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds"))

foo_cr = c(neg_cr, pos_cr)
foo_cr$direction = c(rep("neg", length(neg_cr)), rep("pos", length(pos_cr)))

sample_id = attr(neg_cr, "sample_id")

add_mean_methylation = function(gr) {
	gr2 = GRanges()
	mean_meth = NULL
	for(chr in unique(as.vector(seqnames(gr)))) {
		
		methylation_hooks$set(chr)

		sub_gr = gr[seqnames(gr) == chr]
		gr_cpg = methylation_hooks$GRanges()
		m = methylation_hooks$meth(col_index = sample_id)

		mtch = as.matrix(findOverlaps(sub_gr, gr_cpg))
		mean_m = do.call("rbind", tapply(mtch[, 2], mtch[, 1], function(ind) colMeans(m[ind, , drop = FALSE], na.rm = TRUE)))
		
		mean_meth = rbind(mean_meth, mean_m)
		gr2 = c(gr2, sub_gr)
	}

	rownames(mean_meth) = NULL
	colnames(mean_meth) = paste0("mean_meth_", colnames(mean_meth))
	mcols(gr2) = cbind(mcols(gr2), as.data.frame(mean_meth))
	return(gr2)
}

## add mean methylation matrix to the `GRanges` object
rdata_file = qq("@{OUTPUT_DIR}/rds/sig_cr_mean_methylation_fdr_@{cutoff}_methdiff_@{meandiff}.rds")
if(file.exists(rdata_file) && !rerun) {
	foo_cr2 = readRDS(rdata_file)
} else {
	foo_cr2 = add_mean_methylation(foo_cr)
	saveRDS(foo_cr2, file = rdata_file)
}

meth_mat = as.matrix(mcols(foo_cr2)[, grep("mean_meth", colnames(mcols(foo_cr2)))])
expr_mat = EXPR[foo_cr2$gene_id, sample_id]
meth_diff = rowMeans(meth_mat[, SAMPLE$subgroup == "subgroup1"]) - 
			rowMeans(meth_mat[, SAMPLE$subgroup == "subgroup2"])


gm = genes(TXDB)
gl = width(gm)
names(gl) = names(gm)

## since there are multiple samples in a subgroup, this function
## returns the common regions that cover regions in most of the samples
get_chromatin_states = function(name, sample_id) {

	fn = dir(qq("@{BASE_DIR}/data/chromatin_states/"))
	nm = gsub("^(E\\d+?)_.*$", "\\1", fn)
	fn = fn[nm %in% sample_id]

	all_sample_chromatin_states = lapply(fn, function(x) {
		x = qq("@{BASE_DIR}/data/chromatin_states/@{x}")
		qqcat("reading @{x}...\n")
		gr = read.table(x, sep = "\t")
		gr = gr[gr[[1]] %in% CHROMOSOME, ]
		GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]] + 1, gr[[3]]), states = gr[[4]])
	})

	names(all_sample_chromatin_states) = gsub("^(E\\d+)_.*$", "\\1", fn)

	gf_list = lapply(all_sample_chromatin_states, function(gf) {
		gf[gf$states == name]
	})

	# a given region should cover at least 50% of the samples
	epic::common_regions(gf_list, min_width = 1000, min_coverage = ceiling(0.5*length(gf_list)), gap = 0)
}

subgroup = SAMPLE[sample_id, "subgroup"]

rdata_file = qq("@{OUTPUT_DIR}/rds/genomic_features_list_fdr_@{cutoff}_methdiff_@{meandiff}.rds")
if(file.exists(rdata_file) && !rerun) {
	gr_list = readRDS(rdata_file)
} else {
	df = read.table("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/data/chromatin_states/E099_15_coreMarks_mnemonics.bed.gz", sep = "\t", stringsAsFactors = FALSE)
	all_states = sort(unique(df[[4]]))

	# separate by subgroups
	cs_list_1 = lapply(all_states, get_chromatin_states, sample_id[subgroup == "subgroup1"])
	names(cs_list_1) = paste0(all_states, "_1")
	cs_list_2 = lapply(all_states, get_chromatin_states, sample_id[subgroup == "subgroup2"])
	names(cs_list_2) = paste0(all_states, "_2")


	tfbs = read.table("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/encode_uniform_tfbs_merged_1kb.bed", sep = "\t", stringsAsFactors = FALSE)
	tfbs = GRanges(seqnames = tfbs[[1]], ranges = IRanges(tfbs[[2]], tfbs[[3]]))


	gr_list = c(list(CGI = CGI, shore = CGI_SHORE, tfbs = tfbs), cs_list_1, cs_list_2)
	saveRDS(gr_list, file = rdata_file)
}
## for each region in `foo_cr2`, how much is covered by regions in `gr_list`
foo_cr2 = epic::annotate_to_genomic_features(foo_cr2, gr_list)

## whether it is at tss, gene body or intergenic regions
ga = ifelse(foo_cr2$gene_tss_dist > -1000 & foo_cr2$gene_tss_dist < 2000, "tss",
	 	ifelse(foo_cr2$gene_tss_dist > 2000 & foo_cr2$gene_tss_dist < gl[foo_cr2$gene_id], "gene", "intergenic"))

# a matrix for the overlapping of CGI/shore/tfbs
overlap_mat_0 = as.matrix(mcols(foo_cr2)[, grep("CGI|shore|tfbs", colnames(mcols(foo_cr2)))])
colnames(overlap_mat_0) = gsub("overlap_to_", "", colnames(overlap_mat_0))

# difference for the chromHMM segmentation overlapping in the two subgroups
overlap_mat_1 = as.matrix(mcols(foo_cr2)[, grep("overlap_to.*_1$", colnames(mcols(foo_cr2)))])
colnames(overlap_mat_1) = gsub("overlap_to_", "", colnames(overlap_mat_1))

overlap_mat_2 = as.matrix(mcols(foo_cr2)[, grep("overlap_to.*_2$", colnames(mcols(foo_cr2)))])
colnames(overlap_mat_2) = gsub("overlap_to_", "", colnames(overlap_mat_2))

overlap_mat_diff = overlap_mat_1 - overlap_mat_2
dim(overlap_mat_diff) = dim(overlap_mat_1)
dimnames(overlap_mat_diff) = dimnames(overlap_mat_1)
colnames(overlap_mat_diff) = gsub("_\\d$", "", colnames(overlap_mat_diff))

# when clustering columns, we cluster samples in each subgroup separately
dend1 = as.dendrogram(hclust(dist(t(meth_mat[, SAMPLE$subgroup == "subgroup1"]))))
hc1 = as.hclust(reorder(dend1, colMeans(meth_mat[, SAMPLE$subgroup == "subgroup1"])))
expr_col_od1 = hc1$order
dend2 = as.dendrogram(hclust(dist(t(meth_mat[, SAMPLE$subgroup == "subgroup2"]))))
hc2 = as.hclust(reorder(dend2, colMeans(meth_mat[, SAMPLE$subgroup == "subgroup2"])))
expr_col_od2 = hc2$order
expr_col_od = c(which(SAMPLE$subgroup == "subgroup1")[expr_col_od1], which(SAMPLE$subgroup == "subgroup2")[expr_col_od2])


abs_tss_dist = abs(foo_cr2$gene_tss_dist)
q = quantile(abs_tss_dist, 0.9); q = 5e4
abs_tss_dist[abs_tss_dist > q] = q


# rows are split into four slices for neg_cr and pos_cr separately and ordered by mean value
set.seed(123)
km_meth1 = kmeans(meth_mat[foo_cr2$direction == "neg", SAMPLE$subgroup == "subgroup1"], centers = 4)$cluster
x = tapply(rowMeans(meth_mat[foo_cr2$direction == "neg", ]), km_meth1, mean)
od = structure(rank(x), names = names(x))
km_meth1 = od[as.character(km_meth1)]
km_meth2 = kmeans(meth_mat[foo_cr2$direction == "pos", SAMPLE$subgroup == "subgroup1"], centers = 4)$cluster
x = tapply(rowMeans(meth_mat[foo_cr2$direction == "pos", ]), km_meth2, mean)
od = structure(rank(x), names = names(x))
km_meth2 = od[as.character(km_meth2)]
split = numeric(nrow(meth_mat))
split[foo_cr2$direction == "neg"] = paste0("neg", km_meth1)
split[foo_cr2$direction == "pos"] = paste0("pos", km_meth2)


## now we concatenate heatmaps 

## 1. a one-column heatmap shows row slices
ht_list = Heatmap(split, name = "split", show_row_names = FALSE, show_column_names = FALSE, width = unit(5, "mm"),
	col = c(neg1 = "darkgreen", neg2 = "darkgreen", neg3 = "darkgreen", neg4 = "darkgreen", 
		    pos1 = "red", pos2 = "red", pos3 = "red", pos4 = "red"), show_heatmap_legend = FALSE) +
## 2. methylation for the CRs
Heatmap(meth_mat, name = "methylation", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
	show_row_names = FALSE, show_column_names = FALSE, cluster_columns = FALSE, column_order = expr_col_od, 
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, subgroup = subgroup, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, subgroup = COLOR$subgroup)),
	column_title = "methylation", show_row_dend = FALSE, combined_name_fun = NULL,
	use_raster = TRUE, raster_quality = 2) + 
Heatmap(meth_diff, name = "meth_diff", col = colorRamp2(c(-0.3, 0, 0.3), c("green", "white", "red")),
	show_row_names = FALSE, width = unit(5, "mm"), show_heatmap_legend = FALSE) + 
## 3. expression matrix
Heatmap(expr_mat, name = "expr", show_row_names = FALSE,
	show_column_names = FALSE, cluster_columns = FALSE, column_order = expr_col_od,
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, subgroup = subgroup, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, subgroup = COLOR$subgroup), show_legend = FALSE, show_annotation_name = TRUE),
	column_title = "Expression", show_row_dend = FALSE,
	use_raster = TRUE, raster_quality = 2) + 
## 4. overlapping matrix for CGI/shore/tfbs
Heatmap(overlap_mat_0, name = "overlap0", show_row_names = FALSE, col = colorRamp2(c(0, 1), c("white", "orange")),
	show_column_names = TRUE, cluster_columns = FALSE,
	column_title = "overlap to gf", show_row_dend = FALSE) +
## 5. overlapping matrix for the chromHMM segmentations
Heatmap(overlap_mat_diff, col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")), 
	name = "overlap_diff", show_row_names = FALSE,
	show_column_names = TRUE, cluster_columns = FALSE,
	column_title = "overlap diff", show_row_dend = FALSE,
	use_raster = TRUE, raster_quality = 2) + 
## 6. dist to tss
rowAnnotation(tss_dist = row_anno_points(abs_tss_dist, size = unit(1, "mm"), gp = gpar(col = "#00000020"), axis = TRUE), 
	width = unit(2, "cm")) +
## 7. annotation to genes
Heatmap(ga, name = "anno", col = c("tss" = "red", "gene" = "blue", "intergenic" = "green"), show_row_names = FALSE,
	width = unit(5, "mm"))


pdf(qq("@{OUTPUT_DIR}/plots/sig_cr_heatmap_fdr_@{cutoff}_methdiff_@{meandiff}.pdf"), width = 20, height = 16)
draw(ht_list, main_heatmap = "methylation", split = split,
	column_title = qq("@{length(foo_cr)} cr, width=@{sum(width(foo_cr))}"))
decorate_annotation("tss_dist", slice = length(unique(split)), {
	grid.text("tss_dist", 0.5, unit(0, "npc") - unit(1, "cm"), gp = gpar(fontsize = 10))
})
all_levels = sort(unique(split))
for(i in seq_along(all_levels)) {
	decorate_heatmap_body("split", slice = i, {
		grid.text(all_levels[i], rot = 90, gp = gpar(col = "white"))
	})
	decorate_heatmap_body("methylation", slice = i, {
		grid.rect(gp = gpar(col = "black", fill = "transparent"))
	})
	decorate_heatmap_body("expr", slice = i, {
		grid.rect(gp = gpar(col = "black", fill = "transparent"))
	})
	decorate_heatmap_body("overlap0", slice = i, {
		grid.rect(gp = gpar(col = "black", fill = "transparent"))
	})
	decorate_heatmap_body("overlap_diff", slice = i, {
		grid.rect(gp = gpar(col = "black", fill = "transparent"))
	})
}
dev.off()

sample_id_subgroup1 = intersect(sample_id, rownames(SAMPLE[SAMPLE$subgroup == "subgroup1", ]))
sample_id_subgroup2 = intersect(sample_id, rownames(SAMPLE[SAMPLE$subgroup == "subgroup2", ]))

## barplots or boxplots for the annotation matrix in the eight row slices
pdf(qq("@{OUTPUT_DIR}/plots/sig_cr_heatmap_annotation_barplots_fdr_@{cutoff}_methdiff_@{meandiff}.pdf"), width = 20, height = 12)

par(mfrow = c(3, 5), mar = c(4, 4, 4, 1))
x1 = tapply(rowMeans(meth_mat[, paste0("mean_meth_", sample_id_subgroup1)]), split, mean)
x2 = tapply(rowMeans(meth_mat[, paste0("mean_meth_", sample_id_subgroup2)]), split, mean)
plot(1:8, x1, ylim = c(0, 1), type = "l", main = "mean methylation", axes = FALSE, ylab = "mean methylation")
points(1:8, x1, cex = 1.5, bg = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))(x1), pch = 21)
lines(1:8, x2, lty = 2)
points(1:8, x2, cex = 1.5, bg = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))(x2), pch = 21)
axis(side = 1, at = 1:8, labels = names(x1))
axis(side = 2)

for(i in 1:ncol(overlap_mat_0)) {
	x = tapply(overlap_mat_0[, i], split, mean)
	if(max(abs(x)) > 0.05) {
		barplot(x, main = colnames(overlap_mat_0)[i], col = "orange", ylim = c(0, 1))
	}
}
for(i in 1:ncol(overlap_mat_diff)) {
	x = tapply(overlap_mat_diff[, i], split, function(x) {
			c(sum(x[x < 0])/length(x), sum(x[x > 0]/length(x)))
		})
	x = do.call("cbind", x)
	if(max(abs(x)) > 0.05) {
		barplot(abs(x), main = colnames(overlap_mat_diff)[i], col = c("green", "red"), 
			ylim = c(-0.3, 0.3), offset = x[1, ])
	}
}
boxplot(split(abs(foo_cr2$gene_tss_dist), split), outline = FALSE, main = "dist2tss")
m = do.call("cbind", tapply(ga, split, table))
m = apply(m, 2, function(x) x/sum(x))
barplot(m, main = "annotation to genes", col = c("tss" = "red", "gene" = "blue", "intergenic" = "green")[rownames(m)])
dev.off()


# for(cutoff in c(0.1, 0.05, 0.01)) {
#     for(meandiff in c(0, 0.1, 0.2, 0.3)) {      
# 		cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/03.correlated_regions_sig_heatmap.R --cutoff @{cutoff} --meandiff @{meandiff} --no-rerun")
# 		cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=30:00:00,mem=10G -N correlated_regions_sig_heatmap_fdr_@{cutoff}_meandiff_@{meandiff}' '@{cmd}'")
# 		system(cmd)
# 	}
# }
