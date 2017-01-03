library(methods)
suppressPackageStartupMessages(library(GetoptLong))
type = "neg"
min_reduce = 1
min_width = 1000
name = "12_EnhBiv"
cutoff = 0.05
meandiff = 0.1
rerun = FALSE
nearest_by = "tss"
GetoptLong("type=s", "neg|pos",
	       "name=s", "name",
	       "min_reduce=i", "minimal for reducing",
	       "min_width=i", "minimal width of gf",
	       "cutoff=f", "0.05",
	       "meandiff=f", "0",
	       "rerun!", "rerun",
	       "nearest_by=s", "tss")
name = gsub("/", "_", name)
BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))


km = readRDS(qq("@{OUTPUT_DIR}/rds/mat_all_cr_enriched_to_gene_extend_50000_target_0.33_row_order_km3_and_km4.rds"))[[2]]
km_col = structure(brewer.pal(9, "Set1")[c(3,4,5,1)], names = c(1:4))

neg_cr_all = readRDS(qq("@{OUTPUT_DIR}/rds/all_neg_cr_w6s3.rds"))
pos_cr_all = readRDS(qq("@{OUTPUT_DIR}/rds/all_pos_cr_w6s3.rds"))
cr_all = c(neg_cr_all, pos_cr_all)
cr_all = copy_cr_attribute(neg_cr_all, cr_all)

if(type == "neg") {
	cr = readRDS(qq("@{OUTPUT_DIR}/rds/all_neg_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds"))
	cr_vice = readRDS(qq("@{OUTPUT_DIR}/rds/all_pos_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds"))
} else if(type == "pos") {
	cr = readRDS(qq("@{OUTPUT_DIR}/rds/all_pos_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds"))
	cr_vice = readRDS(qq("@{OUTPUT_DIR}/rds/all_neg_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds"))
}


sample_id = attr(cr, "sample_id")
sample_id_subgroup1 = intersect(sample_id, rownames(SAMPLE[SAMPLE$subgroup == "subgroup1", ]))
sample_id_subgroup2 = intersect(sample_id, rownames(SAMPLE[SAMPLE$subgroup == "subgroup2", ]))

if(name == "encode_tfbs") {
	gf = read.table("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/encode_uniform_tfbs_merged_1kb.bed", sep = "\t", stringsAsFactors = FALSE)
	gf = GRanges(seqnames = gf[[1]], ranges = IRanges(gf[[2]], gf[[3]]))
} else {

	fn = dir(qq("@{BASE_DIR}/data/chromatin_states/"))
	nm = gsub("^(E\\d+?)_.*$", "\\1", fn)
	fn = fn[nm %in% sample_id]

	all_sample_chromatin_states = lapply(fn, function(x) {
		x = qq("@{BASE_DIR}/data/chromatin_states/@{x}")
		qqcat("reading @{x}...\n")
		gr = read.table(x, sep = "\t", stringsAsFactors = FALSE)
		gr = gr[gr[[1]] %in% CHROMOSOME, ]
		gr[[4]] = gsub("/", "_", gr[[4]])
		GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]] + 1, gr[[3]]), states = gr[[4]])
	})

	names(all_sample_chromatin_states) = gsub("^(E\\d+)_.*$", "\\1", fn)

	gf_list = lapply(all_sample_chromatin_states, function(gf) {
		gf[gf$states == name]
	})

	# in each group, in at least 50% samples
	l1 = names(gf_list) %in% sample_id_subgroup1
	gf1 = epic::common_regions(gf_list[l1], min_width = 1000, min_coverage = ceiling(0.5*sum(l1)), gap = 0)
	l2 = names(gf_list) %in% sample_id_subgroup2
	gf2 = epic::common_regions(gf_list[l2], min_width = 1000, min_coverage = ceiling(0.5*sum(l2)), gap = 0)
	seqlengths(gf1) = read.chromInfo()$chr.len[seqlevels(gf1)]
	seqlengths(gf2) = read.chromInfo()$chr.len[seqlevels(gf2)]
	mcols(gf1) = NULL
	mcols(gf2) = NULL

	gf = c(gf1, gf2)
}

# overlap gf to gene extended regions
gm = genes(TXDB)
gm = gm[intersect(names(gm), names(km))]

gm = gm[seqnames(gm) %in% CHROMOSOME]
gm = gm[gm$gene_id %in% unique(cr$gene_id)]
tss = promoters(gm, upstream = 1, downstream = 0)
gl = width(gm)
g = gm
strand(g) = "*"
start(g) = start(g) - 50000
end(g) = end(g) + 50000
start(g) = ifelse(start(g) > 1, start(g), 1)

mtch = as.matrix(findOverlaps(g, gf))
gf = gf[unique(mtch[, 2])]

if(min_reduce >= 0) {
	gf = reduce(gf, min = min_reduce)
}

gf = gf[width(gf) >= min_width, ]

# find associated gene (by tss for by gene)
if(nearest_by == "tss") {
	d = distanceToNearest(gf, tss, select = "all")
	subjectHits = d@subjectHits
	ind = tapply(seq_len(length(d)), d@queryHits, function(ind) {
		ind[which.max(gl[subjectHits[ind]])][1]
	})
	d = d[as.vector(ind)]
	gf2 = gf[d@queryHits]
	gf2$distanceToNearest = d@elementMetadata@listData$distance
	gf2$nearestGene = gm$gene_id[d@subjectHits]
	gf2$nearestGeneStrand = strand(tss[d@subjectHits])
	l = gf2$nearestGeneStrand == "+" & start(gf2) < start(tss[d@subjectHits]) |
	    gf2$nearestGeneStrand == "-" & end(gf2) > end(tss[d@subjectHits])
	l = as.vector(l)
	gf2$distanceToNearest[l] = -gf2$distanceToNearest[l]
} else {
	d = distanceToNearest(gf, gm, select = "all")
	subjectHits = d@subjectHits
	ind = tapply(seq_len(length(d)), d@queryHits, function(ind) {
		ind[which.max(gl[subjectHits[ind]])][1]
	})
	d = d[as.vector(ind)]
	gf2 = gf[d@queryHits]
	gf2$distanceToNearest = d@elementMetadata@listData$distance
	gf2$nearestGene = gm$gene_id[d@subjectHits]
	gf2$nearestGeneStrand = strand(gm[d@subjectHits])
	# if the gf is overlapped to gene body, how much of it is overlapped by the gene
	gg = pintersect(gf2, gm[gf2$nearestGene], resolve.empty = "max")
	gf2$overlapGenePercent = width(gg)/width(gf2)
	l = gf2$nearestGeneStrand == "+" & start(gf2) < start(gm[d@subjectHits]) |
	    gf2$nearestGeneStrand == "-" & end(gf2) > end(gm[d@subjectHits])
	l = as.vector(l)
	gf2$distanceToNearest[l] = -gf2$distanceToNearest[l]
}

strand(gf2) = strand(gm[gf2$nearestGene])
names(gf2) = gf2$nearestGene
gf2$gene_id = gf2$nearestGene

extend = 5000
target_ratio = mean(width(gf2))/(extend*2 + mean(width(gf2)))

mat = normalizeToMatrix(cr, gf2, extend = extend, w = 50, trim = 0, mean_mode = "absolute",
	mapping_column = "gene_id", target_ratio = target_ratio)
if(type == "neg") {
	mat[mat == 1] = -1
}
mat_vice = normalizeToMatrix(cr_vice, gf2, extend = extend, w = 50, trim = 0, mean_mode = "absolute",
	mapping_column = "gene_id", target_ratio = target_ratio)
if(type == "pos") {
	mat_vice[mat_vice == 1] = -1
}
l = rowSums(abs(mat)) > 0
mat = mat[l, ]
gf2 = gf2[l]
mat_vice = mat_vice[l, ]

qqcat("@{length(gf2)} @{name}s remains\n")
if(length(gf2) < 10) {
	q(save = "no")
}


mat_corr = normalizeToMatrix(cr_all, gf2, mapping_column = "gene_id", value_column = "corr",
	extend = extend, mean_mode = "absolute", w = 50, target_ratio = target_ratio, trim = 0, empty_value = 0)

# n_tss = countOverlaps(cgi_extend, tss)

# dist = distanceToNearest(gf2, tss)

strd = as.vector(strand(gf2))
strd = factor(paste0(strd, "strand"), levels = c("-strand", "+strand"))

km = km[names(gf2)]

gf_width = width(gf2)
gf_width[gf_width > quantile(gf_width, 0.99)] = quantile(gf_width, 0.99)

rdata_file = qq("@{OUTPUT_DIR}/rds/cr_enrichedheatmap_@{name}_nearest_by_@{nearest_by}_@{type}_fdr_@{cutoff}_methdiff_@{meandiff}.RData")
if(file.exists(rdata_file) && !rerun) {
	load(rdata_file)
} else {
	meth_mat = enrich_with_methylation(gf2, sample_id, target_ratio = target_ratio, extend = extend)
	meth_mat[attr(meth_mat, "failed_rows"), ] = 0.5

	meth_mat_1 = enrich_with_methylation(gf2, sample_id_subgroup1, target_ratio = target_ratio, extend = extend)
	failed_rows = attr(meth_mat_1, "failed_rows")
	qqcat("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets.\n")
	meth_mat_1[failed_rows, ] = 0.5

	meth_mat_2 = enrich_with_methylation(gf2, sample_id_subgroup2, target_ratio = target_ratio, extend = extend)
	failed_rows = attr(meth_mat_2, "failed_rows")
	qqcat("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets.\n")
	meth_mat_2[failed_rows, ] = 0.5

	meth_mat_diff = meth_mat_1 - meth_mat_2

	cor_mat_list = list()
	hist_mat_list = list()
	hist_mat_list_subgroup1 = list()
	hist_mat_list_subgroup2 = list()
	hist_mat_list_diff = list()

	
	for(k in seq_along(MARKS)) {
		hm_sample = intersect(sample_id, chipseq_hooks$sample_id(MARKS[k]))
		# applied to each sample, each mark
		lt = enrich_with_histone_mark(gf2, sample_id = sample_id, mark = MARKS[k], return_arr = TRUE, target_ratio = target_ratio, extend = extend)
		arr = lt[[1]]

		# only calculate the correlation when there are enough samples
		if(length(hm_sample) >= 5) {
			# detect regions that histone MARKS correlate to expression
			expr2 = EXPR[gf2$gene_id, intersect(colnames(EXPR), hm_sample)]
			cor_mat = matrix(nrow = nrow(expr2), ncol = ncol(mat))
			cor_p_mat = cor_mat

			counter = set_counter(nrow(cor_mat))
			for(i in seq_len(nrow(cor_mat))) {
				counter()
			    for(j in seq_len(ncol(cor_mat))) {
			        x = cor(arr[i, j, ], expr2[i, ], method = "spearman")
			        cor_mat[i, j] = x
			        # cor_p_mat[i, j] = cor.test(arr[i, j, ], expr2[i, ], method = "spearman")$p.value
			    }
			}
			cat("\n")
			cor_mat[is.na(cor_mat)] = 0
			# cor_fdr_mat = p.adjust(cor_p_mat, method = "BH")
			# l1 = cor_fdr_mat < 0.1 & cor_mat > 0
			# cor_mat[l1] = 1
			# l2 = cor_fdr_mat < 0.1 & cor_mat < 0
			# cor_mat[l2] = -1 
			# cor_mat[!(l1 | l2)] = 0
			cor_mat = copyAttr(mat, cor_mat)
			cor_mat_list[[k]] = cor_mat

			if(sum(abs(cor_mat)) == 0) {
				cor_mat_list[[k]] = NA
			}
		} else {
			cor_mat_list[[k]] = NA
		}
		hist_mat_list[[k]] = lt[[2]]
		hist_mat_list_subgroup1[[k]] = apply(arr[, , intersect(hm_sample, sample_id_subgroup1)], c(1, 2), mean, na.rm = TRUE)
		hist_mat_list_subgroup1[[k]] = copyAttr(mat, hist_mat_list_subgroup1[[k]])
		hist_mat_list_subgroup2[[k]] = apply(arr[, , intersect(hm_sample, sample_id_subgroup2)], c(1, 2), mean, na.rm = TRUE)
		hist_mat_list_subgroup2[[k]] = copyAttr(mat, hist_mat_list_subgroup2[[k]])

		hist_mat_list_diff[[k]] = hist_mat_list_subgroup1[[k]] - hist_mat_list_subgroup2[[k]]
	}


	save(meth_mat, meth_mat_1, meth_mat_2, meth_mat_diff, cor_mat_list, 
		hist_mat_list, hist_mat_list_subgroup1, hist_mat_list_subgroup2, hist_mat_list_diff, file = rdata_file)
}


add_boxplot_of_gf_length = function(ht_list) {
	gl = gf_width
	anno_name = "gf_width"

	row_order_list = row_order(ht_list)
	lt = lapply(row_order_list, function(ind) gl[ind])
	bx = boxplot(lt, plot = FALSE)$stats
	n = length(row_order_list)
	x_ind = (seq_len(n) - 0.5)/n
	w = 1/n*0.5
	decorate_annotation(anno_name, slice = 1, {
		rg = range(bx)
		rg[1] = rg[1] - (rg[2] - rg[1])*0.1
		rg[2] = rg[2] + (rg[2] - rg[1])*0.1
		pushViewport(viewport(y = unit(1, "npc") + unit(1, "mm"), just = "bottom", height = unit(2, "cm"), yscale = rg))
		grid.rect(gp = gpar(col = "black"))
		grid.segments(x_ind - w/2, bx[5, ], x_ind + w/2, bx[5, ], default.units = "native", gp = gpar(lty = 1:2))
		grid.segments(x_ind - w/2, bx[1, ], x_ind + w/2, bx[1, ], default.units = "native", gp = gpar(lty = 1:2))
		grid.segments(x_ind, bx[1, ], x_ind, bx[5, ], default.units = "native", gp = gpar(lty = 1:2))
		grid.rect(x_ind, colMeans(bx[c(4, 2), ]), width = w, height = bx[4, ] - bx[2, ], default.units = "native", gp = gpar(fill = "white", lty = 1:2))
		grid.segments(x_ind - w/2, bx[3, ], x_ind + w/2, bx[3, ], default.units = "native", gp = gpar(lty = 1:2))
		grid.yaxis(main = FALSE, gp = gpar(fontsize = 8))
		grid.text(anno_name, y = unit(1, "npc") + unit(2.5, "mm"), gp = gpar(fontsize = 14), just = "bottom")
		upViewport()
	})
}


expr = EXPR[gf2$nearestGene, sample_id, drop = FALSE]
## calculate row orders
expr_mean = rowMeans(expr[, SAMPLE[sample_id, ]$subgroup == "subgroup1"]) - 
			rowMeans(expr[, SAMPLE[sample_id, ]$subgroup == "subgroup2"])
expr_split = ifelse(expr_mean > 0, "high", "low")
expr_split = factor(expr_split, levels = c("high", "low"))

set.seed(123)
target_index = attr(meth_mat, "target_index")
meth_split = kmeans(meth_mat[, target_index], centers = 2)$cluster
x = tapply(rowMeans(meth_mat[, target_index]), meth_split, mean)
od = structure(order(x), names = names(x))
meth_split = paste0("cluster", od[as.character(meth_split)])

combined_split = paste(meth_split, expr_split, sep = "|")

merge_row_order = function(l_list) {
	do.call("c", lapply(l_list, function(l) {
		if(sum(l) == 0) return(integer(0))
		if(sum(l) == 1) return(which(l))
		dend1 = as.dendrogram(hclust(dist_by_closeness2(mat[l, ])))
		dend1 = reorder(dend1, rowMeans(mat[l, ]))
		od = order.dendrogram(dend1)
		which(l)[od]
	}))
}

row_order = merge_row_order(list(
	combined_split == "cluster1|high" & km == 1,
	combined_split == "cluster1|high" & km == 2,
	combined_split == "cluster1|high" & km == 3,
	combined_split == "cluster1|high" & km == 4,
	combined_split == "cluster1|low" & km == 1,
	combined_split == "cluster1|low" & km == 2,
	combined_split == "cluster1|low" & km == 3,
	combined_split == "cluster1|low" & km == 4,
	combined_split == "cluster2|high" & km == 1,
	combined_split == "cluster2|high" & km == 2,
	combined_split == "cluster2|high" & km == 3,
	combined_split == "cluster2|high" & km == 4,
	combined_split == "cluster2|low" & km == 1,
	combined_split == "cluster2|low" & km == 2,
	combined_split == "cluster2|low" & km == 3,
	combined_split == "cluster2|low" & km == 4
))

n_heatmap = 0
col = c("-1" = "darkgreen", "0" = "white", "1" = "red")

## heatmap for expression
# columns are clustered for each subgroup
dend1 = as.dendrogram(hclust(dist(t(expr[, SAMPLE$subgroup == "subgroup1"]))))
hc1 = as.hclust(reorder(dend1, colMeans(expr[, SAMPLE$subgroup == "subgroup1"])))
expr_col_od1 = hc1$order
dend2 = as.dendrogram(hclust(dist(t(expr[, SAMPLE$subgroup == "subgroup2"]))))
hc2 = as.hclust(reorder(dend2, colMeans(expr[, SAMPLE$subgroup == "subgroup2"])))
expr_col_od2 = hc2$order
expr_col_od = c(which(SAMPLE$subgroup == "subgroup1")[expr_col_od1], which(SAMPLE$subgroup == "subgroup2")[expr_col_od2])

cor_col_fun = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))
axis_name = c("-5kb", "Start", "End", "5kb")

ht_list = Heatmap(expr, name = "expr", show_row_names = FALSE,
	show_column_names = FALSE, width = unit(5, "cm"), show_column_dend = FALSE, cluster_columns = FALSE, column_order = expr_col_od,
	top_annotation = HeatmapAnnotation(group = SAMPLE[sample_id, ]$group, sample_type = SAMPLE[sample_id, ]$sample_type, subgroup = SAMPLE[sample_id, ]$subgroup, 
		col = list(group = COLOR$group, sample_type = COLOR$sample_type, subgroup = COLOR$subgroup), show_annotation_name = TRUE, annotation_name_side = "left"),
	column_title = "Expression", show_row_dend = FALSE,
	use_raster = TRUE, raster_quality = 2)
gap = unit(0.3, "cm")
n_heatmap = n_heatmap + 1

ht_list = ht_list + Heatmap(km[names(gf2)], name = "km_groups", col = km_col, show_row_names = FALSE,
	width = unit(1, "cm"))

gap = unit.c(gap, unit(1, "cm"))

mat_mix = mat
mat_mix[mat == 0] = mat_vice[mat == 0]
ht_list = ht_list +
			EnrichedHeatmap(mat_mix, name = qq("@{type}CR"), col = col, split = strd,
                top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(neg_col = "darkgreen", pos_col = "red", lty = 1:2))), 
                top_annotation_height = unit(2, "cm"), column_title = qq("sig@{type}CR"),
                use_raster = TRUE, raster_quality = 2, combined_name_fun = NULL, axis_name = axis_name)+
            rowAnnotation(gf_width = row_anno_points(gf_width, axis = TRUE, gp = gpar(col = "#00000040")),
                width = unit(1, "cm"))
n_heatmap = n_heatmap + 1
gap = unit(c(0.3, 1, 1), "cm")



ht_list = ht_list + EnrichedHeatmap(mat_corr, col = cor_col_fun, name = qq("correlation"), 
      top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(pos_col = "red", neg_col = "darkgreen", lty = 1:2))), 
      top_annotation_height = unit(2, "cm"), column_title = qq("corr_meth"), axis_name = axis_name,
      use_raster = TRUE, raster_quality = 2)
gap = unit.c(gap, unit(1, "cm"))
n_heatmap = n_heatmap + 1

# methylation
ht_list = ht_list + EnrichedHeatmap(meth_mat, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
	name = "methylation", column_title = qq("meth"), axis_name = axis_name,
		heatmap_legend_param = list(title = "methylation"),
		top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "red", lty = 1:2))),
		use_raster = TRUE, raster_quality = 2)
gap = unit.c(gap, unit(1, "cm"))
n_heatmap = n_heatmap + 1

generate_diff_color_fun = function(x) {
	q = quantile(x, c(0, 0.01, 0.05, 0.95, 0.99, 1))
	if(abs(max(q[3:4])) == 0) {
		if(abs(max(q[c(2, 5)])) == 0) {
			if(abs(max(q[c(1, 6)])) == 0) {
				max_q = 1
			} else {
				max_q = abs(max(q[c(1, 6)]))
			}
		} else {
			max_q = abs(max(q[c(2, 5)]))
		}
	} else {
		max_q = abs(max(q[3:4]))
	}
	
	colorRamp2(c(-max_q, 0, max_q), c("#3794bf", "#FFFFFF", "#df8640"))
}

ht_list = ht_list + EnrichedHeatmap(meth_mat_diff, col = generate_diff_color_fun(meth_mat_diff),
	name = "methylation_diff", column_title = qq("meth_diff"), axis_name = axis_name,
		heatmap_legend_param = list(title = "methylation_diff"),
		top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(pos_col = "#df8640", neg_col = "#3794bf", lty = 1:2))),
		use_raster = TRUE, raster_quality = 2)
gap = unit.c(gap, unit(1, "cm"))
n_heatmap = n_heatmap + 1

ht_list2 = NULL
ht_list1 = NULL
# correlation to histone marks
for(i in seq_along(cor_mat_list)) {
	if(i == 3) {
		ht_list1 = ht_list
		ht_list = NULL
	}
	if(length(cor_mat_list[[i]]) > 1) {
		anno_line_col = ifelse(mean(cor_mat_list[[i]], na.rm = TRUE) > 0, "red", "darkgreen")
		ht_list = ht_list + EnrichedHeatmap(cor_mat_list[[i]], col = cor_col_fun, name = qq("corr_@{MARKS[i]}"), 
	          top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(pos_col = "red", neg_col = "darkgreen", lty = 1:2))), 
              top_annotation_height = unit(2, "cm"), column_title = qq("corr_@{MARKS[i]}"), axis_name = axis_name,
              use_raster = TRUE, raster_quality = 2)
	    gap = unit.c(gap, unit(1, "cm"))
   	 	n_heatmap = n_heatmap + 1
   	 }


   	hist_col_fun = function(x) {
   		q = quantile(hist_mat_list[[i]], c(0, 0.95, 0.99, 1))
   		if(q[2] == 0) {
   			if(q[3] == 0) {
   				if(q[4] == 0) {
   					colorRamp2(c(0, 1), c("white", "purple"))
   				} else {
   					colorRamp2(q[c(1, 4)], c("white", "purple"))
   				}
   			} else {
   				colorRamp2(q[c(1, 3)], c("white", "purple"))
   			}
   		} else {
   			colorRamp2(q[c(1, 2)], c("white", "purple"))
   		}
   	}

    ht_list = ht_list + EnrichedHeatmap(hist_mat_list[[i]], col = hist_col_fun(hist_mat_list[[i]]), name = qq("@{MARKS[i]}_1"),
		column_title = qq("@{MARKS[i]}"), axis_name = axis_name,
		heatmap_legend_param = list(title = qq("@{MARKS[i]}_density")),
		top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "purple", lty = 1:2))),
		use_raster = TRUE, raster_quality = 2)
	gap = unit.c(gap, unit(1, "cm"))
	n_heatmap = n_heatmap + 1

	ht_list = ht_list + EnrichedHeatmap(hist_mat_list_diff[[i]], name = qq("@{MARKS[i]}_diff"), col = generate_diff_color_fun(hist_mat_list_diff[[i]]),
		column_title = qq("@{MARKS[i]}_diff"), axis_name = axis_name,
		heatmap_legend_param = list(title = qq("@{MARKS[i]}_diff")),
		top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(pos_col = "#df8640", neg_col = "#3794bf", lty = 1:2))),
		use_raster = TRUE, raster_quality = 2)
	gap = unit.c(gap, unit(1, "cm"))
	n_heatmap = n_heatmap + 1
}
ht_list2 = ht_list

lines_lgd = Legend(at = c("high", "low"), title = "Lines", legend_gp = gpar(lty = 1:2), type = "lines")

# following chunk is necessary because Legend() needs to open a new graphic device
for(i in seq_along(dev.list())) {
	dev.off()
}

pdf(qq("@{OUTPUT_DIR}/plots/cr_enrichedheatmap_@{name}_nearest_by_@{nearest_by}_@{type}_fdr_@{cutoff}_methdiff_@{meandiff}.pdf"), width = n_heatmap + 4, height = 10)
ht_list = Heatmap(expr_split, show_row_names = FALSE, name = "expr_split", col = c("high" = "red", "low" = "darkgreen"), width = unit(5, "mm")) + ht_list1
foo = draw(ht_list,  annotation_legend_list = list(lines_lgd),
		cluster_rows = FALSE, row_order = row_order, show_row_dend = FALSE,
		column_title = qq("cluster by methylation, @{nrow(mat)} rows"), split = meth_split, row_sub_title_side = "left",
		show_heatmap_legend = FALSE, annotation_legend_side = "right")
add_boxplot_of_gf_length(foo)
i = 0
for(f in names(ht_list1@ht_list)) {
	if(grepl("expr|annotation|CGI|km|CR|gf", f)) next
	decorate_column_title(f, {
		grid.rect(height = unit(0.8, "npc"), gp = gpar(fill = brewer.pal(8, "Set2")[as.integer(i/3)+1], col = NA))
		grid.text(ht_list1@ht_list[[f]]@column_title, gp = gpar(fontsize = 14))
	})
	i = i + 1
}

ht_list = Heatmap(expr_split, show_row_names = FALSE, name = "expr_split", col = c("high" = "red", "low" = "darkgreen"), width = unit(5, "mm")) + ht_list2
foo = draw(ht_list, annotation_legend_list = list(lines_lgd),
		cluster_rows = FALSE, row_order = row_order, show_row_dend = FALSE,
		column_title = qq("cluster by methylation, @{nrow(mat)} rows"), split = meth_split, row_sub_title_side = "left",
		show_heatmap_legend = FALSE)
for(f in names(ht_list2@ht_list)) {
	decorate_column_title(f, {
		grid.rect(height = unit(0.8, "npc"), gp = gpar(fill = brewer.pal(8, "Set2")[as.integer(i/3)+1], col = NA))
		grid.text(ht_list2@ht_list[[f]]@column_title, gp = gpar(fontsize = 14))
	})
	i = i + 1
}

dev.off()

# df = read.table("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/data/chromatin_states/E099_15_coreMarks_mnemonics.bed.gz", sep = "\t", stringsAsFactors = FALSE)
# for(type in c("pos", "neg")) {
# 	for(name in sort(c(unique(df[[4]]), "encode_tfbs"))) {
#       name = gsub("/", "_", name)
# 	    cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/04.correlated_regions_enrichedheatmap_at_genomicfeatures.R --no-rerun --type @{type} --name @{name} --min_reduce 1 --nearest_by tss")
# 	    cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=15G -N cr_enrichedheatmap_at_genomicfeature_@{type}_@{name}_nearest_by_tss' '@{cmd}'")
# 	    system(cmd)
# 		cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/04.correlated_regions_enrichedheatmap_at_genomicfeatures.R --no-rerun --type @{type} --name @{name} --min_reduce 1 --nearest_by gene")
# 	    cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=15G -N cr_enrichedheatmap_at_genomicfeature_@{type}_@{name}_nearest_by_genes' '@{cmd}'")
# 	    system(cmd)
# 	}
# }
