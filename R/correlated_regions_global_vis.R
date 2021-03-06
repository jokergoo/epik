# == title
# Visualize global correlation by Hilbert curve
#
# == param
# -cr correlated regions
# -chromosome a vector of chromosome names
# -merge_chr whether merge all chromosomes in one plot
# -add_chr_name whether add chromosome names to the plot
# -title title of the plot
# -legend legend generated from `ComplexHeatmap::Legend`
# -... pass to `HilbertCurve::GenomicHilbertCurve`
#
# == value
# no value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
cr_hilbert_curve = function(cr, chromosome = paste0("chr", 1:22), 
	merge_chr = TRUE, add_chr_name = TRUE, title = "cr", legend = lgd, level = 10, ...) {

	cr_param = metadata(cr)$cr_param
	species = cr_param$genome

	chr_len = read.chromInfo(species = species)$chr.len

	## all cr windows
	col_fun = colorRamp2(c(-1, 0, 1), c("#4DAF4A", "white", "#E41A1C"))
	cm = ColorMapping(col_fun = col_fun)
	lgd = color_mapping_legend(cm, title = "type", plot = FALSE)
	if(merge_chr) {
		hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = level, title = title, legend = legend, ...)
	    hc_layer(hc, cr, col = col_fun(cr$corr), mean_mode = "absolute")
	    hc_map(hc, add = TRUE, fill = NA, border = "#808080", show_labels = add_chr_name)
	} else {
		chromosome = intersect(chromosome, unique(as.vector(seqnames(cr))))
		for(i in seq_along(chromosome)) {
		    chr = chromosome[i]
		    message(qq("make hilbert curve for cr for @{chr}"))
		    cr2 = cr[seqnames(cr) == chr]
		    hc = HilbertCurve(s = 1, e = max(chr_len), mode = "pixel", level = level, title = qq("@{title}@{ifelse(title == '', '', ', ')}@{chr}"), legend = legend)
		    hc_layer(hc, ranges(cr2), col = col_fun(cr2$corr), mean_mode = "absolute")
		}
	}
}

# == title
# Visualize landscape of genome-wide correlations
#
# == param
# -cr correlated regions
# -txdb transcriptome annotation which was used in `correlated_regions`
# -expr expression matrix which was used in `correlated_regions`
# -expr_ha a `ComplexHeatmap::HeatmapAnnotation` object for the expression heatmap
#
# == details
# The landscape of genome-wide correlations is visualized by a list of heatmaps.
# Each row corresponds to a single gene:
#
# - an enriched heatmap in which correlation signals are normalized at gene bodies
# - a point plot showing the length of genes
# - a heatmap of gene expression
# - a heatmap showing the mean methylation in the extended gene regions.
# - a heatmap showing the methylation difference in the extended gene regions.
#
# K-means clustering with four groups is applied on the correlation matrix which has been normalized.
# The four row subclusters are ordered by mean correlation. So basically, the four groups correspond to
# negative gene body correlation, weak negative gene body correlation, weak positive gene body correlation
# and positive gene body correlation. For each subcluster, rows are clustered by the mean methylation matrix.
#
# There is another plot which shows quantitative statistics in the four groups:
#
# - mean gene body correlation
# - gene length
# - expression difference
# - methylation difference
#
# == value
# An updated ``cr`` that includes the partitioning, this information is important for many downstream analysis.
# Users should update it by ``cr = cr_enriched_heatmap(cr, ...)``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
cr_enriched_heatmap = function(cr, txdb, expr, expr_ha, p_na = 0.25) {

	cr_param = metadata(cr)$cr_param
	sample_id = cr_param$sample_id
	subgroup = cr_param$subgroup
	n_subgroup = length(unique(subgroup))
	subgroup_level = unique(subgroup)
	extend = cr_param$extend

	gm = genes(txdb)
	g = gm[gm$gene_id %in% unique(cr$gene_id)]
	tss = promoters(g, upstream = 1, downstream = 0)

	CHROMOSOME = unique(as.vector(seqnames(g)))
	message("normalize cr to genes")
	mat = normalizeToMatrix(cr, g, 
		mapping_column = "gene_id", value_column = "corr",
	    extend = extend, mean_mode = "absolute", w = 200, target_ratio = 1/3,
	    background = NA)
	
	target_index = attr(mat, "target_index")
	n_target = length(target_index)
	p = apply(mat[, target_index], 1, function(x) sum(is.na(x))/n_target)
	l = p < p_na
	message(qq("there are @{length(l) - sum(l)} genes that have too many NAs (> @{p_na}) in the body."))
	mat2 = mat[l, ]
	mat2[is.na(mat2)] = 0

	cr_gi = rownames(mat2)
	g = g[cr_gi]

	w = width(g)
	names(w) = cr_gi
	w[w > quantile(w, 0.8)] = quantile(w, 0.8)

	target_index = attr(mat2, "target_index")

	mat_target = mat2
	
	message("apply kmeans clustering with 4 groups")
	km4 = kmeans(mat_target, centers = 4)$cluster
	x = tapply(rowMeans(mat_target), km4, mean)
	od = structure(rank(x), names = names(x))
	km4 = structure(od[as.character(km4)], names = names(km4))

	expr = expr[cr_gi, sample_id, drop = FALSE]
	expr = t(apply(expr, 1, function(x) {
		qu = quantile(x, c(0.05, 0.95), na.rm = TRUE)
	    x[x < qu[1]] = qu[1]
	    x[x > qu[2]] = qu[2]
	    x
	}))
	expr = t(scale(t(expr)))

	expr_split = NULL
	if(n_subgroup == 2) {
		subgroup1_sample_id = sample_id[subgroup == subgroup_level[1]]
		subgroup2_sample_id = sample_id[subgroup == subgroup_level[2]]
		expr_mean = rowMeans(expr[, subgroup1_sample_id]) - 
					rowMeans(expr[, subgroup2_sample_id])
		expr_split = ifelse(expr_mean > 0, "high", "low")
		expr_split = factor(expr_split, levels = c("high", "low"))
	}

	message("normalize methylation to genes")
	meth_mat = enrich_with_methylation(g, sample_id, target_ratio = 1/3, extend = extend, w = 200)

	message("normalize methylation difference to genes")
	meth_diff_column = ifelse(n_subgroup == 0 | n_subgroup == 1, "meth_IQR", ifelse(n_subgroup == 2, "meth_diff", "meth_diameter"))
	meth_mat_diff = normalizeToMatrix(cr, g, 
		mapping_column = "gene_id", value_column = meth_diff_column,
	    extend = extend, mean_mode = "absolute", w = 200, target_ratio = 1/3,
	    background = 0)
	meth_mat_diff[is.na(meth_mat_diff)] = 0
	
	if(n_subgroup == 2) {
		combined_split = paste(km4, expr_split, sep = ",")
	} else {
		combined_split = as.character(km4)
	}

	row_order = NULL
	if(n_subgroup == 2) {
		for(i1 in 1:4) {
			for(i2 in c("high", "low")) {
				label = paste(i1, i2, sep = ",")
				l = combined_split == label
				message(qq("cluster @{label}, @{sum(l)} rows"))
				if(sum(l) == 1) {
					row_order = c(row_order, which(l))
				} else {
					dend1 = as.dendrogram(hclust(dist(meth_mat[l, ])))
					dend1 = reorder(dend1, rowMeans(meth_mat[l, ]))
					row_od1 = order.dendrogram(dend1)

					row_order = c(row_order, which(l)[row_od1])
				}
			}
		}
	} else {
		for(i1 in 1:4) {
			l = combined_split == i1
			message(qq("cluster km = @{i1}, @{sum(l)} rows"))
			if(sum(l) == 1) {
				row_order = c(row_order, which(l))
			} else {
				dend1 = as.dendrogram(hclust(dist(meth_mat[l, ])))
				dend1 = reorder(dend1, rowMeans(meth_mat[l, ]))
				row_od1 = order.dendrogram(dend1)

				row_order = c(row_order, which(l)[row_od1])
			}
		}
	}

	expr_col_od = do.call("c", lapply(subgroup_level, function(le) {
		dend1 = as.dendrogram(hclust(dist(t(expr[, subgroup == le, drop = FALSE]))))
		hc1 = as.hclust(reorder(dend1, colMeans(expr[, subgroup == le, drop = FALSE])))
		col_od1 = hc1$order
		which(subgroup == le)[col_od1]
	}))

	group_mean_col = structure(brewer.pal(9, "Set1")[c(3,4,5,1)], names = c(1:4))
	if(missing(expr_ha)) {
		if(n_subgroup >= 2) {
			expr_ha = HeatmapAnnotation(subgroup = subgroup, col = list(subgroup = structure(rand_color(n_subgroup), names = subgroup_level)), 
				show_annotation_name = TRUE, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 10))
		}
	}

	message("making heatmaps")
	ht_list = Heatmap(km4, name = "cluster", col = group_mean_col, show_row_names = FALSE, width = 2*grobHeight(textGrob("1", gp = gpar(fontsize = 10))))
	if(n_subgroup == 2) {
		ht_list = ht_list + Heatmap(expr_split, name = "expr_direction", show_row_names = FALSE, width = unit(5, "mm"), col = c("high" = "Red", "low" = "darkgreen"))
	}
	ht_list = ht_list + EnrichedHeatmap(mat2, name = "correlation", col = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red")),
		use_raster = TRUE, raster_quality = 2, cluster_rows = TRUE, show_row_dend = FALSE, column_title = "Correlation",
		top_annotation = HeatmapAnnotation(enrich = anno_enriched(yaxis_side = "left", gp = gpar(col = group_mean_col))), top_annotation_height = unit(3, "cm"),
		axis_name = c("-50kb", "TSS", "TES", "50kb"), combined_name_fun = NULL) +
	rowAnnotation(gene_length = row_anno_points(w, size = unit(0.5, "mm"), gp = gpar(col = "#00000020"), axis = FALSE), width = unit(2, "cm")) +
	Heatmap(expr, name = "expr", show_row_names = FALSE, col = colorRamp2(quantile(expr, c(0, 0.5, 0.95)), c("blue", "white", "red")),
		use_raster = TRUE, raster_quality = 2,
		show_column_names = FALSE, width = unit(5, "cm"), cluster_columns = FALSE, column_order = expr_col_od,
		top_annotation = expr_ha,
		column_title = "Expression", show_row_dend = FALSE) +
	EnrichedHeatmap(meth_mat, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
		name = "methylation", column_title = qq("Methylation"),
		use_raster = TRUE, raster_quality = 2,
		top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = group_mean_col))),
		axis_name = c("-50kb", "TSS", "TES", "50kb"))

	ht_list = ht_list + EnrichedHeatmap(meth_mat_diff, col = generate_diff_color_fun(meth_mat_diff), 
		name = "methylation_diff", column_title = qq("meth_diff"), axis_name = c("-50kb", "TSS", "TES", "50kb"),
		heatmap_legend_param = list(title = "methylation_diff"),
		top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = group_mean_col))),
		use_raster = TRUE, raster_quality = 2)

	# draw(ht_list, split = km3)
	# decorate_annotation("gene_length", slice = 3, {grid.text("gene_length", 0.5, unit(0, "npc") - unit(1, "cm"))})
	lines_lgd = Legend(at = c("cluster1", "cluster2", "cluster3", "cluster4"), title = "Lines", legend_gp = gpar(col = group_mean_col), type = "lines")

	ht_list = draw(ht_list, row_order = row_order, split = km4, cluster_rows = FALSE, main_heatmap = "correlation", annotation_legend_list = list(lines_lgd))
	decorate_annotation("gene_length", slice = 4, {grid.xaxis(at = c(0, 40000, 80000), label = c("0kp", "40kb", "80kb"), gp = gpar(fontsize = 8));
		grid.text("gene length", 0.5, unit(0, "npc") - unit(1, "cm"), gp = gpar(fontsize = 10))})
	for(i in 1:4) {
		decorate_heatmap_body("cluster", slice = i, {grid.text(i, rot = 90)})
		if(n_subgroup == 2) {
			decorate_heatmap_body("expr_direction", slice = i, {grid.rect(gp = gpar(fill = "transparent"))})
		}
	}

	cr_param$km = km4
	cr_param$group_mean_col = group_mean_col
	cr_param$order = row_order
	cr_param$combined_split = combined_split
	cr_param$expr_split = expr_split
	metadata(cr) = list(cr_param = cr_param)

	if(!dev.interactive()) {
		
		device_type = .Device
		filepath = attr(.Device, "filepath")
		if(device_type == "pdf") {
			filepath = gsub("\\.pdf$", ".stat.pdf", filepath, ignore.case = TRUE)
			stat_plot_width = ifelse(n_subgroup == 2, 12, 9)
			stat_plot_height = 3
			pdf(filepath, width = stat_plot_width, height = stat_plot_height)
		} else if(device_type == "png") {
			filepath = gsub("\\.png$", ".stat.png", filepath, ignore.case = TRUE)
			stat_plot_width = ifelse(n_subgroup == 2, 12, 9)*200
			stat_plot_height = 3*200
			png(filepath, width = stat_plot_width, height = stat_plot_height)
		} else {
			filepath = gsub("\\.[^.]+$", ".stat.png", filepath, ignore.case = TRUE)
			stat_plot_width = ifelse(n_subgroup == 2, 12, 9)*200
			stat_plot_height = 3*200
			png(filepath, width = stat_plot_width, height = stat_plot_height)
		}	
			
		message("making stat plot")
		gene_length = w
		corr_col_fun = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))

		layout(matrix(1:4, nrow = 1), widths = c(1, 1, ifelse(n_subgroup == 2, 2, 1), 1))
		lt = tapply(seq_len(nrow(mat2)), km4, function(ind) mat2[ind, target_index])
		group_mean = sapply(lt, mean)
		boxplot(lt, col = group_mean_col, outline = FALSE, xlab = "cluster",
			ylab = "correlation", main = "gene body correlation")

		lt_gene_length = split(gene_length, km4)
		boxplot(lt_gene_length, col = group_mean_col, outline = FALSE,
			xlab = "cluster", ylab = "gene width", main = "gene length")

		if(n_subgroup == 2) {
			lt1 = split(rowMeans(expr[, subgroup1_sample_id]), combined_split)
			lt2 = split(rowMeans(expr[, subgroup2_sample_id]), combined_split)

			lt = c(lt1, lt2)
			group_mean_col2 = paste0(group_mean_col, "80")
			foo = boxplot(lt, col = c(rep(group_mean_col, each = 2), rep(group_mean_col2, each = 2)), outline = FALSE,
				xlab = NULL, "ylab" = "expression", main = "mean expression", names = rep("", 16))
			par(xpd = NA)
			text(seq_along(lt), -1, paste(names(lt), rep(c("group1", "group2"), 4), sep = ","), adj = c(1, 0.5), srt = 45)
			par(xpd = FALSE)
		} else {
			lt = split(rowMeans(expr), combined_split)

			foo = boxplot(lt, col = group_mean_col, outline = FALSE,
				xlab = NULL, "ylab" = "expression", main = "mean expression", names = rep("", 4))
			par(xpd = NA)
			text(seq_along(lt), -1, paste(names(lt), rep(c("group1", "group2"), 4), sep = ","), adj = c(1, 0.5), srt = 45)
			par(xpd = FALSE)
		}

		lt_diff = split(rowMeans(abs(meth_mat_diff[, target_index])), km4)
		lt = sapply(1:4, function(i) sum(lt_diff[[i]]*lt_gene_length[[i]])/sum(lt_gene_length[[i]]))
		names(lt) = 1:4
		barplot(lt, col = group_mean_col, ylim = c(0, max(lt)*1.1),
			xlab = "cluster", ylab = "diff", main = "weighted abs methy diff\nin gene body")
		box()
		dev.off()
	}

	return2(cr, invisible = TRUE)
}


# == title
# Visualize significantly correlated regions
#
# == param
# -cr correlated regions, should be returned by `cr_enriched_heatmap`.
# -txdb transcriptome annotation which was used in `correlated_regions`
# -fdr_cutoff cutoff for fdr
# -meth_diff_cutoff cutoff for methylation difference. If there are no subgroup information or only one subgroup,
#             ``meth_IQR`` column is used for filtering. If there are more than one subgroups, ``meth_diameter``
#             column is used for filtering.
#
# == details
# There are two heatmaps which corresponds to significant negative correlated regions and positive
# correlated regions. Rows are same as the heatmaps produced by `cr_enriched_heatmap`.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
sig_cr_enriched_heatmap = function(cr, txdb, fdr_cutoff = 0.05, meth_diff_cutoff = 0.1) {

	cr_param = metadata(cr)$cr_param
	subgroup = cr_param$subgroup
	n_subgroup = length(unique(subgroup))

	km4 = cr_param$km
	if(is.null(km4)) {
		stop("`cr` should be returned by `cr_enrichedheatmap()`.")
	}
	group_mean_col = cr_param$group_mean_col
	row_order = cr_param$row_order
	cr_gi = names(km4)

	cr = cr[cr$gene_id %in% cr_gi]

	gm = genes(txdb)
	g = gm[cr_gi]

	if(n_subgroup %in% c(0, 1)) {
		l = cr$corr_fdr < fdr_cutoff & cr$meth_IQR > meth_diff_cutoff
	} else {
		l = cr$corr_fdr < fdr_cutoff & cr$meth_anova_fdr < fdr_cutoff & abs(cr$meth_diameter) > meth_diff_cutoff
	}
	
	cr2 = cr[l]
	mat_neg_cr = normalizeToMatrix(cr2[cr2$corr < 0], g, mapping_column = "gene_id",
	    extend = 50000, mean_mode = "absolute", w = 200, target_ratio = 1/3, background = 0)
	mat_pos_cr = normalizeToMatrix(cr2[cr2$corr > 0], g, mapping_column = "gene_id",
	    extend = 50000, mean_mode = "absolute", w = 200, target_ratio = 1/3, background = 0)
	
	ymax = max(c(tapply(seq_len(nrow(mat_neg_cr)), km4, function(ind) max(colMeans(mat_neg_cr[ind, ]))), 
		  tapply(seq_len(nrow(mat_pos_cr)), km4, function(ind) max(colMeans(mat_pos_cr[ind, ])))))

	ht = Heatmap(km4, name = "cluster", col = group_mean_col, show_row_names = FALSE, width = 2*grobHeight(textGrob("1", gp = gpar(fontsize = 10)))) +
	EnrichedHeatmap(mat_neg_cr, name = "neg_cr", col = c("0" = "white", "1" = "darkgreen"),
		use_raster = TRUE, raster_quality = 4, cluster_rows = TRUE, show_row_dend = FALSE,
		top_annotation = HeatmapAnnotation(enrich = anno_enriched(ylim = c(0, ymax), gp = gpar(col = group_mean_col))), top_annotation_height = unit(3, "cm"),
		axis_name = c("-50kb", "TSS", "TES", "50kb"), combined_name_fun = NULL) +
	EnrichedHeatmap(mat_pos_cr, name = "pos_cr", col = c("0" = "white", "1" = "red"),
		use_raster = TRUE, raster_quality = 4, cluster_rows = TRUE, show_row_dend = FALSE,
		top_annotation = HeatmapAnnotation(enrich = anno_enriched(ylim = c(0, ymax), gp = gpar(col = group_mean_col))), top_annotation_height = unit(3, "cm"),
		axis_name = c("-50kb", "TSS", "TES", "50kb"))
	draw(ht, row_order = row_order, split = km4, cluster_rows = FALSE, main_heatmap = "neg_cr", column_title = qq("fdr_cutoff_@{fdr_cutoff}_meth_diff_cutoff_@{meth_diff_cutoff}"))
	for(i in 1:4) {
		decorate_heatmap_body("cluster", slice = i, {grid.text(i, rot = 90)})
	}
}


# == title
# Compare cutoff for determining significant correlated regions
#
# == param
# -cr correlated regions, should be returned by `cr_enriched_heatmap`
# -txdb transcriptome annotation which was used in `correlated_regions`
# -fdr_cutoff a list of cutoffs to compare
# -meth_diff_cutoff  a list of cutoffs to compare. If there are no subgroup information or only one subgroup,
#             ``meth_IQR`` column is used for filtering. If there are more than one subgroups, ``meth_diameter``
#             column is used for filtering.
#
# == details
# It simply plots how correlated signals are enriched at extended gene regions for
# negative correlated regions and positive correlated regions under different combination
# of cutoffs. The basic rule for selecting the cutoffs is they should not be too loose while
# the pattern of the enrichment (on tss/gene body) should still be obvious.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
sig_cr_compare_cutoff = function(cr, txdb, fdr_cutoff = c(0.1, 0.05, 0.01), 
	meth_diff_cutoff = c(0, 0.1, 0.2, 0.3)) {

 	cr_param = metadata(cr)$cr_param
	subgroup = cr_param$subgroup
	n_subgroup = length(unique(subgroup))

	km4 = cr_param$km
	if(is.null(km4)) {
		stop("`cr` should be returned by `cr_enrichedheatmap()`.")
	}
	group_mean_col = cr_param$group_mean_col
	row_order = cr_param$row_order
	cr_gi = names(km4)

	cr = cr[cr$gene_id %in% cr_gi]

	gm = genes(txdb)
	g = gm[cr_gi]

	mat_neg_cr_list = list()
	mat_pos_cr_list = list()

	ymax = 0
	for(i in seq_along(fdr_cutoff)) {
	     for(j in seq_along(meth_diff_cutoff)) {
	     	
	     	message(qq("for cutoff_@{fdr_cutoff[i]}_meandiff_@{meth_diff_cutoff[j]}"))
	     	
	     	if(n_subgroup %in% c(0, 1)) {
				l = cr$corr_fdr < fdr_cutoff[i] & cr$meth_IQR > meth_diff_cutoff[j]
			} else {
				l = cr$corr_fdr < fdr_cutoff[i] & cr$meth_anova_fdr < fdr_cutoff[i] & abs(cr$meth_diameter) > meth_diff_cutoff[j]
			}

			cr2 = cr[l]
			mat_neg_cr = normalizeToMatrix(cr2[cr2$corr < 0], g, mapping_column = "gene_id",
			    extend = 50000, mean_mode = "absolute", w = 200, target_ratio = 1/3, background = 0)
			mat_pos_cr = normalizeToMatrix(cr2[cr2$corr > 0], g, mapping_column = "gene_id",
			    extend = 50000, mean_mode = "absolute", w = 200, target_ratio = 1/3, background = 0)
			
			mat_neg_cr_list[[qq("fdr_cutoff_@{fdr_cutoff[i]}_meth_diff_cutoff_@{meth_diff_cutoff[j]}")]] = mat_neg_cr
			mat_pos_cr_list[[qq("fdr_cutoff_@{fdr_cutoff[i]}_meth_diff_cutoff_@{meth_diff_cutoff[j]}")]] = mat_pos_cr
			line_list_neg = tapply(seq_len(nrow(mat_neg_cr)), km4, function(ind) colMeans(mat_neg_cr[ind, ]))
			line_list_pos = tapply(seq_len(nrow(mat_pos_cr)), km4, function(ind) colMeans(mat_pos_cr[ind, ]))
			ymax = max(c(ymax, unlist(line_list_pos), unlist(line_list_neg)))
		}
	}

	n_meth_diff = length(meth_diff_cutoff)
	n_fdr = length(fdr_cutoff)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = n_meth_diff+3, ncol = n_fdr+2, 
		heights = unit.c(3*grobHeight(textGrob("A")), 3*grobHeight(textGrob("A")), unit(rep(1, n_meth_diff), "null"), unit(1, "cm")),
		widths = unit.c(3*grobHeight(textGrob("A")), unit(rep(1, n_fdr), "null"), unit(1.5, "cm")))))
	for(i in seq_along(fdr_cutoff)) {
	     for(j in seq_along(meth_diff_cutoff)) {
	     	message(qq("add enrich lines for fdr_cutoff_@{fdr_cutoff[i]}_meth_diff_cutoff_@{meth_diff_cutoff[j]}"))
			
			mat_neg_cr = mat_neg_cr_list[[qq("fdr_cutoff_@{fdr_cutoff[i]}_meth_diff_cutoff_@{meth_diff_cutoff[j]}")]]
			mat_pos_cr = mat_pos_cr_list[[qq("fdr_cutoff_@{fdr_cutoff[i]}_meth_diff_cutoff_@{meth_diff_cutoff[j]}")]]

			upstream_index = attr(mat_neg_cr, "upstream_index")
	        downstream_index = attr(mat_neg_cr, "downstream_index")
	        target_index = attr(mat_neg_cr, "target_index")
	        n1 = length(upstream_index)
	        n2 = length(target_index)
	        n3 = length(downstream_index)

			line_list_neg = tapply(seq_len(nrow(mat_neg_cr)), km4, function(ind) colMeans(mat_neg_cr[ind, ]))
			line_list_pos = tapply(seq_len(nrow(mat_pos_cr)), km4, function(ind) colMeans(mat_pos_cr[ind, ]))
			n = length(line_list_neg[[1]])
			ylim = c(0, ymax*1.1)

			pushViewport(viewport(layout.pos.row = j+2, layout.pos.col = i+1))
			pushViewport(viewport(x = 0, y = 1, width = unit(1, "npc") - unit(2, "mm"), height = unit(1, "npc") - unit(2, "mm"),
				just = c("left", "top")))
			pushViewport(viewport(x = 0, y = 0, width = unit(0.5, "npc") - unit(1, "mm"),
				just = c("left", "bottom"), xscale = c(0, n), yscale = ylim))
			for(k in seq_along(line_list_neg)) {
				grid.lines(rep((n1 - 0.5)/n, 2), c(0, 1), gp = gpar(lty = 2, col = "grey"))
	            grid.lines(rep((n1 + n2 - 0.5)/n, 2), c(0, 1), gp = gpar(lty = 2, col = "grey"))
				grid.lines(seq_len(n) - 0.5, line_list_neg[[k]], default.units = "native", gp = gpar(col = group_mean_col[k]))
			}
			if(j == 4) {
				grid.text(c("TSS", "TES"), x = c(n1, n1+n2), y = unit(-1, "mm"), default.units = "native", just = "right", rot = 90, gp = gpar(fontsize = 8))
				grid.text(c("-50kb"), x = 0, y = unit(-1, "mm"), default.units = "native", just = c("right", "top"), rot = 90, gp = gpar(fontsize = 8))
				grid.text(c("50kb"), x = n, y = unit(-1, "mm"), default.units = "native", just = c("right", "bottom"), rot = 90, gp = gpar(fontsize = 8))
			}
			grid.rect(gp = gpar(fill = "transparent"))
			upViewport()

			pushViewport(viewport(x = unit(0.5, "npc") + unit(1, "mm"), y = 0, width = unit(0.5, "npc") - unit(1, "mm"), 
				just = c("left", "bottom"), xscale = c(0, n), yscale = ylim))
			for(k in seq_along(line_list_pos)) {
				grid.lines(rep((n1 - 0.5)/n, 2), c(0, 1), gp = gpar(lty = 2, col = "grey"))
	            grid.lines(rep((n1 + n2 - 0.5)/n, 2), c(0, 1), gp = gpar(lty = 2, col = "grey"))
				grid.lines(seq_len(n) - 0.5, line_list_pos[[k]], default.units = "native", gp = gpar(col = group_mean_col[k]))
			}
			if(j == 4) {
				grid.text(c("TSS", "TES"), x = c(n1, n1+n2), y = unit(-1, "mm"), default.units = "native", just = "right", rot = 90, gp = gpar(fontsize = 8))
				grid.text(c("-50kb"), x = 0, y = unit(-1, "mm"), default.units = "native", just = c("right", "top"), rot = 90, gp = gpar(fontsize = 8))
				grid.text(c("50kb"), x = n, y = unit(-1, "mm"), default.units = "native", just = c("right", "bottom"), rot = 90, gp = gpar(fontsize = 8))
			}

			if(i == 3) {
				grid.yaxis(main = FALSE, gp = gpar(fontsize = 8))
			}
			# grid.yaxis(main = FALSE, gp = gpar(fontsize = 8))
			grid.rect(gp = gpar(fill = "transparent"))
			upViewport()
			upViewport()
			upViewport()
		}
	}
	for(i in seq_along(fdr_cutoff)) {
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i + 1))
		grid.text(qq("fdr < @{fdr_cutoff[i]}"))
		upViewport()

		pushViewport(viewport(layout.pos.row = 2, layout.pos.col = i + 1))
		grid.text("negCR", x = 0.25)
		grid.text("posCR", x = 0.75)
		upViewport()
	}
	for(j in seq_along(meth_diff_cutoff)) {
		pushViewport(viewport(layout.pos.row = j+2, layout.pos.col = 1))
		grid.text(qq("meth_diff > @{meth_diff_cutoff[j]}"), rot = 90)
		upViewport()
	}
	pushViewport(viewport(layout.pos.row = 3:6, layout.pos.col = 5))
	grid.text("proportion", rot = 90, x = unit(1, "cm"))
	upViewport()
}

# == title
# Visualize CR genes in gtrellis layout
#
# == param
# -cr correlated regions, should be returned by `cr_enriched_heatmap`
# -txdb transcriptome annotation which was used in `correlated_regions`
# -expr expression matrix which was used in `correlated_regions`
#
# == details
# CR genes in k-means group 1 and 4 (which correspond to negative correlated gene body
# and positive correlated gene body) are visualized in gtrellis layout. Cytobands
# which are significantly overlapped by CR genes are highlighted. In gtrellis layout, there
# are following tracks:
#
# - rainfall plot for genes in k-means cluster 1. The y-axis corresponds to the minimal distance
#   to neighbouring genes. The color of points corresponds to the epxression value (blue is low expression
#   and red is high expression) and size of points corresponds to gene length.
# - a one row heatmap showing how much each cytoband is covered by CR genes in cluster 1.
# - rainfall plot for genes in k-means cluster 4
# - a one row heatmap showing how much each cytoband is covered by CR genes in cluster 4.
# - cytoband
#
# == value
# A list of two elements which shows how each cytoband is overlapped by CR genes
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
cr_genes_gtrellis = function(cr, txdb, expr) {
	
	cr_param = metadata(cr)$cr_param
	subgroup = cr_param$subgroup
	n_subgroup = length(unique(subgroup))
	species = cr_param$genome

	km4 = cr_param$km
	if(is.null(km4)) {
		stop("`cr` should be returned by `cr_enrichedheatmap()`.")
	}
	group_mean_col = cr_param$group_mean_col
	row_order = cr_param$row_order
	cr_gi = names(km4)

	cr = cr[cr$gene_id %in% cr_gi]

	gm = genes(txdb)
	g = gm[cr_gi]
	w = width(gm)
	names(w) = names(gm)

	den_list = lapply(1:4, function(i) {
		gi = names(km4[km4 == i])
		g2 = g[g$gene_id %in% gi]
		genomicDensity(as.data.frame(g2), window.size = 5e6)
	})
	all_max = sapply(den_list, function(den) max(den[[4]]))
	# den_col_fun = colorRamp2(seq(0, max, length = 11), rev(brewer.pal(11, "RdYlBu")))

	len_range = quantile(width(g), c(0, 0.9))
	len_fun = function(x, size = c(0, 1)) {
		x = (size[2] - size[1])/(len_range[2] - len_range[1])*(x - len_range[1]) + size[1]
		x[x < size[1]] = size[1]
		x[x > size[2]] = size[2]
		x
	}

	expr_fun = colorRamp2(quantile(expr, c(0, 0.5, 0.95)), c("blue", "white", "red"), transparency = 0.2)

	rainfall_list = lapply(1:4, function(i) {
		gi = names(km4[km4 == i])
		g2 = g[g$gene_id %in% gi]
		df = rainfallTransform(as.data.frame(g2))
		df
	})

	cytoband = read.cytoband(species = species)$df
	cytoband = GRanges(seqnames = cytoband[[1]], ranges = IRanges(cytoband[[2]], cytoband[[3]]), band = cytoband[[4]])
	cytoband = annotate_to_genomic_features(cytoband, g, name = "all_genes", prefix = "percent_")
	for(i in 1:4) {
		gi = names(km4[km4 == i])
		g2 = g[g$gene_id %in% gi]

		cytoband = annotate_to_genomic_features(cytoband, g2, name = paste0("cluster", i), prefix = "percent_")
		x = mcols(cytoband)[, paste0("percent_cluster", i)]/mcols(cytoband)[, "percent_all_genes"]
		x[is.na(x)] = 0
		mcols(cytoband)[, paste0("percent_cluster", i)] = x
	}
	for(i in 1:4) {
		gi = names(km4[km4 == i])
		g2 = g[g$gene_id %in% gi]

		cytoband = annotate_to_genomic_features(cytoband, g2, type = "number", name = paste0("cluster", i), prefix = "number_")
	}

	gtrellis_layout(n_track = 4, track_ylim = c(-0.5, 7, 0, 1, -0.5, 7, 0, 1), 
		nrow = 3, compact = TRUE, track_ylab = c("rainfall", "", "rainfall", ""),
		track_axis = c(TRUE, FALSE, TRUE, FALSE),
		track_height = rep(unit.c(unit(1, "null"), unit(2, "mm")), 2),
		add_ideogram_track = TRUE, add_name_track = TRUE, category = paste0("chr", 1:22))

	gene_list = list()
	cytoband_list = list()
	cytoband_random_list = list()
	for(i in c(1,4)) {
		# rainfall track
		y = log10(rainfall_list[[i]]$dist+1)
		y[y == 0] = y[y==0] + runif(sum(y == 0), min = -0.3, max = 0.3)
		ggg = rownames(rainfall_list[[i]])
		add_points_track(rainfall_list[[i]], y, size = unit(len_fun(w[ggg], size = c(1, 4)), "mm"),
			gp = gpar(col = expr_fun(rowMeans(expr)[ggg])))
		
		# percent track
		l = rep(TRUE, length(cytoband))
		for(j in setdiff(c(1,4), i)) {
			l = l & mcols(cytoband)[, paste0("percent_cluster", i)] - mcols(cytoband)[, paste0("percent_cluster", j)] > 0.5
		}
		l = l & mcols(cytoband)[, "percent_all_genes"] > 0.4
		l[is.na(l)] = FALSE
		add_rect_track(cytoband[l], unit(0, "npc"), unit(1, "npc"), gp = gpar(col = NA, fill = "#00FF0040"), track = current_track())
		
		mtch = as.matrix(findOverlaps(cytoband[l], g))
		gene_list[[as.character(i)]] = names(g[unique(mtch[, 2])])
		cytoband_list[[as.character(i)]] = cytoband[l]
		
		# cytoband
		add_track(cytoband[l], track = current_track(), panel_fun = function(gr) {
			n_gr = length(gr)
			n = mcols(gr)[, paste0("number_cluster", i)]
			p = mcols(gr)[, paste0("percent_cluster", i)]
			p0 = mcols(gr)[, "percent_all_genes"]
			is_neighbouring = abs(end(gr)[-n_gr] - start(gr)[-1]) < 2
			is_neighbouring = c(is_neighbouring, FALSE)
			number_mat = as.matrix(mcols(gr)[, paste0("number_cluster", 1:4)])
			# grid.text(qq("@{n}/@{rowSums(number_mat)}, @{sprintf('%.2f', p)},@{sprintf('%.2f', p0)}", collapse = FALSE), x = (start(gr) + end(gr))/2, y = unit(0.1, "npc"), just = "left",
			grid.text(qq("@{sprintf('%.2f', p)}", collapse = FALSE), x = (start(gr) + end(gr))/2, y = unit(ifelse(is_neighbouring, 0.35, 0.1), "npc"), just = "left",
			rot = 90, default.units = "native", gp = gpar(fontsize = 8))
		}, clip = FALSE)

		den_col_fun = colorRamp2(seq(0, 1, length = 11), rev(brewer.pal(11, "RdYlBu")))

		add_track(NULL, clip = FALSE, category = "chr22", track = current_track(), panel_fun = function(gr) {
			oxlim = get_cell_meta_data("extended_xlim")
			oylim = get_cell_meta_data("extended_ylim")
			grid.text(qq("#gene: @{sum(km4 == i)}"), unit(oxlim[2], "native") + unit(2, "cm"), mean(oylim), default.units = "native", just = "left")
		})

		add_heatmap_track(cytoband, as.matrix(mcols(cytoband)[, paste0("percent_cluster", i)]), fill = den_col_fun, border = "black")
		# add_lines_track(den_list[[i]], den_list[[i]][[4]], area = TRUE, gp = gpar(fill = "#B0B0B0"))
	}
	return(invisible(cytoband_list))
}

# == title
# Visualize correlations in cytoband
#
# == param
# -cr correlated regions returned from `cr_enriched_heatmap`
# -txdb transcriptome annotation which was used in `correlated_regions`
# -cytoband_list a list of cytoband returned by `cr_genes_gtrellis`
# -color_head internal use
#
# == details
# This function visualizes significant cytobands which have been found in `cr_genes_gtrellis`.
# 
# For each cytoband, there are several tracks:
#
# - cytoband name. Green corresponds to significant cytobands in group 1 and red corresponds to group 4
# - points showing correlations
# - a one row heamtap showing mean correlation in 50kb window
# - genes, green lines are genes in cluster 1, red lines are genes in cluster 4.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
sig_cytoband_gtrellis = function(cr, txdb, cytoband_list, color_head = TRUE) {
	
	cr_param = metadata(cr)$cr_param
	subgroup = cr_param$subgroup
	n_subgroup = length(unique(subgroup))

	km4 = cr_param$km
	if(is.null(km4)) {
		stop("`cr` should be returned by `cr_enrichedheatmap()`.")
	}
	group_mean_col = cr_param$group_mean_col
	row_order = cr_param$row_order
	cr_gi = names(km4)

	cr = cr[cr$gene_id %in% cr_gi]

	gm = genes(txdb)
	g = gm[cr_gi]

	cytoband_list = lapply(cytoband_list, function(ct) {
		names(ct) = paste(as.vector(seqnames(ct)), ct$band, sep = ".")
		ct
	})

	ct = c(cytoband_list[["1"]], cytoband_list[["4"]])
	ct$group = c(rep("group1", length(cytoband_list[["1"]])), rep("group4", length(cytoband_list[["4"]])))
	ct = sort(ct)
	# ct = ct[!duplicated(names(ct))]
	gl = width(g)
	names(gl) = names(g)
	mtch1 = as.matrix(findOverlaps(ct, cr))
	mtch2 = as.matrix(findOverlaps(ct, g))
	ct_df = as.data.frame(ct)
	ct_df[[1]] = rownames(ct_df)

	col_fun = colorRamp2(c(-0.8, 0, 0.8), c("green", "white", "red"))
	gtrellis_layout(data = ct_df, n_track = 3, track_ylim = c(-1, 1, 0, 1, 0, 1), 
		nrow = 6, compact = TRUE, track_ylab = c("correlation", "", "gene"),
		track_axis = c(TRUE, FALSE, FALSE),
		track_height = unit.c(unit(1, "null"), unit(2, "mm"), unit(1, "cm")),
		add_name_track = TRUE, name_track_fill = NA)
	ct_name = names(ct)
	for(i in seq_len(length(ct))) {
		ind = unique(mtch1[mtch1[, 1] == i, 2])
		gr = cr[ind]
		l = gr$gene_tss_dist > 0 & gr$gene_tss_dist <= gl[gr$gene_id]
		l[is.na(l)] = FALSE
		gr = gr[l]
		gr = GRanges(seqnames = ct_name[i], ranges = ranges(gr), corr = gr$corr)
		add_points_track(gr, gr$corr, category = ct_name[i], size = unit(0.5, "mm"), track = 2, gp = gpar(col = col_fun(gr$corr)))
		
		ct_foo = GRanges(seqnames = ct_name[i], ranges = ranges(ct[i]))
		wd = makeWindows(ct_foo, w = 50000)
		mat = normalizeToMatrix(gr, ct_foo, k = length(wd), value_column = "corr", extend = 0)
		add_heatmap_track(wd, as.vector(mat), fill = col_fun, category = ct_name[i], track = 3)
		add_track(NULL, category = ct_name[i], track = 3, panel_fun = function(gr) grid.rect(gp = gpar(fill = "transparent", col = "black")))

		ind = unique(mtch2[mtch2[, 1] == i, 2])
		gr = g[ind]
		gr = GRanges(seqnames = ct_name[i], ranges = ranges(gr), gene_id = names(gr))
		if(color_head) {
			segment_col = ifelse(gr$gene_id %in% names(km4[km4 == 4]), "red", ifelse(gr$gene_id %in% names(km4[km4 == 1]), "green", "#AAAAAA"))
		} else {
			segment_col = "black"
		}

		add_segments_track(gr, runif(length(gr), min = 0.05, max = 0.95), category = ct_name[i], track = 4,
			gp = gpar(col = segment_col, lwd = ifelse(gr$gene_id %in% names(km4[km4 %in% c(1, 4)]), 2, 1)))

		if(color_head) {
			add_track(NULL, category = ct_name[i], track = 1, panel_fun = function(gr) {
				grid.rect(gp = gpar(fill = ifelse(ct_df[i, "group"] == "group4", "#FF000040", "#00FF0040")))
			})
		}
	}
}

# == title
# CR coverage on genes
#
# == param
# -cr correlated regions
# -txdb transcriptome annotation which was used in `correlated_regions`
#
cr_coverage_on_genes = function(cr, txdb) {
	
	cr_param = metadata(cr)$cr_param
	extend = cr_param$extend

	gm = genes(txdb)
	gm_start = start(gm)
	gm_end = end(gm)
	names(gm_start) = names(gm)
	names(gm_end) = names(gm)


	cr_list = split(cr, cr$gene_id)
	if(length(cr_list) > 2000) {
		cr_list = cr_list[sample(seq_len(length(cr_list)), 2000)]
	}

	pct_list = lapply(cr_list, function(gr) {
		gi = gr$gene_id[1]
		gr = reduce(ranges(gr))
		start = start(gr)
		end = end(gr)

		g_start = gm_start[gi] - extend
		g_end = gm_end[gi] + extend
		g_range = g_end - g_start + 1
		start = (start - g_start + 1)/g_range
		end = (end - g_start)/g_range
		start[start < 0] = 0
		data.frame(start = start, end = end)
	})

	pct = sapply(pct_list, function(x) sum(x$end - x$start))
	od = order(pct)
	pct_list = pct_list[od]
	pct = pct[od]
	gene_length = gm_end[names(pct_list)] - gm_start[names(pct_list)]

	for(i in seq_along(pct_list)) {
		pct_list[[i]] = cbind(pct_list[[i]], i = i)
	}
	pct_df = do.call("rbind", pct_list)

	n = length(pct_list)
	grid.newpage()
	pushViewport(viewport(x = unit(5, "mm"), width = unit(1, "npc") - unit(5, "cm"), 
		height = 0.85, xscale = c(0, 1), yscale = c(0, n),
		just = "left"))
	grid.rect(x = pct_df$start, y = pct_df$i - 1, width = pct_df$end - pct_df$start, height = 1, 
		default.units = "native", gp = gpar(fill = "black", col = NA), just = c("left", "bottom"))
	grid.xaxis(at = c(0, 1), label = c(-extend, extend))
	upViewport()

	pushViewport(viewport(x = unit(1, "npc") - unit(5, "mm"), width = unit(3.5, "cm"), height = 0.85,
		just = "right", xscale = c(0, 1), yscale = c(0, n)))
	grid.lines(pct, seq_len(n), default.units = "native")
	grid.xaxis()
	grid.rect(gp = gpar(fill = "transparent"))
	upViewport() 

	fake_mat = matrix(0, nrow = n, ncol = 10)

	ht = Heatmap(fake_mat, name = "fake", col = c("0" = "white"), cluster_rows = FALSE, cluster_columns = FALSE, 
		rect_gp = gpar(type = "none"), show_heatmap_legend = FALSE, column_title = "CR coverage on genes") +
	rowAnnotation(pct = row_anno_points(pct, axis = TRUE), width = unit(3, "cm")) +
	rowAnnotation(len = row_anno_points(gene_length, axis = TRUE), width = unit(3, "cm"))
	draw(ht, padding = unit(c(2, 1, 0.5, 0.5), "cm"))
	decorate_heatmap_body("fake", {
		pushViewport(viewport(xscale = c(0, 1), yscale = c(0, n)))
		grid.rect(x = pct_df$start, y = n - pct_df$i - 1, width = pct_df$end - pct_df$start, height = 1, 
			default.units = "native", gp = gpar(fill = "black", col = NA), just = c("left", "bottom"))
		grid.xaxis(at = c(0, 1), label = c(-extend, extend))
		upViewport()
	})
	decorate_annotation("pct", {grid.text("Percent", x = 0.5, y = unit(1, "npc") + unit(5, "mm"))})
	decorate_annotation("len", {grid.text("Gene length", x = 0.5, y = unit(1, "npc") + unit(5, "mm"))})
}
