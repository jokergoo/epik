normalize_epigenomic_signals = function(cr, target, marks = NULL, expr = NULL, include_correlation_matrix = TRUE,
	extend = 5000, target_ratio = 0.1) {

	cr_param = metadata(cr)$cr_param
	sample_id = cr_param$sample_id
	subgroup = cr_param$subgroup
	subgroup_level = unique(subgroup)
	n_subgroup = length(subgroup_level)

	target_name = deparse(substitute(target))

	############ enriched to methylations ##################
	message(qq("normalizing methylation to @{target_name}"))

	if(include_correlation_matrix) {
		meth_mat_corr = normalizeToMatrix(cr, target, mapping_column = "gene_id", value_column = "corr",
			extend = extend, mean_mode = "absolute", w = 50, target_ratio = target_ratio, trim = 0, empty_value = 0)
	} else {
		meth_mat_corr = NULL
	}

	meth_mat_mean = enrich_with_methylation(target, sample_id, target_ratio = target_ratio, extend = extend)
	meth_mat_mean[attr(meth_mat_mean, "failed_rows"), ] = 0.5

	if(n_subgroup <= 1) {
		meth_mat_diff = enrich_with_methylation(target, sample_id, target_ratio = target_ratio, extend = extend, mode = rowIQRs)
		meth_mat_diff[attr(meth_mat_diff, "failed_rows"), ] = 0
	} else if(n_subgroup == 2) {
		meth_mat_mean_1 = enrich_with_methylation(target, sample_id[subgroup == subgroup_level[1]], target_ratio = target_ratio, extend = extend)
		failed_rows = attr(meth_mat_mean_1, "failed_rows")
		message(qq("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets."))
		meth_mat_mean_1[failed_rows, ] = 0.5

		meth_mat_mean_2 = enrich_with_methylation(target, sample_id[subgroup == subgroup_level[2]], target_ratio = target_ratio, extend = extend)
		failed_rows = attr(meth_mat_mean_2, "failed_rows")
		message(qq("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets."))
		meth_mat_mean_2[failed_rows, ] = 0.5
		meth_mat_diff = meth_mat_mean_1 - meth_mat_mean_2
	} else {
		meth_mat_mean_list = lapply(subgroup_level, function(le) {
			meth_mat_mean = enrich_with_methylation(target, sample_id[subgroup == le], target_ratio = target_ratio, extend = extend)
			failed_rows = attr(meth_mat_mean, "failed_rows")
			message(qq("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets."))
			meth_mat_mean[failed_rows, ] = 0.5
			meth_mat_mean
		})
		meth_array_mean = array(dim = c(dim(meth_mat_mean), n_subgroup))
		for(i in seq_along(meth_mat_mean_list)) {
			meth_array_mean[, , i] = meth_mat_mean_list[[i]]
		}
		meth_mat_diff = apply(meth_array_mean, c(1, 2), max)  - apply(meth_array_mean, c(1, 2), min)
		meth_mat_diff = copyAttr(meth_mat_mean, meth_mat_diff)
	}


	################# enrich to histome modifications #################
	hist_mat_corr_list = list()
	hist_mat_mean_list = list()
	hist_mat_diff_list = list()
	
	for(k in seq_along(marks)) {
		
		message(qq("normalizing @{marks[k]} signals to @{target_name}"))

		hm_sample_id = intersect(sample_id, chipseq_hooks$sample_id(marks[k]))
		hm_subgroup = subgroup[sample_id %in% hm_sample_id]
		hm_subgroup_level = unique(hm_subgroup)
		n_hm_subgroup = length(hm_subgroup_level)

		# applied to each sample, each mark
		lt = enrich_with_histone_mark(target, sample_id = hm_sample_id, mark = marks[k], return_arr = TRUE, 
			target_ratio = target_ratio, extend = extend)
		hist_mat_mean_list[[k]] = lt[[2]]
		arr = lt[[1]]

		# only calculate the correlation when there are enough samples
		if(length(hm_sample_id) >= 5 && include_correlation_matrix) {
			# detect regions that histone MARKS correlate to expression
			expr2 = expr[target$gene_id, intersect(colnames(expr), hm_sample_id)]
			hist_mat_corr = matrix(nrow = nrow(expr2), ncol = ncol(meth_mat_mean))

			counter = set_counter(nrow(hist_mat_corr), fmt = "  calculate correlation for %s rows.")
			for(i in seq_len(nrow(hist_mat_corr))) {
				counter()
			    for(j in seq_len(ncol(hist_mat_corr))) {
			        suppressWarnings(x <- cor(arr[i, j, ], expr2[i, ], method = "spearman"))
			        hist_mat_corr[i, j] = x
			    }
			}
			hist_mat_corr[is.na(hist_mat_corr)] = 0
			hist_mat_corr = copyAttr(meth_mat_mean, hist_mat_corr)
			hist_mat_corr_list[[k]] = hist_mat_corr
		} else {
			hist_mat_corr_list[k] = list(NULL)
		}

		if(n_hm_subgroup <= 1) {
			hist_mat_diff_list[[k]] = apply(arr, c(1, 2), IQR, na.rm = TRUE)
			hist_mat_diff_list[[k]] = copyAttr(meth_mat_mean, hist_mat_diff_list[[k]])
		} else if(n_hm_subgroup == 2) {
			hm_sample_id_subgroup1 = hm_sample_id[hm_subgroup == hm_subgroup_level[1]]
			hm_sample_id_subgroup2 = hm_sample_id[hm_subgroup == hm_subgroup_level[2]]
			h1 = apply(arr[, , hm_sample_id_subgroup1], c(1, 2), mean, na.rm = TRUE)
			h2 = apply(arr[, , hm_sample_id_subgroup2], c(1, 2), mean, na.rm = TRUE)
			hist_mat_diff_list[[k]] = h1 - h2
			hist_mat_diff_list[[k]] = copyAttr(meth_mat_mean, hist_mat_diff_list[[k]])
		} else {
			h_list = lapply(hm_subgroup_level, function(le) {
				apply(arr[, , hm_sample_id[hm_subgroup == le]], c(1, 2), mean, na.rm = TRUE)
			})
			hm_array_mean = array(dim = c(dim(meth_mat_mean), n_hm_subgroup))
			for(i in seq_along(h_list)) {
				hm_array_mean[, , i] = h_list[[i]]
			}
			hist_mat_diff_list[[k]] = apply(hm_array_mean, c(1, 2), max)  - apply(hm_array_mean, c(1, 2), min)
			hist_mat_diff_list[[k]] = copyAttr(meth_mat_mean, hist_mat_diff_list[[k]])
		}
	}
	names(hist_mat_corr_list) = marks
	names(hist_mat_mean_list) = marks
	names(hist_mat_diff_list) = marks

	return(list(meth_mat_corr = meth_mat_corr,
		        meth_mat_mean = meth_mat_mean,
		        meth_mat_diff = meth_mat_diff,
		        hist_mat_corr_list = hist_mat_corr_list,
		        hist_mat_mean_list = hist_mat_mean_list,
		        hist_mat_diff_list = hist_mat_diff_list))
}

if(is.memoised(normalize_epigenomic_signals)) {
	normalize_epigenomic_signals = memoise(normalize_epigenomic_signals)
}

merge_row_order = function(mat, l_list) {
	do.call("c", lapply(l_list, function(l) {
		if(sum(l) == 0) return(integer(0))
		if(sum(l) == 1) return(which(l))
		dend1 = as.dendrogram(hclust(dist_by_closeness2(mat[l, ])))
		dend1 = reorder(dend1, rowMeans(mat[l, ]))
		od = order.dendrogram(dend1)
		which(l)[od]
	}))
}


add_boxplot_as_column_annotation = function(ht_list, width, anno_name, anno_title = anno_name) {
	gl = width

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
		grid.text(anno_title, y = unit(1, "npc") + unit(2.5, "mm"), gp = gpar(fontsize = 14), just = "bottom")
		upViewport()
	})
}

# == title
# Visualizing enrichment for epigenomic signals at TSS-CGIs
#
# == param
# -cr correalted regions
# -txdb transcriptome annotation which was used in `correlated_regions`
# -expr expression matrix which was used in `correlated_regions`
# -cgi CpG island
# -fdr_cutoff cutoff for fdr of correlation p-value and anove p-values
# -meth_diff_cutoff cutoff for methylation difference
# -marks names of histone marks
# -type use negative correlated regions or positive correlated regions
# -extend base pairs extended to upstream and downstream
# -expr_annotation a `ComplexHeatmap::HeatmapAnnotation` class object
#
# == details
# There are several heatmaps visualize various signals enriched at TSS-CGIs.
#
# - heatmap for gene expression
# - If ``cr`` is returned form `cr_enrichedheatmap`, there is a one column heatmap which
#   shows the k-means cluters genes belong to
# - heatmap for correlated regions
# - a point plot 
# 
#
# == value
# no value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
cr_enriched_heatmap_at_cgi = function(cr, txdb, expr, cgi,
	fdr_cutoff = 0.05, meth_diff_cutoff = 0.1, marks = NULL, type = "neg", extend = 5000,
	expr_annotation) {
	
	eval(SNIPPET_ATTACH_CR_PARAM)

	message(qq("filter sigCR by fdr_cutoff = @{fdr_cutoff}, meth_diff_cutoff = @{meth_diff_cutoff}"))
	eval(SNIPPET_FILTER_SIG_CR)
		
	gene = genes(txdb)

	message("extracting gene tss")
	tss = promoters(gene, upstream = 1, downstream = 0)
	tss = tss[names(tss) %in% cr$gene_id]

	if(length(extend) == 1) extend = rep(extend, 2)

	cgi_extend = cgi
	start(cgi_extend) = start(cgi) - extend[1]
	end(cgi_extend) = end(cgi) + extend[2]
	
	# there should only be one tss in +-5kb of CGI and one tss should
	# only overlaps to one extended CGI
	mtch = as.matrix(findOverlaps(cgi_extend, tss))
	t1 = table(mtch[, 1])
	t2 = table(mtch[, 2])
	s1 = as.numeric(names(t1[t1 == 1]))
	s2 = as.numeric(names(t2[t2 == 1]))
	l = mtch[, 1] %in% s1 & mtch[, 2] %in% s2
	mtch = mtch[l, ]
	cgi2 = cgi[mtch[, 1]]
	cgi2$gene_id = names(tss[mtch[, 2]])
	names(cgi2) = cgi2$gene_id
	strand(cgi2) = strand(tss[mtch[, 2]])
	message(qq("@{length(cgi2)} left filtered by one-to-one mapping between tss and cgi"))

	target_ratio = mean(width(cgi2))/(sum(extend) + mean(width(cgi2)))
	target = cgi2
	target_name = "cgi"

	message("normalize to sigCR")
	eval(SNIPPET_NORMALIZE_SIG_CR)

	if(is.not.null(km)) {
		km = km[names(cgi2)]
	}

	cgi_width = width(cgi2)
	cgi_width[cgi_width > quantile(cgi_width, 0.99)] = quantile(cgi_width, 0.99)
	width_anno_name = "cgi_width"
	width_anno = cgi_width

	message("normalize to epi signals")
	eval(SNIPPET_NORMALIZE_EPI_SIGNALS)

	########### prepare the order of rows and columns
	message("determing row and colum orders")
	eval(SNIPPET_ROW_ORDER_AND_COLUMN_ORDER)

	cor_col_fun = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))
	meth_col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
	cr_col = c("-1" = "darkgreen", "0" = "white", "1" = "red")

	fixed_heatmap = 2

	eval(SNIPPET_HEATMAP_PAGE)

	if(missing(expr_annotation)) {
		if(n_subgroup >= 2) {
			expr_annotation = HeatmapAnnotation(subgroup = subgroup, col = list(subgroup = structure(rand_color(n_subgroup), names = subgroup_level)), 
				show_annotation_name = TRUE, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 10))
		}
	}

	######### construct heatmap list  ############
	## if there are too many heatmaps, they will be put in two pages.
	epi_color = c(brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
	epi_mark_list = list()
	epi_title_list = list()

	n_heatmap = 0
	n_row_group = 2
	ht_list = Heatmap(expr, name = "expr", show_row_names = FALSE,
		show_column_names = FALSE, width = unit(5, "cm"), show_column_dend = FALSE, cluster_columns = FALSE, column_order = expr_col_od,
		top_annotation = expr_annotation, column_title = "Expression", show_row_dend = FALSE,
		use_raster = TRUE, raster_quality = 2)
	n_heatmap = n_heatmap + 1
	
	if(is.not.null(km)) {
		ht_list = ht_list + Heatmap(km[names(target)], name = "km_groups", col = km_col, show_row_names = FALSE,
			width = unit(0.5, "cm"))
	}

	ht_list = ht_list + EnrichedHeatmap(mat_mix, name = qq("@{type}CR"), col = cr_col,
	                top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(neg_col = "darkgreen", pos_col = "red", lty = 1:n_row_group))), 
	                top_annotation_height = unit(2, "cm"), column_title = qq("sig@{type}CR"),
	                use_raster = TRUE, raster_quality = 2, combined_name_fun = NULL)
	n_heatmap = n_heatmap + 1

	ht_list = ht_list + rowAnnotation(foo_width = row_anno_points(width_anno, axis = TRUE, gp = gpar(col = "#00000040")),
	                width = unit(1, "cm"))

	message("append epi heatmaps")
	eval(SNIPPET_APPEND_EPI_HEATMAP)

	message("draw heatmaps")
	eval(SNIPPET_DRAW_HEATMAP)

	return(invisible(NULL))
}


# == title
# Visualizing enrichment for epigenomic signals at TSS
#
# == param
# -cr
# -txdb
# -expr
# -cgi
# -fdr_cutoff
# -meth_diff_cutoff
# -marks
# -type
# -extend
# -expr_annotation 
#
# == details
#
# == value
# no value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
cr_enriched_heatmap_at_tss = function(cr, txdb, expr, cgi,
	fdr_cutoff = 0.05, meth_diff_cutoff = 0.1, marks = NULL, type = "neg", extend = c(5000, 10000),
	expr_annotation) {
	
	eval(SNIPPET_ATTACH_CR_PARAM)

	message(qq("filter sigCR by fdr_cutoff = @{fdr_cutoff}, meth_diff_cutoff = @{meth_diff_cutoff}"))
	eval(SNIPPET_FILTER_SIG_CR)
		
	gene = genes(txdb)

	message("extracting gene tss")
	tss = promoters(gene, upstream = 1, downstream = 0)
	tss = tss[names(tss) %in% sig_cr$gene_id]
	target = tss
	axis_name = c("-5KB", "TSS", "5KB")
	target_ratio = 0.1
	target_name = "tss"
	
	message("normalize to sigCR")
	eval(SNIPPET_NORMALIZE_SIG_CR)

	if(is.not.null(km)) {
		km = km[names(target)]
	}

	gl = width(gene[tss$gene_id])
	gl[gl > quantile(gl, 0.95)] = quantile(gl, 0.95)
	width_anno_name = "gene_width"
	width_anno = gl

	# normalize to CGI
	mat_cgi = normalizeToMatrix(cgi, target, extend = extend, mean_mode = "absolute", w = 50)

	message("normalize to epi signals")
	eval(SNIPPET_NORMALIZE_EPI_SIGNALS)

	########### prepare the order of rows and columns
	message("determing row and colum orders")
	eval(SNIPPET_ROW_ORDER_AND_COLUMN_ORDER)

	cor_col_fun = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))
	meth_col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
	cr_col = c("-1" = "darkgreen", "0" = "white", "1" = "red")

	fixed_heatmap = 3
	eval(SNIPPET_HEATMAP_PAGE)

	if(missing(expr_annotation)) {
		if(n_subgroup >= 2) {
			expr_annotation = HeatmapAnnotation(subgroup = subgroup, col = list(subgroup = structure(rand_color(n_subgroup), names = subgroup_level)), 
				show_annotation_name = TRUE, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 10))
		}
	}

	######### construct heatmap list  ############
	## if there are too many heatmaps, they will be put in two pages.
	epi_color = c(brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
	epi_mark_list = list()
	epi_title_list = list()

	n_heatmap = 0
	n_row_group = 2
	ht_list = Heatmap(expr, name = "expr", show_row_names = FALSE,
		show_column_names = FALSE, width = unit(5, "cm"), show_column_dend = FALSE, cluster_columns = FALSE, column_order = expr_col_od,
		top_annotation = expr_annotation, column_title = "Expression", show_row_dend = FALSE,
		use_raster = TRUE, raster_quality = 2)
	n_heatmap = n_heatmap + 1
	
	ht_list = ht_list + rowAnnotation(foo_width = row_anno_points(width_anno, axis = TRUE, gp = gpar(col = "#00000040")),
	                width = unit(1, "cm"))

	if(is.not.null(km)) {
		ht_list = ht_list + Heatmap(km[names(tss)], name = "km_groups", col = km_col, show_row_names = FALSE,
			width = unit(1, "cm"))
	}

	ht_list = ht_list + EnrichedHeatmap(mat_cgi, col = c("white", "darkorange"), name = "CGI",
	  top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "darkorange", lty = 1:n_row_group))), 
	  top_annotation_height = unit(2, "cm"), column_title = "CGI", axis_name = axis_name,
	  use_raster = TRUE, raster_quality = 2) 


	ht_list = ht_list + EnrichedHeatmap(mat_mix, name = qq("@{type}CR"), col = cr_col, 
	                top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(neg_col = "darkgreen", pos_col = "red", lty = 1:n_row_group))), 
	                top_annotation_height = unit(2, "cm"), column_title = qq("sig@{type}CR"),
	                use_raster = TRUE, raster_quality = 2, combined_name_fun = NULL, axis_name = axis_name)
	n_heatmap = n_heatmap + 1

	message("append epi heatmaps")
	eval(SNIPPET_APPEND_EPI_HEATMAP)

	message("draw heatmaps")
	eval(SNIPPET_DRAW_HEATMAP)

	return(invisible(NULL))

}

# == title
# Visualizing enrichment for epigenomic signals at gene body
#
# == param
# -cr
# -txdb
# -expr
# -cgi
# -K
# -marks
# -type
# -extend
# -expr_annotation 
#
# == details
#
# == value
# no value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
cr_enriched_heatmap_at_gene = function(cr, txdb, expr, cgi, K = 1, marks = NULL, type = "neg", extend = 5000,
	expr_annotation) {
	
	eval(SNIPPET_ATTACH_CR_PARAM)

	if(is.null(km)) {
		stop("`cr` should be returned by ``.")
	}

	gi = names(km[km == K])
	cr = cr[cr$gene_id %in% gi]
		
	gene = genes(txdb)
	gene = gene[names(gene) %in% cr$gene_id]

	target = gene
	target_name = "gene"
	target_ratio = 0.6
	axis_name = c("-5KB", "TSS", "TES", "5KB")

	expr = expr[names(gene), , drop = FALSE]
	
	km = km[names(gene)]

	# normalize to CGI
	mat_cgi = normalizeToMatrix(cgi, gene, extend = extend, target_ratio = target_ratio, mean_mode = "absolute", w = 50)

	message("normalize to epi signals")
	eval(SNIPPET_NORMALIZE_EPI_SIGNALS)

	gl = width(gene[gene$gene_id])
	gl[gl > quantile(gl, 0.95)] = quantile(gl, 0.95)
	width_anno_name = "gene_width"
	width_anno = gl

	expr = expr[names(gene), sample_id, drop = FALSE]
	if(n_subgroup == 2) {
		expr_mean = rowMeans(expr[, subgroup == subgroup_level[1], drop = FALSE]) - 
					rowMeans(expr[, subgroup == subgroup_level[1], drop = FALSE])
		expr_split = ifelse(expr_mean > 0, "high", "low")
		expr_split = factor(expr_split, levels = c("high", "low"))
	} else {
		expr_split = NULL
	}

	n_upstream_index = length(attr(meth_mat_mean, "upstream_index"))
	meth_split = kmeans(meth_mat_mean[, seq(round(n_upstream_index*0.8), round(n_upstream_index*1.2))], centers = 3)$cluster
	x = tapply(rowMeans(meth_mat_mean[, seq(round(n_upstream_index*0.8), round(n_upstream_index*7/5))]), meth_split, mean)
	od = structure(order(x), names = names(x))
	meth_split = paste0("cluster", od[as.character(meth_split)])

	if(n_subgroup == 2) {
		combined_split = paste(meth_split, expr_split, sep = "|")
		row_order = merge_row_order(meth_mat_mean, list(
			combined_split == "cluster1|high",
			combined_split == "cluster1|low",
			combined_split == "cluster2|high",
			combined_split == "cluster2|low",
			combined_split == "cluster3|high",
			combined_split == "cluster3|low"
		))
	} else {
		combined_split = meth_split
		row_order = merge_row_order(meth_mat_mean, list(
			combined_split == "cluster1",
			combined_split == "cluster2",
			combined_split == "cluster3"
		))
	}

	expr_col_od = do.call("c", lapply(subgroup_level, function(le) {
		dend1 = as.dendrogram(hclust(dist(t(expr[, subgroup == le, drop = FALSE]))))
		hc1 = as.hclust(reorder(dend1, colMeans(expr[, subgroup == le, drop = FALSE])))
		col_od1 = hc1$order
		which(subgroup == le)[col_od1]
	}))

	cor_col_fun = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))
	meth_col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
	cr_col = c("-1" = "darkgreen", "0" = "white", "1" = "red")

	fixed_heatmap = 2
	eval(SNIPPET_HEATMAP_PAGE)

	if(missing(expr_annotation)) {
		if(n_subgroup >= 2) {
			expr_annotation = HeatmapAnnotation(subgroup = subgroup, col = list(subgroup = structure(rand_color(n_subgroup), names = subgroup_level)), 
				show_annotation_name = TRUE, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 10))
		}
	}

	######### construct heatmap list  ############
	## if there are too many heatmaps, they will be put in two pages.
	epi_color = c(brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
	epi_mark_list = list()
	epi_title_list = list()

	n_heatmap = 0
	n_row_group = 3
	ht_list = Heatmap(expr, name = "expr", show_row_names = FALSE,
		show_column_names = FALSE, width = unit(5, "cm"), show_column_dend = FALSE, cluster_columns = FALSE, column_order = expr_col_od,
		top_annotation = expr_annotation, column_title = "Expression", show_row_dend = FALSE,
		use_raster = TRUE, raster_quality = 2)
	n_heatmap = n_heatmap + 1
	
	ht_list = ht_list + rowAnnotation(foo_width = row_anno_points(width_anno, axis = TRUE, gp = gpar(col = "#00000040")),
	                width = unit(1, "cm"))

	ht_list = ht_list + Heatmap(km[names(gene)], name = "km_groups", col = km_col, show_row_names = FALSE,
			width = unit(1, "cm"))

	ht_list = ht_list + EnrichedHeatmap(mat_cgi, col = c("white", "darkorange"), name = "CGI",
	  top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "darkorange", lty = 1:n_row_group))), 
	  top_annotation_height = unit(2, "cm"), column_title = "CGI", axis_name = axis_name,
	  use_raster = TRUE, raster_quality = 2) 
	n_heatmap = n_heatmap + 1

	message("append epi heatmaps")
	eval(SNIPPET_APPEND_EPI_HEATMAP)

	message("draw heatmaps")
	eval(SNIPPET_DRAW_HEATMAP)

	return(invisible(NULL))

}

# == title
# Visualizing enrichment for epigenomic signals at TSS-CGIs
#
# == param
# -cr
# -txdb
# -expr
# -gf
# -fdr_cutoff
# -meth_diff_cutoff
# -marks
# -type
# -extend
# -min_reduce
# -min_width
# -nearest_by
# -expr_annotation 
#
# == details
#
# == value
# no value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
cr_enriched_heatmap_at_genomic_features = function(cr, txdb, expr, gf,
	fdr_cutoff = 0.05, meth_diff_cutoff = 0.1, marks = NULL, type = "neg", extend = 5000,
	min_reduce = 1, min_width = 1000, nearest_by = "tss", expr_annotation) {

	eval(SNIPPET_ATTACH_CR_PARAM)

	message(qq("filter sigCR by fdr_cutoff = @{fdr_cutoff}, meth_diff_cutoff = @{meth_diff_cutoff}"))
	eval(SNIPPET_FILTER_SIG_CR)
	
	gf_origin = gf

	# overlap gf to gene extended regions
	message("extracting gene tss")
	gm = genes(TXDB)
	gm = gm[gm$gene_id %in% unique(cr$gene_id)]
	tss = promoters(gm, upstream = 1, downstream = 0)
	gl = width(gm)
	g = gm
	strand(g) = "*"
	start(g) = start(g) - gm_extend[1]
	end(g) = end(g) + ifelse(length(gm_extend) == 2, gm_extend[2], gm_extend[1])
	start(g) = ifelse(start(g) > 1, start(g), 1)

	mtch = as.matrix(findOverlaps(g, gf))
	gf = gf[unique(mtch[, 2])]
	message(qq("@{length(gf)} are in extended gene regios."))

	if(min_reduce >= 0) {
		gf = reduce(gf, min = min_reduce)
		message(qq("@{length(gf)} remain after merging by min_reduce <= @{min_reduce}"))
	}

	gf = gf[width(gf) >= min_width, ]
	message(qq("@{length(gf)} regions remain after removing regions with width <= @{min_width}"))
	
	# find associated gene (by tss for by gene)
	if(nearest_by == "tss") {
		message("look for nearest tss")
		d = distanceToNearest(gf, tss, select = "all")
		subjectHits = subjectHits(d)
		ind = tapply(seq_len(length(d)), queryHits(d), function(ind) {
			ind[which.max(gl[subjectHits[ind]])][1]
		})
		d = d[as.vector(ind)]
		gf2 = gf[queryHits(d)]
		gf2$distanceToNearest = mcols(d)$distance
		gf2$nearestGene = gm$gene_id[subjectHits(d)]
		gf2$nearestGeneStrand = strand(tss[subjectHits(d)])
		l = gf2$nearestGeneStrand == "+" & start(gf2) < start(tss[subjectHits(d)]) |
		    gf2$nearestGeneStrand == "-" & end(gf2) > end(tss[subjectHits(d)])
		l = as.vector(l)
		gf2$distanceToNearest[l] = -gf2$distanceToNearest[l]
	} else {
		message("look for nearest gene body")
		d = distanceToNearest(gf, gm, select = "all")
		subjectHits = subjectHits(d)
		ind = tapply(seq_len(length(d)), queryHits(d), function(ind) {
			ind[which.max(gl[subjectHits[ind]])][1]
		})
		d = d[as.vector(ind)]
		gf2 = gf[queryHits(d)]
		gf2$distanceToNearest = mcols(d)$distance
		gf2$nearestGene = gm$gene_id[subjectHits(d)]
		gf2$nearestGeneStrand = strand(gm[subjectHits(d)])
		# if the gf is overlapped to gene body, how much of it is overlapped by the gene
		gg = pintersect(gf2, gm[gf2$nearestGene], resolve.empty = "max")
		gf2$overlapGenePercent = width(gg)/width(gf2)
		l = gf2$nearestGeneStrand == "+" & start(gf2) < start(gm[subjectHits(d)]) |
		    gf2$nearestGeneStrand == "-" & end(gf2) > end(gm[subjectHits(d)])
		l = as.vector(l)
		gf2$distanceToNearest[l] = -gf2$distanceToNearest[l]
	}
	message(qq("@{length(gf2)} regions remain after overlapping to genes"))

	strand(gf2) = strand(gm[gf2$nearestGene])
	names(gf2) = gf2$nearestGene
	gf2$gene_id = gf2$nearestGene

	target_ratio = mean(width(gf2))/(sum(extend) + mean(width(gf2)))
	target = gf2
	target_name = "gf"

	message("normalize to sigCR")
	eval(SNIPPET_NORMALIZE_SIG_CR)

	message("normalize to gf")
	mat_gf = normalizeToMatrix(gf_origin, target, extend = extend, w = 50, trim = 0, target_ratio = target_ratio)
	
	if(length(target) < 10) {
		message("too few regions left, just quit.")
		return(invisible(NULL))
	}

	if(is.not.null(km)) {
		km = km[names(target)]
	}

	gf_width = width(gf2)
	gf_width[gf_width > quantile(gf_width, 0.99)] = quantile(gf_width, 0.99)
	width_anno_name = "gf_width"
	width_anno = gf_width

	message("normalize to epi signals")
	eval(SNIPPET_NORMALIZE_EPI_SIGNALS)

	########### prepare the order of rows and columns
	message("determing row and colum orders")
	eval(SNIPPET_ROW_ORDER_AND_COLUMN_ORDER)

	cor_col_fun = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red"))
	meth_col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
	cr_col = c("-1" = "darkgreen", "0" = "white", "1" = "red")
	gf_col = colorRamp2(c(0, 1), c("white", "blue"))

	fixed_heatmap = 3
	eval(SNIPPET_HEATMAP_PAGE)

	if(missing(expr_annotation)) {
		if(n_subgroup >= 2) {
			expr_annotation = HeatmapAnnotation(subgroup = subgroup, col = list(subgroup = structure(rand_color(n_subgroup), names = subgroup_level)), 
				show_annotation_name = TRUE, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 10))
		}
	}

	######### construct heatmap list  ############
	## if there are too many heatmaps, they will be put in two pages.
	epi_color = c(brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
	epi_mark_list = list()
	epi_title_list = list()

	n_heatmap = 0
	n_row_group = 2
	ht_list = Heatmap(expr, name = "expr", show_row_names = FALSE,
		show_column_names = FALSE, width = unit(5, "cm"), show_column_dend = FALSE, cluster_columns = FALSE, column_order = expr_col_od,
		top_annotation = expr_annotation, column_title = "Expression", show_row_dend = FALSE,
		use_raster = TRUE, raster_quality = 2)
	n_heatmap = n_heatmap + 1
	
	if(is.not.null(km)) {
		ht_list = ht_list + Heatmap(km[names(target)], name = "km_groups", col = km_col, show_row_names = FALSE,
			width = unit(0.5, "cm"))
	}

	ht_list = ht_list + EnrichedHeatmap(mat_gf, name = "gf", col = gf_col,
	                top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "blue", lty = 1:n_row_group))), 
	                top_annotation_height = unit(2, "cm"), column_title = qq("gf"),
	                use_raster = TRUE, raster_quality = 2, combined_name_fun = NULL)
	n_heatmap = n_heatmap + 1

	ht_list = ht_list + EnrichedHeatmap(mat_mix, name = qq("@{type}CR"), col = cr_col,
	                top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(neg_col = "darkgreen", pos_col = "red", lty = 1:n_row_group))), 
	                top_annotation_height = unit(2, "cm"), column_title = qq("sig@{type}CR"),
	                use_raster = TRUE, raster_quality = 2, combined_name_fun = NULL)
	n_heatmap = n_heatmap + 1

	ht_list = ht_list + rowAnnotation(foo_width = row_anno_points(width_anno, axis = TRUE, gp = gpar(col = "#00000040")),
	                width = unit(1, "cm"))

	message("append epi heatmaps")
	eval(SNIPPET_APPEND_EPI_HEATMAP)

	message("draw heatmaps")
	eval(SNIPPET_DRAW_HEATMAP)

	return(invisible(NULL))
}


