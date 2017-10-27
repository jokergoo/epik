
SNIPPET_ATTACH_CR_PARAM = expression({
	cr_param = metadata(cr)$cr_param
	sample_id = cr_param$sample_id
	km = cr_param$km
	if(!is.null(km)) {
		km = km[names(km) %in% unique(cr$gene_id)]
		cr = cr[cr$gene_id %in% names(km)]
	}
	km_col = cr_param$group_mean_col
	subgroup = cr_param$subgroup
	subgroup_level = unique(subgroup)
	n_subgroup = length(subgroup_level)
	gm_extend = cr_param$extend
	if(length(gm_extend) == 1) gm_extend = rep(gm_extend, 2)
})

SNIPPET_FILTER_SIG_CR = expression({
	if(!any(grepl("fdr", colnames(mcols(cr))))) {
		stop("use `cr_add_fdr_column()` first")
	}
	if(n_subgroup %in% c(0, 1)) {
		l = cr$corr_fdr < fdr_cutoff & cr$meth_IQR > meth_diff_cutoff
	} else {
		l = cr$corr_fdr < fdr_cutoff & cr$meth_anova_fdr < fdr_cutoff & abs(cr$meth_diameter) > meth_diff_cutoff
	}
	sig_cr = cr[l]

	if(type == "neg") {
		sig_cr_vice = sig_cr[sig_cr$corr > 0]
		sig_cr = sig_cr[sig_cr$corr < 0]
	} else if(type == "pos") {
		sig_cr_vice = sig_cr[sig_cr$corr < 0]
		sig_cr = sig_cr[sig_cr$corr > 0]
	}	
})

SNIPPET_NORMALIZE_SIG_CR = expression({
	mat_cr = normalizeToMatrix(sig_cr, target, extend = extend, mean_mode = "absolute",
		mapping_column = "gene_id", target_ratio = target_ratio)
	if(type == "neg") {
		mat_cr[mat_cr == 1] = -1
	}
	mat_cr_vice = normalizeToMatrix(sig_cr_vice, target, extend = extend, mean_mode = "absolute",
		mapping_column = "gene_id", target_ratio = target_ratio)
	if(type == "pos") {
		mat_cr_vice[mat_cr_vice == 1] = -1
	}
	l = rowSums(abs(mat_cr)) > 0
	mat_cr = mat_cr[l, ]
	target = target[l]
	mat_cr_vice = mat_cr_vice[l, ]

	message(qq("@{length(target)}/@{length(l)} targets remained"))

	mat_mix = mat_cr
	mat_mix[mat_cr == 0] = mat_cr_vice[mat_cr == 0]

})

SNIPPET_NORMALIZE_EPI_SIGNALS = expression({
	epi_mat_list = normalize_epigenomic_signals(cr, target, marks, expr, extend = extend, target_ratio = target_ratio)

	meth_mat_corr = epi_mat_list$meth_mat_corr
	meth_mat_mean = epi_mat_list$meth_mat_mean
	meth_mat_diff = epi_mat_list$meth_mat_diff
	hist_mat_corr_list = epi_mat_list$hist_mat_corr_list
	hist_mat_mean_list = epi_mat_list$hist_mat_mean_list
	hist_mat_diff_list = epi_mat_list$hist_mat_diff_list

})

SNIPPET_ROW_ORDER_AND_COLUMN_ORDER = expression({
	expr = expr[target$gene_id, sample_id, drop = FALSE]
	if(n_subgroup == 2) {
		expr_mean = rowMeans(expr[, subgroup == subgroup_level[1], drop = FALSE]) - 
					rowMeans(expr[, subgroup == subgroup_level[1], drop = FALSE])
		expr_split = ifelse(expr_mean > 0, "high", "low")
		expr_split = factor(expr_split, levels = c("high", "low"))
	} else {
		expr_split = NULL
	}

	target_index = attr(meth_mat_mean, "target_index")
	n_upstream_index = length(attr(meth_mat_mean, "upstream_index"))
	if(target_name == "tss") {
		meth_mat_for_clustering = meth_mat_mean[, seq(round(n_upstream_index*0.8), round(n_upstream_index*7/5))]
		meth_split = kmeans(meth_mat_for_clustering, centers = 2)$cluster
	} else if(target_name == "gene") {
		meth_mat_for_clustering = meth_mat_mean[, seq(round(n_upstream_index*0.8), round(n_upstream_index*7/5))]
		meth_split = kmeans(meth_mat_for_clustering, centers = 3)$cluster
	} else {
		meth_mat_for_clustering = meth_mat_mean[, target_index]
		meth_split = kmeans(meth_mat_for_clustering, centers = 2)$cluster
	}
	x = tapply(rowMeans(meth_mat_for_clustering), meth_split, mean)
	od = structure(order(x), names = names(x))
	meth_split = paste0("cluster", od[as.character(meth_split)])

	if(n_subgroup == 2) {
		combined_split = paste(meth_split, expr_split, sep = "|")
		if(is.null(km)) {
			row_order = merge_row_order(meth_mat_mean, list(
				combined_split == "cluster1|high",
				combined_split == "cluster1|low",
				combined_split == "cluster2|high",
				combined_split == "cluster2|low"
			))
		} else {
			row_order = merge_row_order(meth_mat_mean, list(
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
		}
	} else {
		combined_split = meth_split
		if(is.null(km)) {
			row_order = merge_row_order(meth_mat_mean, list(
				combined_split == "cluster1",
				combined_split == "cluster2"
			))	
		} else {
			row_order = merge_row_order(meth_mat_mean, list(
				combined_split == "cluster1" & km == 1,
				combined_split == "cluster1" & km == 2,
				combined_split == "cluster1" & km == 3,
				combined_split == "cluster1" & km == 4,
				combined_split == "cluster2" & km == 1,
				combined_split == "cluster2" & km == 2,
				combined_split == "cluster2" & km == 3,
				combined_split == "cluster2" & km == 4
			))
		}
	}

	expr_col_od = do.call("c", lapply(subgroup_level, function(le) {
		dend1 = as.dendrogram(hclust(dist(t(expr[, subgroup == le, drop = FALSE]))))
		hc1 = as.hclust(reorder(dend1, colMeans(expr[, subgroup == le, drop = FALSE])))
		col_od1 = hc1$order
		which(subgroup == le)[col_od1]
	}))
})

SNIPPET_HEATMAP_PAGE = expression({
	total_heatmap = fixed_heatmap + is.not.null(meth_mat_corr) + 2 + sum(sapply(hist_mat_corr_list, is.not.null)) + 
	                                                 sum(sapply(hist_mat_mean_list, is.not.null)) +
	                                                 sum(sapply(hist_mat_diff_list, is.not.null))

	n_heatmap_first_page = total_heatmap
	if(total_heatmap > 12) {
		k = fixed_heatmap + (is.not.null(meth_mat_corr)) + 2
		half_cutoff = round(total_heatmap/2)
		for(i in seq_along(marks)) {
			three_n = c(k + (is.not.null(hist_mat_corr_list[[i]])), 
				        k + (is.not.null(hist_mat_corr_list[[i]])) + (is.not.null(hist_mat_mean_list[[i]])),
				        k + (is.not.null(hist_mat_corr_list[[i]])) + (is.not.null(hist_mat_mean_list[[i]])) + (is.not.null(hist_mat_diff_list[[i]])))
			if(!(all(three_n - half_cutoff > 0) || all(three_n - half_cutoff < 0) )) {
				n_heatmap_first_page = k
				break
			}
			k = k + (is.not.null(hist_mat_corr_list[[i]])) + (is.not.null(hist_mat_mean_list[[i]])) + (is.not.null(hist_mat_diff_list[[i]]))
		}
		message(qq("There are two pages of heatmaps: (@{n_heatmap_first_page}, @{total_heatmap - n_heatmap_first_page})"))
	}

})

SNIPPET_APPEND_EPI_HEATMAP = expression({
	if(is.not.null(meth_mat_corr)) {
		ht_list = ht_list + EnrichedHeatmap(meth_mat_corr, col = cor_col_fun, name = qq("corr_meth"), 
		      top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(pos_col = "red", neg_col = "darkgreen", lty = 1:n_row_group))), 
		      top_annotation_height = unit(2, "cm"), column_title = qq("corr_meth"),
		      use_raster = TRUE, raster_quality = 2)
		n_heatmap = n_heatmap + 1
		epi_mark_list[[1]] = "corr_meth"
		epi_title_list[[1]] = "corr_meth"
	}

	# methylation
	ht_list = ht_list + EnrichedHeatmap(meth_mat_mean, col = meth_col_fun, 
		name = "methylation", column_title = qq("meth"),
			heatmap_legend_param = list(title = "methylation"),
			top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "red", lty = 1:n_row_group))),
			use_raster = TRUE, raster_quality = 2)
	n_heatmap = n_heatmap + 1
	epi_mark_list[[1]] = c(epi_mark_list[[1]], "methylation")
	epi_title_list[[1]] = c(epi_title_list[[1]], "meth")

	ht_list = ht_list + EnrichedHeatmap(meth_mat_diff, col = generate_diff_color_fun(meth_mat_diff),
		name = "methylation_diff", column_title = qq("meth_diff"),
			heatmap_legend_param = list(title = "methylation_diff"),
			top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(pos_col = "#df8640", neg_col = "#3794bf", lty = 1:n_row_group))),
			use_raster = TRUE, raster_quality = 2)
	n_heatmap = n_heatmap + 1
	epi_mark_list[[1]] = c(epi_mark_list[[1]], "methylation_diff")
	epi_title_list[[1]] = c(epi_title_list[[1]], "meth_diff")
	
	ht_list2 = NULL
	ht_list1 = ht_list
	# correlation to histone marks
	for(i in seq_along(hist_mat_mean_list)) {
		if(n_heatmap == n_heatmap_first_page) {
			ht_list1 = ht_list
			ht_list = NULL
		}

		epi_mark_list[i+1] = list(NULL)
		epi_title_list[i+1] = list(NULL)
		if(is.not.null(hist_mat_corr_list[[i]])) {
			if(sum(hist_mat_mean_list[[i]]) > 0) {
				ht_list = ht_list + EnrichedHeatmap(hist_mat_corr_list[[i]], col = cor_col_fun, name = qq("corr_@{MARKS[i]}"), 
			          top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(pos_col = "red", neg_col = "darkgreen", lty = 1:n_row_group))), 
		              top_annotation_height = unit(2, "cm"), column_title = qq("corr_@{MARKS[i]}"),
		              use_raster = TRUE, raster_quality = 2)
			    n_heatmap = n_heatmap + 1
			    epi_mark_list[[i+1]] = c(epi_mark_list[[i+1]], qq("corr_@{MARKS[i]}"))
			    epi_title_list[[i+1]] = c(epi_title_list[[i+1]], qq("corr_@{MARKS[i]}"))
			}
	   	}

	   	if(sum(hist_mat_mean_list[[i]]) > 0) {
		    ht_list = ht_list + EnrichedHeatmap(hist_mat_mean_list[[i]], col = generate_color_fun(hist_mat_mean_list[[i]]), 
		    	name = qq("@{MARKS[i]}_mean"), column_title = qq("@{MARKS[i]}"),
				heatmap_legend_param = list(title = qq("@{MARKS[i]}_density")),
				top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "purple", lty = 1:n_row_group))),
				use_raster = TRUE, raster_quality = 2)
			n_heatmap = n_heatmap + 1
			epi_mark_list[[i+1]] = c(epi_mark_list[[i+1]], qq("@{MARKS[i]}_mean"))
			epi_title_list[[i+1]] = c(epi_title_list[[i+1]], qq("@{MARKS[i]}"))
		}

		if(sum(abs(hist_mat_diff_list[[i]])) > 0) {
			ht_list = ht_list + EnrichedHeatmap(hist_mat_diff_list[[i]], name = qq("@{MARKS[i]}_diff"), col = generate_diff_color_fun(hist_mat_diff_list[[i]]),
				column_title = qq("@{MARKS[i]}_diff"),
				heatmap_legend_param = list(title = qq("@{MARKS[i]}_diff")),
				top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(pos_col = "#df8640", neg_col = "#3794bf", lty = 1:n_row_group))),
				use_raster = TRUE, raster_quality = 2)
			n_heatmap = n_heatmap + 1
			epi_mark_list[[i+1]] = c(epi_mark_list[[i+1]], qq("@{MARKS[i]}_diff"))
			epi_title_list[[i+1]] = c(epi_title_list[[i+1]], qq("@{MARKS[i]}_diff"))
		}
	}
	if(n_heatmap > n_heatmap_first_page) {
		ht_list2 = ht_list
	}

})

SNIPPET_DRAW_HEATMAP = expression({
	lines_lgd = Legend(at = c("high", "low"), title = "Lines", legend_gp = gpar(lty = 1:2), type = "lines")

	### generate the heatmaps
	if(is.null(expr_split)) {
		ht_list = ht_list1
	} else {
		ht_list = Heatmap(expr_split, show_row_names = FALSE, name = "expr_split", col = c("high" = "red", "low" = "darkgreen"), width = unit(5, "mm")) + ht_list1
	}
	foo = draw(ht_list,  annotation_legend_list = list(lines_lgd),
			cluster_rows = FALSE, row_order = row_order, show_row_dend = FALSE,
			column_title = qq("cluster by methylation, @{nrow(meth_mat_mean)} rows"), split = meth_split, row_sub_title_side = "left",
			show_heatmap_legend = FALSE, annotation_legend_side = "right")
	add_boxplot_as_column_annotation(foo, width_anno, "foo_width", width_anno_name)

	if(is.not.null(km))	{
		row_order_list = row_order(foo)
		lt = lapply(row_order_list, function(ind) table(km[ind])/length(ind))
		lt = lapply(lt, function(x) x[order(names(x))])
		n = length(row_order_list)
		x_ind = (seq_len(n) - 0.5)/n
		w = 1/n*0.8
		decorate_heatmap_body("km_groups", slice = 1, {
			rg = c(0, 1)
			pushViewport(viewport(y = unit(1, "npc") + unit(1, "mm"), just = "bottom", height = unit(2, "cm"), yscale = rg))
			for(i in seq_along(lt)) {
				y = lt[[i]]
				grid.rect(x_ind[i], cumsum(y), width = w, height = y, just = c("center", "top"), gp = gpar(fill = km_col[names(y)], lty = i))
			}
			upViewport()
		})
	}

	for(i in seq_along(epi_mark_list)) {
		for(j in seq_along(epi_mark_list[[i]])) {
			nm = epi_mark_list[[i]][j]
			if(nm %in% names(ht_list@ht_list)) {
				decorate_title(nm, {
					grid.rect(height = unit(0.8, "npc"), gp = gpar(fill = epi_color[i], col = NA))
					grid.text(epi_title_list[[i]][j], gp = gpar(fontsize = 14))
				}, which = "column")
			}
		}
	}

	if(n_heatmap > n_heatmap_first_page) {
		
		if(interactive()) {
			readline("press Enter to generate heatmaps on the second page:")
		}

		if(is.null(expr_split)) {
			ht_list = ht_list2
		} else {
			ht_list = Heatmap(expr_split, show_row_names = FALSE, name = "expr_split", col = c("high" = "red", "low" = "darkgreen"), width = unit(5, "mm")) + ht_list2
		}
		foo = draw(ht_list, annotation_legend_list = list(lines_lgd),
				cluster_rows = FALSE, row_order = row_order, show_row_dend = FALSE,
				column_title = qq("cluster by methylation, @{nrow(meth_mat_mean)} rows"), split = meth_split, row_sub_title_side = "left",
				show_heatmap_legend = FALSE)
		for(i in seq_along(epi_mark_list)) {
			for(j in seq_along(epi_mark_list[[i]])) {
				nm = epi_mark_list[[i]][j]
				if(nm %in% names(ht_list@ht_list)) {
					decorate_title(nm, {
						grid.rect(height = unit(0.8, "npc"), gp = gpar(fill = epi_color[i], col = NA))
						grid.text(epi_title_list[[i]][j], gp = gpar(fontsize = 14))
					}, which = "column")
				}
			}
		}
	}
})
