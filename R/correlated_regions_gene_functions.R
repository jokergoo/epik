
do_david = function(gene_list, email) {
	
	res_list = lapply(gene_list, function(genes) {
		genes = gsub("\\.\\d+", "", genes)
		submit_to_david(genes, email)
	})
}

if(!is.memoised(do_david)) {
	do_david = memoise(do_david)
}

# Functional enrichment for CR genes by DAVID
#
# == param
# -cr correlated regions returned by `cr_enriched_heatmap`
# -email email for DAVID API (https://david.ncifcrf.gov/content.jsp?file=WS.html )
# -count_cutoff minimal number of CR genes in a function term
# -fdr_cutoff cutoff of fdr of the enrichment test
# -pop_count_cutoff maximum number of population genes in a function term.
#
# == details
# Genes in k-means group 1 and 4 are sent to DAVID web server to do functional enrichment. 
# The significant functions are visualized as a heatmap.
#
# Only three Gene Ontology (biological process, molecular function and cellular component) categories are used.
#
# There is also a heatmap which shows the significant enrichment.
#
# == value
# A list of function enrichments.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
cr_genes_function_enrichment = function(cr, email, hit_cutoff = 5, fdr_cutoff = 0.05,
	km_group = 1:4) {
	
	cr_param = metadata(cr)$cr_param
	subgroup = cr_param$subgroup
	n_subgroup = length(unique(subgroup))

	km4 = cr_param$km
	if(is.null(km4)) {
		stop("`cr` should be returned by `cr_enrichedheatmap()`.")
	}
	group_mean_col = cr_param$group_mean_col
	row_order = cr_param$order
	combined_split = cr_param$combined_split
	expr_split = cr_param$expr_split

	l = km4 %in% km_group
	km4 = km4[l]
	row_order = row_order[l]
	combined_split = combined_split[l]
	if(!is.null(expr_split)) {
		expr_split = expr_split[l]
	}

	cr_gi = names(km4)

	cr = cr[cr$gene_id %in% cr_gi]

	gene_list = split(names(km4), cr_param$combined_split)

	df_all_list = lapply(gene_list, function(genes) {
		Sys.sleep(5)
		submit_to_david(genes, email)
	})
	df_list = lapply(df_all_list, function(x) {
		reduce_david_results(x, fdr_cutoff = fdr_cutoff, hit_cutoff = hit_cutoff)
	})
	for(i in names(gene_list)) {
		df_list[[i]] = cbind(df_list[[i]], cluster = as.numeric(i))
	}
	term_df = do.call("rbind", df_list)

	unique_terms = unique(term_df$termName)
	gene_mat = matrix(0, nrow = length(km4), ncol = length(unique_terms))
	rownames(gene_mat) = names(km4)
	colnames(gene_mat) = unique_terms
	for(i in seq_len(nrow(term_df))) {
		g_vector = term_df[, "geneIds"][[i]]
		gene_mat[g_vector, term_df[i, "termName"]] = 1
	}

	ht_list = Heatmap(km4, name = "cluster", col = group_mean_col, show_row_names = FALSE, width = 2*grobHeight(textGrob("1", gp = gpar(fontsize = 10))))
	if(!is.null(expr_split)) {
		ht_list = ht_list + Heatmap(expr_split, name = "expr_direction", show_row_names = FALSE, width = unit(5, "mm"), 
			col = c("high" = "Red", "low" = "darkgreen"))
	}
	unique_category_name = unique(term_df$categoryName)
	category_col = structure(brewer.pal(length(unique_category_name), "Set1"), names = unique_category_name)
	ht_list = ht_list + Heatmap(gene_mat, name = "identity", col = c("0" = "white", "1" = "blue"), show_row_names = FALSE, 
			cluster_columns = FALSE, combined_name_fun = NULL, column_names_max_height = max_text_width(colnames(gene_mat), gp = gpar(fontsize = 12)),
			use_raster = TRUE, raster_quality = 2, show_row_dend = FALSE, split = combined_split, show_heatmap_legend = FALSE,
			top_annotation = HeatmapAnnotation(categoryName = term_df$categoryName, col = list(categoryName = category_col), show_annotation_name = TRUE))
	draw(ht_list, main_heatmap = "identity", row_order = row_order, cluster_rows = FALSE)

	
	for(i in seq_along(df_list)) {
		for(j in seq_len(ncol(gene_mat))) {
			current_term_list = df_list[[i]]$Term
			ind = which(gsub("GO:.*~", "", df_list[[i]]$Term) == term_names[j])
			if(length(ind)) {
				message(qq("@{term_names[j]} is significant in @{names(df_list)[i]}"))
				decorate_heatmap_body("identity", slice = i, {
					grid.rect(x = j/n_term, width = 1/n_term, default.units = "native", just = "right", gp = gpar(fill = "transparent", lwd = 2))
					grid.text(paste(df_list[[i]][ind, "listHits"], sprintf("%.1e", df_list[[i]][ind, "benjamini"]), sep = ", "), 
						x = (j-0.5)/n_term, y = unit(1, "npc") - unit(2, "mm"), rot = 90, just = "right", gp = gpar(fontsize = 8))
				})
			}
		}
	}
	draw(ht_list, main_heatmap = "GO")

	return(invisible(res_list))
}
