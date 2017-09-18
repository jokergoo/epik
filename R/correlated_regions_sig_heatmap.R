

add_mean_methylation = function(gr, sample_id) {
	gr2 = GRanges()
	mean_meth = NULL
	for(chr in unique(as.vector(seqnames(gr)))) {
		
		methylation_hooks$set_chr(chr, verbose = FALSE)

		sub_gr = gr[seqnames(gr) == chr]
		gr_cpg = methylation_hooks$gr
		m = methylation_hooks$meth[, sample_id]

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

mean_overlap = function(gr, gf_list) {
	gr2 = gr
	gr = annotate_to_genomic_features(gr, gf_list, prefix = "__overlap__")
	overlap_mat = as.matrix(mcols(gr)[, grep("^__overlap__", colnames(mcols(gr)))])
	overlap = rowMeans(overlap_mat)
}

mean_chromHMM_overlap = function(gr, cs_list, state) {
	cs_list = lapply(cs_list, function(cs) cs[cs$states == state])
	mean_overlap(gr, cs_list)
}

# == title
# Heatmaps for significant correlated regions
#
# == param
# -sig_cr significant correlated regions, should be processed by `cr_reduce`
# -txdb transcriptome annotation which was used in `correlated_regions`
# -expr expression matrix which was used in `correlated_regions`
# -ha top annotation for the expression heatmap, should be constructed by `ComplexHeatmap::HeatmapAnnotation`
# -gf_list a list of `GenomicRanges::GRanges` objects which are additional genomic features used as annotations
#
# == details
# There are several heatmaps showing associations between difference sources of datasets, each row in the heatmaps is
# a correlated region or other genomic association to this correlated region.
#
# - heatmap for methylation (mean methylation in CR)
# - one column heatmap which shows the methylation difference
# - heatmap for gene expression
# - heatmap describing how genomic features overlap to correlated regions
# - if `chipseq_hooks`$chromHMM is set and there are two subgroups, there is a heamtap showing the difference of
#   the overlapping of different chromatin states in the two groups
# - a point plot showing the distance to nearest TSS
# - overlap to promoter/gene body/intergenic regions
#
# For the list heatmaps, rows are firstly split by negative correlation and positive correlation. In each sub cluster,
# it is split by k-means clustering (four groups), and in each k-means cluster, rows are ordered by hierarchical clustering.
#
# There will also be plots showing general statistics for each annotation.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
sig_cr_heatmap = function(sig_cr, txdb, expr, ha = NULL, gf_list = NULL) {

	cr_param = metadata(sig_cr)$cr_param
	sig_cr$direction = ifelse(sig_cr$corr > 0, "pos", "neg")

	sample_id = cr_param$sample_id
	subgroup = cr_param$subgroup
	subgroup_color = cr_param$col
	subgroup_level = unique(subgroup)
	n_subgroup = length(subgroup_level)

	if(is.null(ha)) ha = HeatmapAnnotation(subgroup = subgroup, col = list(subgroup = subgroup_color))
	
	meth_mat = as.matrix(mcols(sig_cr)[, grepl("mean_meth_", colnames(mcols(sig_cr)))])
	meth_diff = sig_cr$meth_diff
	expr_mat = expr[sig_cr$gene_id, sample_id]
	
	if(!is.null(gf_list)) {
		gr_foo = annotate_to_genomic_features(sig_cr, gf_list)
		overlap_gf = as.matrix(mcols(gr_foo)[, grepl("overlap_to_", colnames(mcols(gr_foo)))])
		colnames(overlap_gf) = gsub("overlap_to_", "", colnames(overlap_gf))
	} else {
		overlap_gf = NULL
	}

	if(!is.null(chipseq_hooks$chromHMM) && n_subgroup == 2) {
		message("overlap to chrome states")
		sample_id_subgroup1 = sample_id[subgroup == subgroup_level[1]]
		sample_id_subgroup2 = sample_id[subgroup == subgroup_level[2]]

		cs_list1 = get_chromHMM_list(sample_id_subgroup1)
		cs_list2 = get_chromHMM_list(sample_id_subgroup2)

		all_states = unique(cs_list1[[1]]$states)
		overlap_diff_cs = matrix(0, nrow = length(sig_cr), ncol = length(all_states))
		colnames(overlap_diff_cs) = all_states
		for(i in seq_along(all_states)) {
			overlap_diff_cs[, i] = mean_chromHMM_overlap(sig_cr, cs_list1, all_states[i]) - mean_chromHMM_overlap(sig_cr, cs_list2, all_states[i])
		}
	} else {
		overlap_diff_cs = NULL
	}

	gm = genes(txdb)
	gl = width(gm)
	names(gl) = names(gm)
	
	## whether it is at tss, gene body or intergenic regions
	ga = ifelse(sig_cr$gene_tss_dist > -1000 & sig_cr$gene_tss_dist < 2000, "tss",
		 	ifelse(sig_cr$gene_tss_dist > 2000 & sig_cr$gene_tss_dist < gl[sig_cr$gene_id], "gene", "intergenic"))

	col_od = do.call("c", lapply(subgroup_level, function(le) {
		dend1 = as.dendrogram(hclust(dist(t(meth_mat[, subgroup == le, drop = FALSE]))))
		hc1 = as.hclust(reorder(dend1, colMeans(meth_mat[, subgroup == le, drop = FALSE])))
		col_od1 = hc1$order
		which(subgroup == le)[col_od1]
	}))

	abs_tss_dist = abs(sig_cr$gene_tss_dist)
	q = quantile(abs_tss_dist, 0.9)
	abs_tss_dist[abs_tss_dist > q] = q

	# rows are split into four slices for neg_cr and pos_cr separately and ordered by mean value
	km_meth1 = kmeans(meth_mat[sig_cr$direction == "neg", ], centers = 4)$cluster
	x = tapply(rowMeans(meth_mat[sig_cr$direction == "neg", ]), km_meth1, mean)
	od = structure(rank(x), names = names(x))
	km_meth1 = od[as.character(km_meth1)]
	km_meth2 = kmeans(meth_mat[sig_cr$direction == "pos", ], centers = 4)$cluster
	x = tapply(rowMeans(meth_mat[sig_cr$direction == "pos", ]), km_meth2, mean)
	od = structure(rank(x), names = names(x))
	km_meth2 = od[as.character(km_meth2)]
	split = numeric(nrow(meth_mat))
	split[sig_cr$direction == "neg"] = paste0("neg", km_meth1)
	split[sig_cr$direction == "pos"] = paste0("pos", km_meth2)

	## now we concatenate heatmaps 

	## 1. a one-column heatmap shows row slices
	ht_list = Heatmap(split, name = "split", show_row_names = FALSE, show_column_names = FALSE, width = unit(5, "mm"),
		col = c(neg1 = "darkgreen", neg2 = "darkgreen", neg3 = "darkgreen", neg4 = "darkgreen", 
			    pos1 = "red", pos2 = "red", pos3 = "red", pos4 = "red"), show_heatmap_legend = FALSE)
	## 2. methylation for the CRs
	meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
	ht_list = ht_list + Heatmap(meth_mat, name = "methylation", col = meth_col_fun,
		show_row_names = FALSE, show_column_names = FALSE, cluster_columns = FALSE, column_order = col_od, 
		top_annotation = ha, column_title = "methylation", show_row_dend = FALSE, combined_name_fun = NULL,
		use_raster = TRUE, raster_quality = 2)

	if(!is.null(meth_diff)) {
		ht_list = ht_list + Heatmap(meth_diff, name = "meth_diff", col = generate_diff_color_fun(meth_diff, col = c("green", "white", "red")),
			show_row_names = FALSE, width = unit(5, "mm"), show_heatmap_legend = FALSE) 
	}
	## 3. expression matrix
	ht_list = ht_list + Heatmap(expr_mat, name = "expr", show_row_names = FALSE,
		show_column_names = FALSE, cluster_columns = FALSE, column_order = col_od,
		top_annotation = ha, column_title = "Expression", show_row_dend = FALSE,
		use_raster = TRUE, raster_quality = 2)
	## 4. overlapping matrix for CGI/shore/tfbs
	ht_list = ht_list + Heatmap(overlap_gf, name = "overlap0", show_row_names = FALSE, col = colorRamp2(c(0, 1), c("white", "orange")),
		show_column_names = TRUE, cluster_columns = FALSE,
		column_title = "overlap to gf", show_row_dend = FALSE)

	if(!is.null(overlap_diff_cs)) {
		ht_list = ht_list + Heatmap(overlap_diff_cs, col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")), 
			name = "overlap_diff", show_row_names = FALSE,
			show_column_names = TRUE, cluster_columns = FALSE,
			column_title = "overlap diff", show_row_dend = FALSE,
			use_raster = TRUE, raster_quality = 2)
	}

	## 6. dist to tss
	rht_list = ht_list + rowAnnotation(tss_dist = row_anno_points(abs_tss_dist, size = unit(1, "mm"), gp = gpar(col = "#00000020"), axis = TRUE), 
		width = unit(2, "cm"), show_annotation_name = TRUE)
	## 7. annotation to genes
	ht_list = ht_list + Heatmap(ga, name = "anno", col = c("tss" = "red", "gene" = "blue", "intergenic" = "green"), show_row_names = FALSE,
		width = unit(5, "mm"))

	draw(ht_list, main_heatmap = "methylation", split = split,
		column_title = qq("@{length(sig_cr)} cr, sum_width=@{sum(width(sig_cr))}"))
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
		if(!is.null(overlap_diff_cs)) {
			decorate_heatmap_body("overlap_diff", slice = i, {
				grid.rect(gp = gpar(col = "black", fill = "transparent"))
			})
		}
	}

	if(!dev.interactive()) {
		
		device_type = .Device
		filepath = attr(.Device, "filepath")

		n_plots = 3
		if(!is.null(overlap_gf)) {
			n_plots = n_plots + ncol(overlap_gf)
		}
		if(!is.null(overlap_diff_cs)) {
			n_plots = n_plots + ncol(overlap_diff_cs)
		}
		nrow_plots = ceiling(sqrt(n_plots))
		ncol_plots = ceiling(n_plots/nrow_plots)
		stat_plot_width = ncol_plots*3
		stat_plot_height = nrow_plots*3
		if(device_type == "pdf") {
			filepath = gsub("\\.pdf$", ".stat.pdf", filepath, ignore.case = TRUE)
			pdf(filepath, width = stat_plot_width, height = stat_plot_height)
		} else if(device_type == "png") {
			filepath = gsub("\\.png$", ".stat.png", filepath, ignore.case = TRUE)
			png(filepath, width = stat_plot_width*200, height = stat_plot_height*200)
		} else {
			filepath = gsub("\\.[^.]+$", ".stat.png", filepath, ignore.case = TRUE)
			png(filepath, width = stat_plot_width*200, height = stat_plot_height*200)
		}
		
		op = par(no.readonly = TRUE)
		on.exit(par(op))

		par(mfrow = c(nrow_plots, ncol_plots))

		plot(1:8, ylim = c(0, 1), type = "n", main = "mean methylation", axes = FALSE, ylab = "mean methylation")
		for(i in seq_along(subgroup_level)) {
			x1 = tapply(rowMeans(meth_mat[, subgroup == subgroup_level[i]]), split, mean)
			lines(1:8, x1, lty = i)
			points(1:8, x1, cex = 1.5, bg = meth_col_fun(x1), pch = 21)
		}
		axis(side = 1, at = 1:8, labels = names(x1))
		axis(side = 2)

		if(!is.null(overlap_gf)) {
			for(i in 1:ncol(overlap_gf)) {
				x = tapply(overlap_gf[, i], split, mean)
				if(max(abs(x)) > 0.05) {
					barplot(x, main = colnames(overlap_gf)[i], col = "orange", ylim = c(0, 1))
				}
			}
		}

		if(!is.null(overlap_diff_cs)) {
			for(i in 1:ncol(overlap_diff_cs)) {
				x = tapply(overlap_diff_cs[, i], split, function(x) {
						c(sum(x[x < 0])/length(x), sum(x[x > 0]/length(x)))
					})
				x = do.call("cbind", x)
				if(max(abs(x)) > 0.05) {
					barplot(abs(x), main = colnames(overlap_diff_cs)[i], col = c("green", "red"), 
						ylim = c(-0.3, 0.3), offset = x[1, ])
				}
			}
		}

		boxplot(split(abs(sig_cr$gene_tss_dist), split), outline = FALSE, main = "dist2tss")
		m = do.call("cbind", tapply(ga, split, table))
		m = apply(m, 2, function(x) x/sum(x))

		barplot(m, main = "annotation to genes", col = c("tss" = "red", "gene" = "blue", "intergenic" = "green")[rownames(m)])

		dev.off()
	}

	return(invisible(NULL))
}
