
# == title
# Visualize methylation by Hilbert curve
#
# == param
# -subgroup subgroup information which corresponds to sample IDs stored in `methylation_hooks`.
#      The value can be a vector with same length as sample IDs or a named vector that names are sample IDs.
# -comparison if there are more than one subgroups, the comparison of two subgroups which shows the methylation difference
#        The value is a vector of length of two and the difference is calculated as subgroup[1] - subgroup[2]
# -chromosome a vector of chromosome
# -species species
# -type Three types of visualization supported, see "details" section
# 
# == details
# Genome is segmented by 1kb window and mean methylation for each 1kb window is calculated, later visualized
# by Hilbert curve.
#
# There are three types of visualization methods:
#
# -global_mean the mean methylation averaged from all samples
# -subgroup_mean the mean methylation averaged in every subgroup
# -difference the difference of methylation in two subgroups
#
# == value
# a `GenomicRanges::GRanges` object which contains mean methylation for the 1kb segmentation and other statistics.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
hilbert_curve_methylation_difference = function(subgroup, comparison, chromosome = paste0("chr", 1:22), 
	species = "hg19", type = c("global_mean", "subgroup_mean", "difference")) {

	message("split genome by 1kb window")
	chromInfo = getChromInfoFromUCSC(species)
	chromInfo = chromInfo[chromInfo$chrom %in% chromosome, ]
	chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
	chromGr_1kb_window = makeWindows(chromGr, w = 1000, short.keep = TRUE)
	mcols(chromGr_1kb_window) = NULL

	methylation_hooks$set_chr(sample(chromosome, 1), verbose = FALSE)
	sample_id = methylation_hooks$sample_id

	if(is.null(names(subgroup))) {
		if(length(subgroup) != length(sample_id)) {
			stop("length of `subgroup` should be same as the number of samples in `methylation_hooks`.")
		}
	} else {
		sample_id = intersect(names(subgroup), sample_id)
		subgroup = subgroup[sample_id]
	}

	le = unique(subgroup)

	## difference of methylation
	# mean methylaiton by 1kb window, compared to histone modifications which are also segmented by 1kb window
	message("calculating mean methylation in every 1kb window")
	gr_meth = get_mean_methylation_in_genomic_features(sample_id, chromosome = chromosome, 
		genomic_features = list(chromGr_1kb_window = chromGr_1kb_window))[[1]]
	mat = mcols(gr_meth)
	mat = as.matrix(mat[, -ncol(mat)])
	mcols(gr_meth) = NULL

	l = apply(mat, 1, function(x) sum(is.na(x))/length(x) < 0.1)
	# mat = mat[l, ]
	# gr_meth = gr_meth[l]

	for(i in seq_along(le)) {
		message(qq("calculating mean methylation in subgroup @{le[i]}"))
		mcols(gr_meth)[, qq("mean_@{le[i]}")] = rowMeans(mat[, subgroup == le[i], drop = FALSE], na.rm = TRUE)
	}
	gr_meth$mean = rowMeans(mat, na.rm = TRUE)

	if("global_mean" %in% type) {
		message("making hilbert curve for the global mean methylation")
		col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
		cm = ColorMapping(col_fun = col_fun)
		lgd = color_mapping_legend(cm, title = "mean_meth", plot = FALSE)

		hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean methylation in all samples"), legend = lgd)
		hc_layer(hc, gr_meth[l], col = col_fun(gr_meth[l]$mean))
		hc_map(hc, add = TRUE, fill = NA, border = "#808080")
		seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
		hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
		hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
	}

	if("subgroup_mean" %in% type) {
		for(i in seq_along(le)) {
			message(qq("making hilbert curve for the mean methylation in subgroup @{le[i]}"))
			col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
			cm = ColorMapping(col_fun = col_fun)
			lgd = color_mapping_legend(cm, title = "mean_meth", plot = FALSE)

			hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean methylation in @{le[i]}"), legend = lgd)
			hc_layer(hc, gr_meth[l], col = col_fun(mcols(gr_meth)[l, qq("mean_@{le[i]}")]))
			hc_map(hc, add = TRUE, fill = NA, border = "#808080")
			seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
			hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
			hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
		}
	}

	if("difference" %in% type) {
		if(missing(comparison)) comparison = le[1:2]
		message(qq("making hilbert curve for the methylation difference between @{comparison[1]} and @{comparison[2]}"))
		diff = mcols(gr_meth)[l, qq("mean_@{comparison[1]}")] - mcols(gr_meth)[l, qq("mean_@{comparison[2]}")]
		col_fun = generate_diff_color_fun(diff)
		cm = ColorMapping(col_fun = col_fun)
		lgd = color_mapping_legend(cm, title = "meth_diff", plot = FALSE)

		hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean methylation difference"), legend = lgd)
		hc_layer(hc, gr_meth[l], col = col_fun(diff))
		hc_map(hc, add = TRUE, fill = NA, border = "#808080")
		seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
		hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
		hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
	}

	return(invisible(gr_meth))
}


average_in_window = function(gr1, gr2, v, empty_value = 0, chromosome = unique(as.vector(seqnames(gr1))), window_size = NULL) {
	x = rep(empty_value, length(gr1))
	chromosome = intersect(chromosome, intersect( unique(as.vector(seqnames(gr1))),  unique(as.vector(seqnames(gr2)))))
	for(chr in chromosome) {
		message(qq("  averaging in window for @{chr}"))
		l = as.vector(seqnames(gr1) == chr)
		l2 = as.vector(seqnames(gr2) == chr)
		if(sum(l)) {
			window = ranges(gr1[l])
			ir = ranges(gr2[l2])
			vx = v[l2]

			mtch = as.matrix(findOverlaps(window, ir))
			v2 = HilbertCurve:::average_in_window(window, ir, mtch, vx, "w0", empty_value)
			x[l][unique(mtch[, 1])] = v2
			
		}
	}
	return(x)
}


# == title
# Visualize ChIP-Seq signals by Hilbert curve
#
# == param
# -mark name of the histone mark, should also be supported in `chipseq_hooks`
# -subgroup subgroup information which corresponds to sample IDs stored in `methylation_hooks`.
#      The value should be a named vector that names are sample IDs.
# -comparison if there are more than one subgroups, the comparison of two subgroups which shows the methylation difference
#        The value is a vector of length of two and the difference is calculated as subgroup[1] - subgroup[2]
# -chromosome a vector fo chromosome
# -species species
# -type four types of visualization supported, see "details" section
# -density_column the column name of the signal defined in `chipseq_hook`$peak
# 
# == details
# Genome is segmented by 1kb window and mean signal for each 1kb window is calculated, later visualized
# by Hilbert curve.
#
# There are four types of visualization methods:
#
# -global_mean the mean signal averaged from all samples
# -subgroup_mean the mean signal averaged in every subgroup
# -abs_difference the absolute difference of signal in two subgroups
# -rel_difference the relative difference in two subgroups. The value is calculated as absolute difference divided by mean singal.
#
# == value
# a `GenomicRanges::GRanges` object which contains mean signal for the 1kb segmentation and other statistics.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
hilbert_curve_chipseq_difference = function(mark, subgroup, comparison, chromosome = paste0("chr", 1:22), 
	species = "hg19", type = c("global_mean", "subgroup_mean", "abs_difference", "rel_difference"),
	density_column = "density") {

	message("split genome by 1kb window")
	chromInfo = getChromInfoFromUCSC(species)
	chromInfo = chromInfo[chromInfo$chrom %in% chromosome, ]
	chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
	chromGr_1kb_window = makeWindows(chromGr, w = 1000, short.keep = TRUE)
	mcols(chromGr_1kb_window) = NULL

	sample_id = chipseq_hooks$sample_id(mark)

	if(is.null(names(subgroup))) {
		stop("`subgroup` must have names to correspond to sample IDs.")
	} else {
		sample_id = intersect(names(subgroup), sample_id)
		subgroup = subgroup[sample_id]
	}

	le = unique(subgroup)

	hm_list = get_peak_list(mark, sample_id, chr = chromosome)

	message(qq("calculating mean @{mark} signal in every 1kb window"))
	gr = chromGr_1kb_window
	mat = do.call("cbind", lapply(hm_list, function(x) average_in_window(chromGr_1kb_window, x, mcols(x)[, density_column], chromosome = chromosome)))
	colnames(mat) = sample_id

	mean_density = sapply(hm_list, function(x) sum(as.numeric(width(x))*mcols(x)[, density_column])/sum(as.numeric(width(x))))
	
	for(i in seq_along(le)) {
		message(qq("calculating mean @{mark} signal in subgroup @{le[i]}"))
		mcols(gr)[, qq("mean_@{le[i]}")] = rowMeans(mat[, subgroup == le[i], drop = FALSE], na.rm = TRUE)
	}
	gr$mean = rowMeans(mat, na.rm = TRUE)

	if("global_mean" %in% type) {
		message(qq("making hilbert curve for the global mean signal for @{mark}"))
		col_fun = generate_diff_color_fun(gr$mean)
		cm = ColorMapping(col_fun = col_fun)
		lgd = color_mapping_legend(cm, title = "mean_signal", plot = FALSE)

		hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean @{mark} signal in all samples"), legend = lgd)
		hc_layer(hc, gr, col = col_fun(gr$mean))
		hc_map(hc, add = TRUE, fill = NA, border = "#808080")
		seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
		hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
		hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
	}

	if("subgroup_mean" %in% type) {
		col_fun = generate_diff_color_fun(mcols(gr)[, qq("mean_@{le[i]}")])
		for(i in seq_along(le)) {
			message(qq("making hilbert curve for the mean signal for @{mark} in subgroup @{le[i]}"))
			cm = ColorMapping(col_fun = col_fun)
			lgd = color_mapping_legend(cm, title = "mean_signal", plot = FALSE)

			hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean @{mark} signal in @{le[i]}"), legend = lgd)
			hc_layer(hc, gr, col = col_fun(mcols(gr)[, qq("mean_@{le[i]}")]))
			hc_map(hc, add = TRUE, fill = NA, border = "#808080")
			seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
			hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
			hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
		}
	}

	if("abs_difference" %in% type) {
		if(missing(comparison)) comparison = le[1:2]
		message(qq("making hilbert curve for the absolute signal difference between @{comparison[1]} and @{comparison[2]}"))
		diff = mcols(gr)[, qq("mean_@{comparison[1]}")] - mcols(gr)[, qq("mean_@{comparison[2]}")]
		col_fun = generate_diff_color_fun(diff)
		cm = ColorMapping(col_fun = col_fun)
		lgd = color_mapping_legend(cm, title = "abs_diff", plot = FALSE)

		hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean absolute @{mark} signal difference"), legend = lgd)
		hc_layer(hc, gr, col = col_fun(diff))
		hc_map(hc, add = TRUE, fill = NA, border = "#808080")
		seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
		hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
		hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
	}

	if("rel_difference" %in% type) {
		if(missing(comparison)) comparison = le[1:2]
		message(qq("making hilbert curve for the relative signal difference between @{comparison[1]} and @{comparison[2]}"))
		diff = mcols(gr)[, qq("mean_@{comparison[1]}")] - mcols(gr)[, qq("mean_@{comparison[2]}")]
		mean = (mcols(gr)[, qq("mean_@{comparison[1]}")] + mcols(gr)[, qq("mean_@{comparison[2]}")])/2
		s0 = quantile(mean[mean > 1e-6], 0.05)
		# diff = diff/(mean + s0)
		diff = diff/mean(mean_density)
		col_fun = generate_diff_color_fun(diff)
		cm = ColorMapping(col_fun = col_fun)
		lgd = color_mapping_legend(cm, title = "rel_diff", plot = FALSE)

		hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean relative @{mark} signal difference"), legend = lgd)
		hc_layer(hc, gr, col = col_fun(diff))
		hc_map(hc, add = TRUE, fill = NA, border = "#808080")
		seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
		hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
		hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
	}

	return(invisible(gr))
}

get_chipseq_association_stat = function(gr_list, q_cutoff) {

	n_gr = length(gr_list)
	gr_name = names(gr_list)
	gr_name_combine = outer(as.vector(outer(gr_name, gr_name, paste, sep = "_vs_")), 
		                    as.vector(outer(c("up", "down"), c("up", "down"), paste, sep = "_")),
		                    paste, sep = "_")
	gr_name_combine = as.vector(gr_name_combine)
	jaccard_coef = matrix(nrow = length(gr_name_combine), ncol = length(q_cutoff), dimnames = list(gr_name_combine, q_cutoff))
	jaccard_intersect = jaccard_coef
	jaccard_union = jaccard_coef
	jaccard_gr1 = jaccard_coef
	jaccard_gr2 = jaccard_coef

	gr_list = lapply(gr_list, function(gr) {
		gr[abs(gr$diff) > 1e-6]
	})

	for(k in seq_along(q_cutoff)) {
		for(i in 1:(n_gr-1)) {
			for(j in (i+1):n_gr) {
				gr1 = gr_list[[i]]
				gr2 = gr_list[[j]]
				
				nm = qq("@{gr_name[i]}_vs_@{gr_name[j]}")
				x1 = abs(gr1$diff)
				c1 = quantile(x1, q_cutoff[k])
				x2 = abs(gr2$diff)
				c2 = quantile(x2, q_cutoff[k])
					
				for(d1 in c("up", "down")) {
					for(d2 in c("up", "down")) {
						nm2 = qq("@{nm}_@{d1}_@{d2}")
						if(d1 == "up") {
							gr_subset1 = gr1[gr1$diff > c1]
						} else {
							gr_subset1 = gr1[gr1$diff < -c1]
						}
						if(d2 == "up") {
							gr_subset2 = gr2[gr2$diff > c2]
						} else {
							gr_subset2 = gr2[gr2$diff < -c2]
						}
						jaccard_coef[nm2, k] = genomic_corr_jaccard(gr_subset1, gr_subset2)
						jaccard_intersect[nm2, k] = sum(as.numeric(width(intersect(gr_subset1, gr_subset2))))
						jaccard_union[nm2, k] = sum(as.numeric(width(union(gr_subset1, gr_subset2))))
						jaccard_gr1[nm2, k] = sum(as.numeric(width(gr_subset1)))
						jaccard_gr2[nm2, k] = sum(as.numeric(width(gr_subset2)))
					}
				}
				
				message(qq("associating @{nm}, at q_@{q_cutoff[k]}"))
			}
		}
	}

	l = is.na(jaccard_coef[, 1])
	jaccard_coef = jaccard_coef[!l, ]
	jaccard_intersect = jaccard_intersect[!l, ]
	jaccard_union = jaccard_union[!l, ]
	jaccard_gr1 = jaccard_gr1[!l, ]
	jaccard_gr2 = jaccard_gr2[!l, ]

	lt = list(jaccard_coef = jaccard_coef,
		      jaccard_intersect = jaccard_intersect,
		      jaccard_union = jaccard_union,
		      jaccard_gr1 = jaccard_gr1,
		      jaccard_gr2 = jaccard_gr2)
	return(lt)
}

if(!is.memoized(get_chipseq_association_stat)) {
	get_chipseq_association_stat = memoise(get_chipseq_association_stat)
}

# == title
# General association between histome marks
#
# == param
# -gr_list a list of `GenomicRanges::GRanges` objects which show signal difference in two groups.
#          Each `GenomicRanges::GRanges` object must have a ``diff`` column. The object is usually from `hilbert_curve_chipseq_difference`.
# -q quantile of difference
#
# == details
# For each pair of histome marks, the Jaccard coefficient for the regions which show higher difference
# than ``q`` is calcualted and visualized.
#
# == value
# no value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
general_chipseq_association = function(gr_list, q = 0.9) {
	
	q_cutoff = seq(0.1, 0.9, by = 0.1)
	
	if(!all(q %in% q_cutoff)) {
		stop("q should be in `seq(0.1, 0.9, by = 0.1)`.")
	}

	lt = get_chipseq_association_stat(gr_list, q_cutoff)
	jaccard_coef = lt$jaccard_coef
	jaccard_union = lt$jaccard_union
	jaccard_intersect = lt$jaccard_intersect
	jaccard_union = lt$jaccard_union
	jaccard_gr1 = lt$jaccard_gr1
	jaccard_gr2 = lt$jaccard_gr2

	single_plot = function(q = 0.9, l) {

		k = which(colnames(jaccard_coef) == as.character(q))

		max_y = max(jaccard_coef[l, ], na.rm = TRUE)*1.1
		max_u = max(jaccard_union[, as.character(q)])
		direction1 = sapply(strsplit(rownames(jaccard_coef[l, ]), "_"), "[", 4)
		direction2 = sapply(strsplit(rownames(jaccard_coef[l, ]), "_"), "[", 5)

		txt = gsub("_(down|up)_(down|up)$", "", rownames(jaccard_coef[l, ]))
		txt_width = max_text_width(txt, gp = gpar(fontsize = 8))*1.1
		txt_width = max(unit.c(min(unit.c(txt_width, unit(0.2, "npc"))), txt_width))
		plot_region_width = unit(1, "npc") - unit(2, "cm") - unit(2, "cm") - txt_width - unit(5, "mm")
		pushViewport(viewport(x = unit(2, "cm"), y = unit(2, "cm"), 
			width = plot_region_width,
			height = unit(1, "npc") - unit(5, "mm") - unit(2, "cm"),
			just = c("left", "bottom"),
			xscale = c(0.09, 0.91), yscale = c(0, max_y)))
		grid.rect()
		for(i in which(l)) {
			grid.lines(q_cutoff, jaccard_coef[i, ], default.units = "native")
		}
		grid.lines(c(q, q), c(0, max_y), default.units = "native", gp = gpar(col = "grey", lty = 2))
		grid.text(qq("q_@{q}"), unit(q, "native") - unit(2, "mm"), unit(max_y, "native") - unit(2, "mm"), default.units = "native", 
			gp = gpar(fontsize = 8, col = "grey"), just = c("right", "top"))
		grid.points(rep(q, sum(l)), jaccard_coef[l, k], default.units = "native", pch = 16, size = unit(2, "mm"))
		grid.xaxis(gp = gpar(fontsize = 8))
		grid.yaxis(gp = gpar(fontsize = 8))

		grid.text("Quantile of absolute difference", x = unit(0.5, "native"), y = unit(-15, "mm"))
		grid.text("Jaccard coefficient", x = unit(-15, "mm"), y = unit(0.5, "npc"), rot = 90)
		
		od = order(jaccard_coef[l, k])
		
		x = jaccard_coef[l, k][od]
		labels = rownames(jaccard_coef)[l][od]
		txt = txt[od]
		text_height = convertHeight(grobHeight(textGrob(txt, gp = gpar(fontsize = 8)))*2, "native", valueOnly = TRUE)
		h1 = x - text_height*0.75
		h2 = x + text_height*0.75
		min_y = convertHeight(unit(-1, "cm"), "native", valueOnly = TRUE)
		pos = smartAlign(h1, h2, c(min_y, max_y))
		h = (pos[, 1] + pos[, 2])/2

		w = txt_width
		w2 = jaccard_union[labels, k]/max_u*w

		s1 = unit(q, "native")
		s2 = unit(0.9, "native") + unit(2, "cm")
		grid.segments(s1, x, (s2 - s1)*0.333 + s1, x, default.units = "native", gp = gpar(col = "#00000040"))
		grid.segments((s2 - s1)*0.333 + s1, x, (s2 - s1)*0.667 + s1, h, default.units = "native", gp = gpar(col = "#00000040"))
		grid.segments((s2 - s1)*0.667 + s1, h, s2, h, default.units = "native", gp = gpar(col = "#00000040"))
		grid.text(txt, x = s2 + w*0.5, y = h, default.units = "native", gp = gpar(fontsize = 8), just = "center")
		
		grid.rect(s2, h - 0.5*text_height, width = w2, height = unit(1, "mm"), default.units = "native", just = c("left", "top"),
			gp = gpar(fill = "#00000060", col = "#00000060"))
		grid.rect(s2, h, 
			width = w*(jaccard_gr1[labels, k]/jaccard_union[labels, k]),
			height = text_height, default.units = "native", just = c("left"),
			gp = gpar(fill = ifelse(direction1 == "up", "#FF000040", "#00FF0040")))
		grid.rect(s2 + w, h, 
			width = w*(jaccard_gr2[labels, k]/jaccard_union[labels, k]),
			height = text_height, default.units = "native", just = c("right"),
			gp = gpar(fill = ifelse(direction2 == "up", "#FF000040", "#00FF0040")))
		
		major.by = round(max_u/4)
	    digits = as.numeric(gsub("^.*e([+-]\\d+)$", "\\1", sprintf("%e", major.by)))
	    major.by = round(major.by, digits = -1 * digits)
		major.at = seq(0, max_u, by = major.by)
		if(major.by > 1e+06) {
			major.tick.labels = paste(major.at/1e+06, "MB", sep = "")
		}
		else if(major.by > 1000) {
			major.tick.labels = paste(major.at/1000, "KB", sep = "")
		}
		else {
			major.tick.labels = paste(major.at, "bp", sep = "")
		}

		grid.segments(s2, unit(-1.2, "cm"), s2 + w, unit(-1.2, "cm"))
		breaks_pos = s2 + w*(major.at/max_u)
		grid.segments(breaks_pos, unit(-1.2, "cm"), breaks_pos, unit(-1.3, "cm"))
		grid.text(major.tick.labels, breaks_pos, unit(-1.4, "cm"), just = "top", gp = gpar(fontsize = 8))
		upViewport()

	}

	for(j in seq_along(q)) {
		message(qq("venn diagrams are marked for q_@{q[j]}"))
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2)))
		pushViewport(viewport(layout.pos.row = 1:2, layout.pos.col = 1))
		single_plot(q[j], l = grepl("down_up|up_down", rownames(jaccard_coef)))
		upViewport()
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
		single_plot(q[j], l = grepl("up_up", rownames(jaccard_coef)))
		upViewport()
		pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
		single_plot(q[j], l = grepl("down_down", rownames(jaccard_coef)))
		upViewport()
		upViewport()
	}
}

# == title
# General association between histone modifications and methylations
#
# == param
# -gr_list a list of `GenomicRanges::GRanges` objects which show signal difference in two groups.
#          Each `GenomicRanges::GRanges` object must have a ``diff`` column. The object is usually from `hilbert_curve_chipseq_difference`.
# -gr_meth a `GenomicRanges::GRanges` object which shows methylation difference in two groups.
#          The object is usually from `hilbert_curve_methylation_difference`
#
# == details
# For each histone mark, the distribution of methylation difference in regions which show
# high histome modification signal difference is visualized.
#
# == value
# no value is returned
# 
# == author
# Zuguang Gu <z.gu@dkfz.de> 
general_chipseq_association_to_methylation = function(gr_list, gr_meth) {
	
	gr_name = names(gr_list)
	p_neg = vector("list", length(gr_list))
	names(p_neg) = names(gr_list)
	m_neg_box = p_neg
	m_pos_box = p_neg

	if(is.null(gr_meth$diff)) {
		stop("`gr_meth` should have a 'diff' column which contains methylation difference between two subgroups.")
	}

	q_cutoff = seq(0.1, 0.9, by = 0.1)

	for(i in seq_along(gr_list)) {

		gr = gr_list[[i]]

		if(is.null(gr$diff)) {
			stop(qq("`gr_list$@{gr_name[i]}` should have a 'diff' column which contains difference between two subgroups."))
		}
		gr = gr[abs(gr$diff) > 1e-6]
		for(q in q_cutoff) {
			message(qq("associate @{gr_name[i]} to methylation with |signal| > q_@{q}"))
			c = quantile(abs(gr$diff), q)

			gr_neg = gr[gr$diff < -c]
			gr_pos = gr[gr$diff > c]
			mtch1 = as.matrix(findOverlaps(gr_neg, gr_meth))
			x1 = gr_meth[unique(mtch1[, 2])]$diff
			mtch2 = as.matrix(findOverlaps(gr_pos, gr_meth))
			x2 = gr_meth[unique(mtch2[, 2])]$diff
			p_neg[[i]][as.character(q)] = sum(as.numeric(width(gr_neg)))/(sum(as.numeric(width(gr_neg))) + sum(as.numeric(width(gr_pos))))
			m_neg_box[[i]][[as.character(q)]] = quantile(x1, c(0.25, 0.5, 0.75), na.rm = TRUE)
			m_pos_box[[i]][[as.character(q)]] = quantile(x2, c(0.25, 0.5, 0.75), na.rm = TRUE)
		}
		m_neg_box[[i]] = do.call("cbind", m_neg_box[[i]])
		m_pos_box[[i]] = do.call("cbind", m_pos_box[[i]])
	}

	## make the plot
	message("making plot")
	single_plot = function(i, title = "") {
		pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 1, 
			height = unit.c(3*grobHeight(textGrob(title)), unit(c(3, 1), "null"), unit(5, "mm") + 2*grobHeight(textGrob("foo", gp = gpar(fontsize = 8)))))))
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
		grid.text(title)
		upViewport()

		max = max(abs(c(m_neg_box[[i]], m_pos_box[[i]])))*1.05
		pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1, xscale = c(0.05, 0.95), yscale = c(-max, max)))
		grid.rect()
		l = p_neg[[i]] > 0.1
		grid.lines(q_cutoff[l], m_neg_box[[i]][2, l], default.units = "native", gp = gpar(col = "green", lwd = 2))
		grid.polygon(c(q_cutoff[l], rev(q_cutoff[l])), c(m_neg_box[[i]][1, l], rev(m_neg_box[[i]][3, l])), 
			default.units = "native", gp = gpar(fill = "#00FF0060", col = NA))
		l = 1 - p_neg[[i]] > 0.1
		grid.lines(q_cutoff[l], m_pos_box[[i]][2, l], default.units = "native", gp = gpar(col = "red", lwd = 2))
		grid.polygon(c(q_cutoff[l], rev(q_cutoff[l])), c(m_pos_box[[i]][1, l], rev(m_pos_box[[i]][3, l])), 
			default.units = "native", gp = gpar(fill = "#FF000060", col = NA))
		grid.yaxis(gp = gpar(fontsize = 8))
		if(ci == 1) {
			grid.text("methylation difference", x = unit(-1.4, "cm"), rot = 90, gp = gpar(fontsize = 10))
		}
		upViewport()

		pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1, xscale = c(0.05, 0.95), yscale = c(0, 1.1)))
		grid.rect(q_cutoff, p_neg[[i]], width = 0.1, height = p_neg[[i]][l], default.units = "native",
			just = c("top"), gp = gpar(fill = "#B0FFB0"))
		grid.rect(q_cutoff, p_neg[[i]], width = 0.1, height = 1 - p_neg[[i]][l], default.units = "native",
			just = c("bottom"), gp = gpar(fill = "#FFB0B0"))
		grid.xaxis(gp = gpar(fontsize = 8))
		grid.yaxis(at = c(0, 0.5, 1), gp = gpar(fontsize = 8))
		if(ci == 1) {
			grid.text("propotion", x = unit(-1.4, "cm"), rot = 90, gp = gpar(fontsize = 10))
		}
		upViewport()

		upViewport()
	}
	n_gr = length(gr_list)
	nr = round(sqrt(n_gr))
	nc = ceiling(n_gr/nr)
	pushViewport(viewport(x = unit(0.5, "cm"), width = unit(1, "npc") - unit(0.5, "cm") - unit(2, "mm"), just = "left"))
	pushViewport(viewport(layout = grid.layout(nrow = nr, ncol = nc)))
	for(i in seq_along(gr_list)) {
		ri = ceiling(i / nc)
		ci = (i - 1) %% nc + 1
		pushViewport(viewport(layout.pos.row = ri, layout.pos.col = ci))
		pushViewport(viewport(x = unit(1.2, "cm"), width = unit(1, "npc") - unit(1.2, "cm"), just = "left"))
		single_plot(i, title = gr_name[i])
		upViewport()
		upViewport()
	}
	upViewport()
	upViewport()
}