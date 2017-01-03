
generate_diff_color_fun = function(x) {
	q = quantile(x, c(0.05, 0.95), na.rm = TRUE)
	max_q = max(abs(q))
	colorRamp2(c(-max_q, 0, max_q), c("#3794bf", "#FFFFFF", "#df8640"))
}


hilbert_curve_methylation_difference = function(subgroup, comparison, chromosome = paste0("chr", 1:22), 
	species = "hg19", type = c("global_mean", "subgroup_mean", "difference")) {

	qqcat("split genome by 1kb window\n")
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
		if(length(sample_id) == 0) {
			if(length(subgroup) != length(sample_id)) {
				stop("length of `subgroup` should be same as the number of samples in `methylation_hooks`.")
			}
		}
		subgroup = subgroup[sample_id]
	}

	le = unique(subgroup)

	## difference of methylation
	# mean methylaiton by 1kb window, compared to histone modifications which are also segmented by 1kb window
	qqcat("calculating mean methylation in every 1kb window\n")
	gr_meth = get_mean_methylation_in_genomic_features(sample_id, chromosome = chromosome, 
		gf_list = list(chromGr_1kb_window = chromGr_1kb_window), filter_fun = function(s) TRUE)[[1]]
	mat = mcols(gr_meth)
	mat = as.matrix(mat[, -ncol(mat)])
	mcols(gr_meth) = NULL

	for(i in seq_along(le)) {
		qqcat("calculating mean methylation in subgroup @{le[i]}\n")
		mcols(gr_meth)[, qq("mean_@{le[i]}")] = rowMeans(mat[, subgroup == le[i], drop = FALSE], na.rm = TRUE)
	}
	gr_meth$mean = rowMeans(mat, na.rm = TRUE)

	if("global_mean" %in% type) {
		qqcat("making hilbert curve for the global mean methylation\n")
		col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
		cm = ColorMapping(col_fun = col_fun)
		lgd = color_mapping_legend(cm, title = "mean_meth", plot = FALSE)

		hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean methylation in all samples"), legend = lgd)
		hc_layer(hc, gr_meth, col = col_fun(gr_meth$mean))
		hc_map(hc, add = TRUE, fill = NA, border = "#808080")
		seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
		hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
		hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
	}

	if("subgroup_mean" %in% type) {
		for(i in seq_along(le)) {
			qqcat("making hilbert curve for the mean methylation in subgroup @{le[i]}\n")
			col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
			cm = ColorMapping(col_fun = col_fun)
			lgd = color_mapping_legend(cm, title = "mean_meth", plot = FALSE)

			hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean methylation in @{le[i]}"), legend = lgd)
			hc_layer(hc, gr_meth, col = col_fun(mcols(gr_meth)[, qq("mean_@{le[i]}")]))
			hc_map(hc, add = TRUE, fill = NA, border = "#808080")
			seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
			hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
			hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
		}
	}

	if("difference" %in% type) {
		qqcat("making hilbert curve for the methylation difference between @{comparison[1]} and @{comparison[2]}\n")
		diff = mcols(gr_meth)[, qq("mean_@{comparison[1]}")] - mcols(gr_meth)[, qq("mean_@{comparison[2]}")]
		col_fun = generate_diff_color_fun(diff)
		cm = ColorMapping(col_fun = col_fun)
		lgd = color_mapping_legend(cm, title = "meth_diff", plot = FALSE)

		hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = qq("mean methylation difference"), legend = lgd)
		hc_layer(hc, gr_meth, col = col_fun(diff))
		hc_map(hc, add = TRUE, fill = NA, border = "#808080")
		seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
		hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
		hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
	}

	return(invisible(gr_meth))
}


average_in_window = function(gr1, gr2, v, empty_value = 0, chromosome = unique(as.vector(seqnames(gr1)))) {
	x = rep(empty_value, length(gr1))
	chromosome = intersect(chromosome, intersect( unique(as.vector(seqnames(gr1))),  unique(as.vector(seqnames(gr2)))))
	for(chr in chromosome) {
		l = as.vector(seqnames(gr1) == chr)
		if(sum(l)) {
			window = ranges(gr1[seqnames(gr1) == chr])
			ir = ranges(gr2[seqnames(gr2) == chr])

			mtch = as.matrix(findOverlaps(window, ir))
			v2 = HilbertCurve:::average_in_window(window, ir, mtch, v, "w0", empty_value)
			x[as.vector(seqnames(gr1) == chr)][unique(mtch[, 1])] = v2
		}
	}
	return(x)
}



hilbert_curve_chipseq_difference = function(mark, subgroup, comparison, chromosome = paste0("chr", 1:22), 
	species = "hg19", type = c("global_mean", "subgroup_mean", "abs_difference", "rel_difference"),
	density_column = "density") {

	qqcat("split genome by 1kb window\n")
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

	hm_list = get_peak_list(mark, sample_id)

	qqcat("calculating mean @{mark} signal in every 1kb window\n")
	gr = chromGr_1kb_window
	mat = do.call("cbind", lapply(hm_list, function(x) average_in_window(chromGr_1kb_window, x, mcols(x)[, density_column], chromosome = chromosome)))
	colnames(mat) = sample_id

	mean_density = sapply(hm_list, function(x) sum(as.numeric(width(x))*mcols(x)[, density_column])/sum(as.numeric(width(x))))
	
	for(i in seq_along(le)) {
		qqcat("calculating mean @{mark} signal in subgroup @{le[i]}\n")
		mcols(gr)[, qq("mean_@{le[i]}")] = rowMeans(mat[, subgroup == le[i], drop = FALSE], na.rm = TRUE)
	}
	gr$mean = rowMeans(mat, na.rm = TRUE)

	if("global_mean" %in% type) {
		qqcat("making hilbert curve for the global mean signal for @{mark}\n")
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
		col_fun = generate_diff_color_fun(mcols(gr))
		for(i in seq_along(le)) {
			qqcat("making hilbert curve for the mean signal for @{mark} in subgroup @{le[i]}\n")
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
		qqcat("making hilbert curve for the absolute signal difference between @{comparison[1]} and @{comparison[2]}\n")
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
		qqcat("making hilbert curve for the relative signal difference between @{comparison[1]} and @{comparison[2]}\n")
		diff = mcols(gr)[, qq("mean_@{comparison[1]}")] - mcols(gr)[, qq("mean_@{comparison[2]}")]
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


