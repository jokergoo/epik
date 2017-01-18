
cr_hilbert = function(cr, chromosome = paste0("chr", 1:22), merge_chr = TRUE) {

	chr_len = read.chromInfo()$chr.len

	## all cr windows
	col_fun = colorRamp2(c(-1, 0, 1), c("#4DAF4A", "white", "#E41A1C"))
	cm = ColorMapping(col_fun = col_fun)
	lgd = color_mapping_legend(cm, title = "type", plot = FALSE)
	if(merge_chr) {
		hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = "cr for all chromosomes", legend = lgd)
	    hc_layer(hc, cr, col = col_fun(cr$corr), mean_mode = "absolute")
	    hc_map(hc, add = TRUE, fill = NA, border = "#808080")

		seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
		hc = GenomicHilbertCurve(chr = chromosome, mode = "normal", level = 6, newpage = FALSE)
		hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))

	} else {
		for(i in seq_along(chromosome)) {
		    chr = chromosome[i]
		    cat(chr, "\n")
		    cr2 = cr[seqnames(cr) == chr]
		    hc = HilbertCurve(s = 1, e = max(chr_len), mode = "pixel", level = 10, title = chr, legend = lgd)
		    hc_layer(hc, ranges(cr2), col = col_fun(cr2$corr), mean_mode = "absolute")
		}
	}
}
