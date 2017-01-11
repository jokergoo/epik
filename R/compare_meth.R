
compare_meth = function(g, cr_smoothed) {
	g = sort(g)
	g = g[g$gene_id %in% cr_smoothed$gene_id]

	cr_param = metadata(cr_smoothed)$cr_param
	sample_id = cr_param$sample_id
	subgroup = cr_param$cr_subgroup
	col = cr_param$col
	if(!is.null(subgroup)) col = col[subgroup]
	for(i in seq_len(length(g))) {
		qqcat("compare smoothed and raw methylation for @{g[i]$gene_id}\n")
		compare_meth_in_one_gene(g[i], cr_smoothed[cr_smoothed$gene_id == g[i]$gene_id], sample_id, col)
	}
}

compare_meth_in_one_gene = function(g, cr_smoothed, sample_id, col) {

	cr_param = metadata(cr)$cr_param
	extend = cr_param$extend

	chr = as.vector(seqnames(g))
	s = start(g) - extend
	e = end(g) + extend
	xrange = c(s, e)

	col = add_transparency(col, 0.2)

	methylation_hooks$set_chr(chr, verbose = FALSE)

	meth_gr = methylation_hooks$gr
	meth_site = start(meth_gr)
	ind = extract_sites(s, e, meth_site, TRUE, 0)
	meth = methylation_hooks$meth[ind, sample_id]
	raw = methylation_hooks$raw[ind, sample_id]
	cov = methylation_hooks$cov[ind, sample_id]
	site = meth_site[ind]

	cov_cutoff_1 = round(quantile(cov, 0.25))
	cov_cutoff_2 = round(quantile(cov, 0.5))

	gtrellis_layout(GRanges(seqnames = chr, ranges = IRanges(s, e)),
		n_track = 6, track_ylab = c("smoothed meth", "corr", "raw meth", qq("raw meth\n(cov >= @{cov_cutoff_1})(q25)"), qq("raw meth\n(cov >= @{cov_cutoff_2})(q50)"), "median cpg cov"),
		track_ylim = c(0, 1, -1, 1, 0, 1, 0, 1, 0, 1, c(0, max(rowMedians(cov)))),
		title = qq("@{g$gene_id}, extend @{extend} bp to both side, @{length(site)} CpG sites"))

	add_track(NULL, panel_fun = function(gr) {
		for(i in seq_len(ncol(meth))) {
			grid.lines(site, meth[, i], gp = gpar(col = col[i]), default.units = "native")
		}
	})

	add_track(NULL, panel_fun = function(gr) {
		l = cr_smoothed$corr > 0
		if(sum(l)) {
			grid.rect(mid(ranges(cr_smoothed[l])), 0, width = width(cr_smoothed[l]), height = abs(cr_smoothed$corr[l]), default.units = "native",
				gp = gpar(fill = "red", col = NA), just = "bottom")
		}
		l = cr_smoothed$corr < 0
		if(sum(l)) {
				grid.rect(mid(ranges(cr_smoothed[l])), 0, width = width(cr_smoothed[l]), height = abs(cr_smoothed$corr[l]), default.units = "native",
				gp = gpar(fill = "green", col = NA), just = "top")
		}
	})
	
	# raw methylation
	add_track(NULL, panel_fun = function(gr) {
		for(i in seq_len(ncol(raw))) {
			grid.lines(site, raw[, i], gp = gpar(col = col[i]), default.units = "native")
		}
	})
	
	# raw methylation with CpG coverage > q25
	add_track(NULL, panel_fun = function(gr) {
		for(i in seq_len(ncol(raw))) {
		    xx = raw[, i]
		    yy = cov[, i]
		    if(sum(yy >= cov_cutoff_1) > 1)
		    	grid.lines(site[yy >= cov_cutoff_1], xx[yy >= cov_cutoff_1], gp = gpar(col = col[i]), default.units = "native")
		}
	})

	# raw methylation with CpG coverage > q50
	add_track(NULL, panel_fun = function(gr) {
		for(i in seq_len(ncol(raw))) {
		    xx = raw[, i]
		    yy = cov[, i]
		    if(sum(yy > cov_cutoff_2) > 1)
		    	grid.lines(site[yy >= cov_cutoff_2], xx[yy >= cov_cutoff_2], gp = gpar(col = col[i]), default.units = "native")
		}
	})

	cov_quantile = apply(cov, 1, quantile, c(0.25, 0.5, 0.75))
	add_track(NULL, panel_fun = function(gr) {
		grid.segments(site, 0, site, cov_quantile[2, ], default.units = "native")
		# grid.segments(site, cov_quantile[1, ], site, cov_quantile[3, ], default.units = "native", gp = gpar(col = "#FF000040"))
	})

}
