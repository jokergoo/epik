
# == title
# Compare raw and smoothed methylation
#
# == param
# -gi a single gene id
# -cr_smoothed correlated regions using smoothed methylation
# -txdb transcriptome annotation if ``start`` and ``end`` are not set
# -start start position of the region of interested (in the extended gene region)
# -end end position of the region of interested (in the extended gene region)
#
# == details
# If ``start`` and ``end`` are not set, the whole extended gene will be plotted.
#
# The aim of this function is see whether smoothing can improve the methylation dataset.
#
# There will be six tracks:
#
# - smoothed methylation
# - correlation between methylation and gene expression
# - raw methylation
# - raw methylation for those CpG sites with coverage larger than 25th percential of all CpG coverage
# - raw methylation for those CpG sites with coverage larger than 50th percential of all CpG coverage
# - CpG coverage, the bottom, middle and top lines correspond to 25th, 50th and 75th percential from all samples
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
compare_meth = function(gi, cr_smoothed, cr_raw, txdb = NULL, start = NULL, end = NULL) {

	cr_param = metadata(cr_smoothed)$cr_param
	extend = cr_param$extend
	sample_id = cr_param$sample_id
	subgroup = cr_param$cr_subgroup
	col = cr_param$col
	if(!is.null(subgroup)) col = col[subgroup]

	cr_smoothed = cr_smoothed[cr_smoothed$gene_id == gi]

	if(is.not.null(start) && is.not.null(end)) {
		chr = as.vector(seqnames(cr_smoothed))[1]
		s = start
		e = end
		cr_smoothed = cr_smoothed[start(cr_smoothed) >= s & end(cr_smoothed) <= e]
	} else {
		g = genes(txdb)
		chr = as.vector(seqnames(g))
		s = start(g) - extend
		e = end(g) + extend
	}

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

	if(is.not.null(start) && is.not.null(end)) {
		title = qq("@{chr}:@{start}-@{end}, @{length(site)} CpG sites")
	} else {
		title = qq("@{gi}, extend @{extend} bp to both side, @{length(site)} CpG sites")
	}

	has_cr_raw = !missing(cr_raw)
	if(has_cr_raw) {
		n_track = 7
		track_ylab = c("smoothed meth", "corr", "raw meth", "corr from raw meth", qq("raw meth\n(cov >= @{cov_cutoff_1})(q25)"), qq("raw meth\n(cov >= @{cov_cutoff_2})(q50)"), "median cpg cov")
		track_ylim = c(0, 1, -1, 1, 0, 1, -1, 1, 0, 1, 0, 1, c(0, max(rowMedians(cov))))
	} else {
		n_track = 6
		track_ylab = c("smoothed meth", "corr", "raw meth", qq("raw meth\n(cov >= @{cov_cutoff_1})(q25)"), qq("raw meth\n(cov >= @{cov_cutoff_2})(q50)"), "median cpg cov")
		track_ylim = c(0, 1, -1, 1, 0, 1, 0, 1, 0, 1, c(0, max(rowMedians(cov))))
	}
	gtrellis_layout(GRanges(seqnames = chr, ranges = IRanges(s, e)),
		n_track = n_track, track_ylab = track_ylab, track_ylim = track_ylim,
		title = title)

	add_track(NULL, panel_fun = function(gr) {
		for(i in seq_len(ncol(meth))) {
			grid.lines(site, meth[, i], gp = gpar(col = col[i]), default.units = "native")
		}
	}, use_raster = TRUE, raster_quality = 2)

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
	}, use_raster = TRUE, raster_quality = 2)
	
	# raw methylation
	add_track(NULL, panel_fun = function(gr) {
		for(i in seq_len(ncol(raw))) {
			grid.lines(site, raw[, i], gp = gpar(col = col[i]), default.units = "native")
		}
	}, use_raster = TRUE, raster_quality = 2)

	if(has_cr_raw) {
		add_track(NULL, panel_fun = function(gr) {
			l = cr_raw$corr > 0
			if(sum(l)) {
				grid.rect(mid(ranges(cr_raw[l])), 0, width = width(cr_raw[l]), height = abs(cr_raw$corr[l]), default.units = "native",
					gp = gpar(fill = "red", col = NA), just = "bottom")
			}
			l = cr_raw$corr < 0
			if(sum(l)) {
					grid.rect(mid(ranges(cr_raw[l])), 0, width = width(cr_raw[l]), height = abs(cr_raw$corr[l]), default.units = "native",
					gp = gpar(fill = "green", col = NA), just = "top")
			}
		}, use_raster = TRUE, raster_quality = 2)
	}
	
	# raw methylation with CpG coverage > q25
	add_track(NULL, panel_fun = function(gr) {
		for(i in seq_len(ncol(raw))) {
		    xx = raw[, i]
		    yy = cov[, i]
		    if(sum(yy >= cov_cutoff_1) > 1)
		    	grid.lines(site[yy >= cov_cutoff_1], xx[yy >= cov_cutoff_1], gp = gpar(col = col[i]), default.units = "native")
		}
	}, use_raster = TRUE, raster_quality = 2)

	# raw methylation with CpG coverage > q50
	add_track(NULL, panel_fun = function(gr) {
		for(i in seq_len(ncol(raw))) {
		    xx = raw[, i]
		    yy = cov[, i]
		    if(sum(yy > cov_cutoff_2) > 1)
		    	grid.lines(site[yy >= cov_cutoff_2], xx[yy >= cov_cutoff_2], gp = gpar(col = col[i]), default.units = "native")
		}
	}, use_raster = TRUE, raster_quality = 2)

	cov_quantile = apply(cov, 1, quantile, c(0.25, 0.5, 0.75))
	add_track(NULL, panel_fun = function(gr) {
		grid.segments(site, 0, site, cov_quantile[2, ], default.units = "native")
		# grid.segments(site, cov_quantile[1, ], site, cov_quantile[3, ], default.units = "native", gp = gpar(col = "#FF000040"))
	}, use_raster = TRUE, raster_quality = 2)

}
