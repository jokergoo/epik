### functions related to QC
## the functions make related plots and output some tables

### mat_density_plot for coverage and methylation, add mean and median lines on the heatmap


# == title
# Basic qc plot for distribution of methylation and CpG coverage
#
# == param
# -sample_id a vector of sample ids
# -chromosome a vector of chromosome names
# -background background regions
#
# == detail
# For each sample id, it will produce five plots:
#
# - mean/median CpG coverage per chromosome
# - histogram of CpG coverage
# - mean/median methylation per chromosome 
# - histogram of methylation
# - mean/median Methylation at each CpG coverage 
#
# == value
# A list of corresponding statistics
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
wgbs_qcplot = function(sample_id, chromosome = paste0("chr", 1:22), background = NULL) {

	# coverage and methylation per chromosome
	data = rep(list(list(cov = NULL, meth = NULL, strand = NULL, cov_count = NULL)), length(sample_id))
	names(data) = sample_id
	
	for(chr in chromosome) {
		
		methylation_hooks$set_chr(chr, verbose = FALSE)
		if(!is.null(background)) {
			mtch = as.matrix(findOverlaps(methylation_hooks$gr, background))
		}
		for(sid in sample_id) {
			qqcat("split data for @{sid}, @{chr}")

			cv = methylation_hooks$cov[, sid]
			ind = which(cv != 0)
			data[[sid]]$cov_count[[chr]] = c("zero" = length(cv) - length(ind), "all" = length(cv))
			if(!is.null(background)) {
				ind = intersect(ind, mtch[, 1])
			}
			if(length(ind) > 10000) {
				cat(", randomly sample 10000 CpG sites.")
				ind = sort(sample(ind, 10000))
			}
			cat("\n")
			cv = cv[ind]

			mh = as.vector(methylation_hooks$meth[ind, sid])
			
			strd = as.vector(strand(methylation_hooks$gr[ind]))
			
			data[[sid]]$cov[[chr]] = cv
			data[[sid]]$meth[[chr]] = mh
			data[[sid]]$strand[[chr]] = strd
		}
	}

	for(sid in sample_id) {

		qqcat("making qc plot for @{sid}\n")

		cov = data[[sid]]$cov
		meth = data[[sid]]$meth
		strand = data[[sid]]$strand
		cov_count = data[[sid]]$cov_count

		par(mfrow = c(2, 3))
		
		# mean coverage per chromosome
		cpg_coverage_mean = sapply(cov, mean)
		cpg_coverage_median = sapply(cov, median)
		plot(c(0, length(cpg_coverage_mean)), c(0, max(c(cpg_coverage_mean, cpg_coverage_median))), axes = FALSE, ann = FALSE, type="n")
		for(i in seq_along(cpg_coverage_mean)) {
			abline(v = i, lty = 3, col = "grey")
			lines(c(i-1, i), c(cpg_coverage_mean[i], cpg_coverage_mean[i]), lwd = 2)
			lines(c(i-1, i), c(cpg_coverage_median[i], cpg_coverage_median[i]), lwd = 2, col = "red")
		}
		abline(v = 0, lty = 3, col = "grey")
		par(las = 3)
		axis(side = 1, at = seq_along(cpg_coverage_mean)-0.5, labels = names(cpg_coverage_mean))
		axis(side = 2)
		box()
		par(las = 0)
		title(main = qq("Coverage per chromosome (@{sid})"), ylab = "mean and median CpG coverage")
		legend("bottomleft", lty=1, col = c("black", "red"), legend = c("mean", "median"))

		
		# coverage distribution
		cov_count = do.call("cbind", cov_count)
		zero_cov_rate = sprintf("%.2f", sum(cov_count[1, ])/sum(cov_count[2, ])*100)
		if(all(unique(unlist(strand)) %in% c("+", "-"))) {
			x = unlist(cov)
			q99 = quantile(unlist(cov), 0.99)
			y = unlist(strand)
			x1 = x[y == "+"]
			x2 = x[y == "-"]
			ta = table(x)
			ta1 = table(x1)
			ta2 = table(x2)
			xlim = range(c(as.numeric(names(ta)), as.numeric(names(ta1)), as.numeric(names(ta2))))
			ylim = range(ta, ta1, ta2)
			plot(as.numeric(names(ta)), ta, xlim = xlim, ylim = ylim, main = qq("histogram of CpG coverage (@{sid})\n@{zero_cov_rate}% have zero coverage"), log = "x", axes = FALSE, type = "h", ylab = "", xlab="CpG coverage")
			axis(side = 1)
			breaks = seq(0, max(ta)/sum(ta), by = 0.02)
			axis(side = 2, at = breaks*sum(ta), labels = breaks)
			box()
			par(new = TRUE)
			plot(as.numeric(names(ta1))+0.2, ta1, xlim = xlim, ylim = ylim, type = "h", col = "red", log = "x", axes = FALSE, ann = FALSE)
			par(new = TRUE)
			plot(as.numeric(names(ta2))+0.4, ta2, xlim = xlim, ylim = ylim, type = "h", col = "blue", log = "x", axes = FALSE, ann = FALSE)
			#axis(side = 2)
			legend("topright", lty = 1, col = c("black", "red", "blue"), legend = c("strand *", "strand +", "strand -"))
			par(new = FALSE)
		} else {
			ta = table(unlist(cov))
			q99 = quantile(unlist(cov), 0.99)
			plot(as.numeric(names(ta)), ta, main = qq("histogram of CpG coverage (@{sid})\n@{zero_cov_rate}% have zero coverage"), log = "x", axes = FALSE, type = "h", ylab = "", xlab="CpG coverage")
			abline(v = q99, lty = 2, col = "blue"); text(q99, 0, "q99", adj = c(0, 0))
			axis(side = 1)
			breaks = seq(0, max(ta)/sum(ta), by = 0.02)
			axis(side = 2, at = breaks*sum(ta), labels = breaks)
			box()
		}

		# ## barplot of zero-coverage and non-zero-coverage
		# plot(c(0, length(cov_count)), c(0, max(unlist(cov_count))), axes = FALSE, ann = FALSE, type = "n")
		# for(i in seq_along(cov_count)) {
		# 	rect(i-1, 0, i, cov_count[[i]]["zero"], col = "orange")
		# 	rect(i-1, cov_count[[i]]["zero"], i, cov_count[[i]]["all"], col = "blue")
		# 	abline(v = i, lty = 3, col = "grey")
		# }
		# abline(v = 0, lty = 3, col = "grey")
		# par(las = 3)
		# axis(side = 1, at = seq_along(cov_count) - 0.5, labels = names(cov_count))
		# axis(side = 2)
		# box()
		# par(las = 0)
		# title(main = qq("%zeor/non_zero CpG coverage (@{sid})"), ylab = "number of CpG sites")
		# legend("topright", pch = 15, col = c("orange", "blue"), legend = c("zero", "non-zero"))

		# mean methylation per chromosome
		cpg_methyrate_mean = sapply(meth, mean)
		cpg_methyrate_median = sapply(meth, median)
		plot(c(0, length(cpg_methyrate_mean)), c(0, 1), axes = FALSE, ann = FALSE, type = "n")
		for(i in seq_along(cpg_methyrate_mean)) {
			abline(v = i, lty = 3, col = "grey")
			lines(c(i-1, i), c(cpg_methyrate_mean[i], cpg_methyrate_mean[i]), lwd = 2)
			lines(c(i-1, i), c(cpg_methyrate_median[i], cpg_methyrate_median[i]), lwd = 2, col = "red")
		}
		abline(v = 0, lty = 3, col = "grey")
		par(las = 3)
		axis(side = 1, at = seq_along(cpg_methyrate_mean) - 0.5, labels = names(cpg_methyrate_mean))
		axis(side = 2)
		box()
		par(las = 0)
		title(main = qq("methylation per chromosome (@{sid})"), ylab = "mean and median methylation")
		legend("bottomleft", lty=1, col = c("black", "red"), legend = c("mean", "median"))

		# distribution of methylation on all chromosomes
		hist(unlist(meth), main = qq("histogram of methylation (@{sid})"), xlab = "methylation")

		
		# methylation to coverage
		coverage2methyrate = tapply(unlist(meth), unlist(cov), mean)
		plot(as.numeric(names(coverage2methyrate)), coverage2methyrate, ylim = c(0, 1), pch=16, col = "#000000A0", log = "x", cex = 0.8, xlab = "CpG coverage", ylab = "mean methylation", main = qq("Mean Methylation for each CpG coverage (@{sid})"))
		coverage2methyrate = tapply(unlist(meth), unlist(cov), median)
		points(as.numeric(names(coverage2methyrate)), coverage2methyrate, pch=16, cex = 0.8, col = "#FF0000A0")
		legend("bottomleft", pch = 16, col = c("black", "red"), legend = c("mean", "median"))
		abline(v = q99, lty = 2, col = "blue"); text(q99, 0, "q99 of cov", adj = c(0, 0))
		par(mfrow = c(1, 1))
	}
	
	return(invisible(data))
}

# == title
# Plot coverage and methylation for a single sample
#
# == param
# -sid a single sample id
# -chromosome a vector of chromosome names
# -species species
# -nw number of windows
# -pct point type
# -pt_gp graphic parameters for points (``col`` will be excluded)
# -transparency transparency of points
# -title title of the plot
# -... pass to `gtrellis::gtrellis_layout`
#
# == details
# The whole genome is segented by ``nw`` windows and mean methylation and mean CpG coverage
# are visualized as two tracks.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
gtrellis_coverage_and_methylation = function(sid, chromosome = paste0("chr", 1:22), 
	species = "hg19", nw = 10000, pch = ".", pt_gp = gpar(), transparency = 0.8, 
	title = qq("Distribution of CpG coverage and methylation for @{sid}"), ...) {

	# window size
	chr_len = read.chromInfo(species = species)$chr.len[chromosome]
	w = round(sum(chr_len)/nw)
	qqcat("number of windows: @{nw}, window size: @{w} bp\n")

	flag = 0
	col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"), transparency = transparency)

	for(chr in chromosome) {
		methylation_hooks$set_chr(chr, verbose = FALSE)

		meth = methylation_hooks$meth[, sid]
		cov = methylation_hooks$cov[, sid]; cov = log10(cov+1)
		gr = methylation_hooks$gr
		chr_gr = GRanges(seqnames = chr, ranges = IRanges(1, chr_len[chr]))
		chr_window = makeWindows(chr_gr, w = w)
		mtch = as.matrix(findOverlaps(chr_window, gr))

		gr2 = chr_window[unique(mtch[,1])]
		meth = tapply(meth[ mtch[,2] ], mtch[,1], mean, na.rm = TRUE)
		cov = tapply(cov[ mtch[,2] ], mtch[,1], mean, na.rm = TRUE)
		
		if(flag == 0) {
			gtrellis_layout(category = chromosome, species = species, 
				n_track = 2, track_ylim = c(0, quantile(cov, 0.99), 0, 1), 
				track_ylab = c("log10(coverage)", "methylation"),
				add_name_track = TRUE, add_ideogram_track = TRUE, ...)
			flag = 1
		}
		qqcat("making plot for @{chr}\n")
		add_points_track(gr2, cov, track = 2, category = chr, pch = pch, 
			gp = gp_c(gpar(col = add_transparency("black", transparency)), pt_gp))
		add_points_track(gr2, meth, track = 3, category = chr, pch = pch, 
			gp = gp_c(gpar(col = col_fun(meth)), pt_gp))
	}
}

gp_c = function(gp1, gp2) {
	gp = c(gp1, gp2)
	class(gp) = "gpar"
	gp
}

add_transparency = function(col, transparency = 0) {
	rgb(t(col2rgb(col)/255), alpha = 1 - transparency)
}


# == title
# Plot methylation for multiple samples as heatmaps
#
# == param
# -sample_id a vector of sample ids
# -subtype annotation of samples (e.g. subtypes)
# -chromosome a vector of chromosome names
# -species species
# -nw number of windows
# -title title of the plot
# -... pass to `gtrellis::gtrellis_layout`
#
# == details
# The whole genome is segented by ``nw`` windows. Methylation in different subgroups are visualized as separated tracks.
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
gtrellis_methylation_for_multiple_samples = function(sample_id, subgroup, 
	chromosome = paste0("chr", 1:22), species = "hg19", nw = 1000, 
	title = qq("genome-wide methylation for @{length(sample_id)} samples"), ...) {
	
	col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
	type = unique(subgroup)
		
	n = table(subgroup)[type]
	ty = numeric(2*length(n))
	ty[seq_along(n)*2-1] = 0.5
	ty[seq_along(n)*2] = n + 0.5

	gtrellis_layout(category = chromosome, species = species, track_axis = FALSE, title = title,
		n_track = length(type), track_ylab = type, track_ylim = ty, track_height = n,
		add_name_track = TRUE, add_ideogram_track = TRUE, ...)

	# window size
	chr_len = read.chromInfo(species = species)$chr.len[chromosome]
	w = round(sum(chr_len)/nw)
	qqcat("number of windows: @{nw}, window size: @{w} bp\n")

	for(chr in chromosome) {
		methylation_hooks$set_chr(chr, verbose = FALSE)

		meth = methylation_hooks$meth[, sample_id]
		gr = methylation_hooks$gr
		chr_gr = GRanges(seqnames = chr, ranges = IRanges(1, chr_len[chr]))
		chr_window = makeWindows(chr_gr, w = w)
		mtch = as.matrix(findOverlaps(chr_window, gr))

		gr2 = chr_window[unique(mtch[,1])]
		meth = tapply(mtch[,2], mtch[,1], function(i) colMeans(meth[i, , drop = FALSE]))

		meth = do.call("rbind", meth)

		for(i in seq_along(type)) {
			qqcat("making plot for @{chr}, @{type[i]}\n")
			# sid = sample_id[subgroup == type[i]]
			# m = meth[, sid, drop = FALSE]
			# add_track(gr2, track = i+1, category = chr, panel.fun = function(gr) {
			# 	x = (start(gr2) + end(gr2))/2
			# 	for(i in seq_along(sid)) {
			# 		y = rep(i, length(x)) + (runif(length(x))-0.5)*0.8
			# 		grid.points(x, y, pch = ".", gp = gpar(col = col_fun(m[, i])))
			# 	}
			# })
			m = meth[, subgroup == type[i], drop = FALSE]
			add_heatmap_track(gr2, m, category = chr, track = i + 1, fill = col_fun)
			add_track(gr2, category = chr, track = i + 1, panel_fun = function(gr) {
				grid.rect(gp = gpar(fill = "transparent"))
			})
		}
	}
}


# == title 
# Visualize distribution of a matrix or a list
#
# == param
# -x a matrix or a list. If it is a matrix, distribution in columns are visualized
# -subgroup subgroup information
# -reorder_column if it is true, samples are first ordered by subgroups and in each subgroup, samples are
#        ordered by median values
# -ha additional annotation can be specified as a `ComplexHeatmap::HeatmapAnnotation` object
# -type three types of plots are supported, see details
# -title title for the plot
# -... pass to `ComplexHeatmap::densityHeatmap`
#
# == details
# Three types of plots for visualizing distributions are supported:
#
# -densityHeatmap: distribution is visualized as heatmaps, use `ComplexHeatmap::densityHeatmap`
# -lineplot: distribution is visualized as normal line plot, use `graphics::matplot`
# -MDS: multiple dimension scaling
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
mat_dist = function(x, subgroup = NULL, reorder_column = TRUE, od = if(is.matrix(x)) seq_len(ncol(x)) else seq_along(x), 
	ha = if(is.null(subgroup)) NULL else HeatmapAnnotation(subgroup = subgroup, show_annotation_name = TRUE), 
	type = c("densityHeatmap", "MDS"), title = title, ...) {

	if("densityHeatmap" %in% type) {
		if(reorder_column) {
			if(inherits(x, "list")) {
				od = order(factor(subgroup, levels = unique(subgroup), ordered = TRUE), sapply(x, median, na.rm = TRUE))
			} else {
				od = order(factor(subgroup, levels = unique(subgroup), ordered = TRUE), apply(x, 2, median, na.rm = TRUE))
			}
		} 

		densityHeatmap(x, anno = ha, title = title, column_order = od, ...)
	}

	if("lineplot" %in% type) {
		# use line plot to represent distributions
		par(xpd = FALSE)
		if(is.matrix(x)) {
			den_x = matrix(nrow = 512, ncol = dim(x)[2])
			den_y = matrix(nrow = 512, ncol = dim(x)[2])
			for(i in seq_len(dim(x)[2])) {
				den = density(x[, i], na.rm = TRUE)
				den_x[, i] = den$x
				den_y[, i] = den$y
			}
		} else {
			den_x = matrix(nrow = 512, ncol = length(x))
			den_y = matrix(nrow = 512, ncol = length(x))
			for(i in seq_along(x)) {
				den = density(x[[i]], na.rm = TRUE)
				den_x[, i] = den$x
				den_y[, i] = den$y
			}
		}
		matplot(den_x, den_y, type = "l", col = col[anno], xlab = "value", ylab = "density", main = qq("density distribution: @{title}"))
	}

	if("MDS" %in% type) {
		## MDS plot
		if(is.data.frame(x) || is.matrix(x)) {
			mat = as.matrix(x)
			loc = cmdscale(dist2(t(mat), pairwise_fun = function(x, y) {l = is.na(x) | is.na(y); x = x[!l]; y = y[!l]; sqrt(sum((x-y)^2))}))
			
			if(is.null(ha)) {
				col = structure(seq_along(unique(subgroup)), names = unique(subgroup))
				col_v = col[subgroup]
			} else {
				col = ha@anno_list$subgroup@color_mapping@colors
				col_v = col[subgroup]
			}
			plot(loc[, 1], loc[, 2], pch = 16, cex = 1, col = col_v, main = qq("MDS:@{title}"), xlab = "dimension 1", ylab = "dimension 2")
			legend("bottomleft", pch = 16, legend = names(col), col = col)

			# plot(loc[, 1], loc[, 2], type = "n", pch = 16, cex = 1, col = col[anno], main = qq("MDS:@{title}"), xlab = "dimension 1", ylab = "dimension 2")
			text(loc[, 1], loc[, 2], colnames(x), col = col_v, cex = 0.8)
		}
	}
	return(od)
}


# == title
# Global methylation distribution
# 
# == param
# -sample_id a vector of sample ids
# -subgroup subgroup information
# -reorder_column if it is true, samples are first ordered by subgroups and in each subgroup, samples are
#        ordered by median values
# -ha additional annotation can be specified as a `ComplexHeatmap::HeatmapAnnotation` object
# -chromosome chromosomes
# -by_chr whether make the plot by chromosome
# -background background to look into. The value can be a single `GenomicRanges::GRanges` object or a list of `GenomicRanges::GRanges` objects.
# -p probability to randomly sample CpG sites
# -meth_range the range of methylation on the plot
#
# == details
# The distribution density is visualized as heatmaps.
#
# == value
# If ``by_chr`` is set to ``FALSE``, it returns a vector of column order.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
global_methylation_distribution = function(sample_id, subgroup, 
	reorder_column = TRUE, 
	ha = HeatmapAnnotation(subgroup = subgroup, show_annotation_name = TRUE), 
	chromosome = paste0("chr", 1:22), by_chr = FALSE, 
	# max_cov = 100,
	background = NULL, p = NULL, meth_range = c(0, 1)) {
	
	###############################################
	# distribution of global methylation
	if(inherits(background, "list")) {
		meth_list = vector("list", length(sample_id))
		# cov_list = NULL

		if(length(background) != length(sample_id)) {
			stop("Since you specified `background` as a list, the length should be same as `sample_id`.")
		}

		for(chr in chromosome) {

			methylation_hooks$set_chr(chr, verbose = FALSE)

			meth_gr = methylation_hooks$gr
			ind_list = lapply(seq_along(sample_id), function(i) {
				mtch = as.matrix(findOverlaps(meth_gr, background[[i]]))
				ind = unique(mtch[, 1])
				nr = length(ind)
				if(is.null(p) || by_chr) p <<- min(c(3000, nr))/nr
				ind = ind[sample(c(FALSE, TRUE), nr, replace = TRUE, prob = c(1-p, p))]

				qqcat("random sampled @{length(ind)} sites from @{nr} sites on @{chr} in @{sample_id[i]} (with p = @{sprintf('%.1e', p)})\n")
				ind
			})

			current_meth_list = lapply(seq_along(ind_list), function(i) {
				m = methylation_hooks$meth[ind_list[[i]], sample_id[i]]
				if(!is.null(methylation_hooks$cov)) {
					cov = methylation_hooks$cov[ind_list[[i]], sample_id[i]]
					l = cov == 0
					l[is.na(l)] = TRUE
					m[l] = NA
				}
				m
			})
			# current_cov_list = lapply(seq_along(ind_list), function(i) {
			# 	cov = methylation_hooks$cov[ind_list[[i]], sample_id[i]]
			# 	cov[cov == 0] = NA
			# 	cov[cov > max_cov] = NA
			# 	log10(cov)
			# })
			qqcat("on average there are @{round(mean(sapply(current_meth_list, function(x) sum(is.na(x)))))} CpG without coverage information.\n")
			
			if(by_chr) {
				# it can be for some chromosomes, no CpG sites are sampled
				try(od <- mat_dist(current_meth_list, reorder_column = reorder_column, subgroup = subgroup, ha = ha, title = qq("methylation:@{chr}"), range = meth_range, ylab = "methylation"))
				# try(mat_dist(current_cov_list, reorder_column = FALSE, od = od, subgroup = subgroup, ha = ha, title = qq("coverage:@{chr}"), range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})")))
			}

			meth_list = lapply(seq_along(current_meth_list), function(i) {
				c(meth_list[[i]], current_meth_list[[i]])
			})
			# cov_list = lapply(seq_along(cov_list), function(i) {
			# 	c(cov_list[[i]], current_cov_list[[i]])
			# })
		}

		if(!by_chr) {
			meth_list = lapply(meth_list, function(meth) {
				n = length(meth)
				if(n > 100000) {
					meth[sample(n, 100000)]
				} else {
					meth
				}
			})
			# cov_list = lapply(cov_list, function(cov) {
			# 	n = length(cov)
			# 	if(n > 100000) {
			# 		cov[sample(n, 100000)]
			# 	}
			# })
			
			od = mat_dist(meth_list, reorder_column = reorder_column, subgroup = subgroup, ha = ha, title = "methylation", range = meth_range, ylab = "methylation")
			# mat_dist(cov_list, reorder_column = FALSE, od = od, subgroup = subgroup, ha = ha, title = "coverage", range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})"))
			
			return(invisible(od))
		}
	} else {
		meth_mat = NULL
		# cov_mat = NULL
		for(chr in chromosome) {

			methylation_hooks$set_chr(chr, verbose = FALSE)

			meth_gr = methylation_hooks$gr
			if(!is.null(background)) {
				mtch = as.matrix(findOverlaps(meth_gr, background))
				ind = unique(mtch[, 1])
			} else {
				ind = seq_len(length(meth_gr))
			}
			
			nr = length(ind)
			if(is.null(p) || by_chr) p = min(c(3000, nr))/nr
			ind = ind[sample(c(FALSE, TRUE), nr, replace = TRUE, prob = c(1-p, p))]
			
			qqcat("random sampled @{length(ind)} sites from @{nr} sites on @{chr} (with p = @{sprintf('%.1e', p)})\n")
			mm = methylation_hooks$meth[ind, sample_id]
			if(!is.null(methylation_hooks$cov)) {
				cm = methylation_hooks$cov[ind, sample_id]
				l = cm == 0
				l[is.na(l)] = TRUE
				mm[l] = NA
			}
			# cm[cm == 0] = NA
			# cm[cm > max_cov] = NA
			
			qqcat("on average there are @{round(mean(apply(mm, 2, function(x) sum(is.na(x)))))} CpG without coverage information.\n")

			meth_mat = rbind(meth_mat, mm)
			# cov_mat = rbind(cov_mat, cm)

			if(by_chr) {
				
				try(od <- mat_dist(mm, reorder_column = reorder_column, subgroup = subgroup, ha = ha, title = qq("methylation:@{chr}"), range = meth_range, ylab = "methylation"))
				# try(mat_dist(log10(cm), reorder_column = FALSE, od = od, subgroup = subgroup, ha = ha, title = qq("coverage:@{chr}"), range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})")))
			}
		}

		if(!by_chr) {
			nr = nrow(meth_mat)
			if(nr > 100000) {
				meth_mat = meth_mat[sample(nr, 100000), ]
				# cov_mat = cov_mat[sample(nr, 100000), ]
			}
			
			od = mat_dist(meth_mat, reorder_column = reorder_column, subgroup = subgroup, ha = ha, title = "methylation", range = meth_range, ylab = "methylation")
			# mat_dist(log10(cov_mat), reorder_column = FALSE, od = od, subgroup = subgroup, ha = ha, title = "coverage", range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})"))
			
			
			return(invisible(od))
		}
	}
}