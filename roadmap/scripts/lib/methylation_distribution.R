
# -x a data frame or a list
.mat_dist = function(x, anno, col, reorder_column = TRUE, od = seq_along(anno), ha = NULL, title = NULL, ...) {

	if(reorder_column) {
		if(inherits(x, "list")) {
			od = order(factor(anno, levels = unique(anno), ordered = TRUE), sapply(x, median, na.rm = TRUE))
		} else {
			od = order(factor(anno, levels = unique(anno), ordered = TRUE), apply(x, 2, median, na.rm = TRUE))
		}
	} 
	if(is.null(ha)) ha = HeatmapAnnotation(df = data.frame(type = anno), col = list(type = col))

	densityHeatmap(x, anno = ha, title = title, column_order = od, ...)
	for(an in sapply(ha@anno_list, function(x) x@name)) {
		decorate_annotation(an, {
			grid.text(an, x = unit(-2, "mm"), just = "right")
		})
	}

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

	## MDS plot
	if(is.data.frame(x) || is.matrix(x)) {
		mat = as.matrix(x)
		loc = cmdscale(dist2(t(mat), pairwise_fun = function(x, y) {l = is.na(x) | is.na(y); x = x[!l]; y = y[!l]; sqrt(sum((x-y)^2))}))
		plot(loc[, 1], loc[, 2], pch = 16, cex = 1, col = col[anno], main = qq("MDS:@{title}"), xlab = "dimension 1", ylab = "dimension 2")
		legend("bottomleft", pch = 16, legend = names(col), col = col)

		plot(loc[, 1], loc[, 2], type = "n", pch = 16, cex = 1, col = col[anno], main = qq("MDS:@{title}"), xlab = "dimension 1", ylab = "dimension 2")
		text(loc[, 1], loc[, 2], colnames(x), col = col[anno], cex = 0.8)
	}

	return(od)
}


# == title
# Global methylation distribution
# 
# == param
# -sample_id a vector of sample ids
# -annotation subtype information
# -annotation_color color for subtypes
# -reorder_column whether reorder the samples
# -ha additional annotation can be specified as a `ComplexHeatmap::HeatmapAnnotation` object
# -chromosome chromosomes
# -by_chr whether make the plot by chromosome
# -max_cov maximum coverage (used to get rid of extremely high coverage which affects visualization of CpG coverage distribution)
# -background background to look into. The value can be a single `GenomicRanges::GRanges` object or a list of `GenomicRanges::GRanges` objects.
# -p probability to randomly sample CpG sites
#
# == details
# It visualize distribution of methylation values and CpG coverages through heatmaps.
#
# == value
# If ``by_chr`` is set to ``FALSE``, it returns a vector of column order.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
global_methylation_distribution = function(sample_id, annotation, 
	annotation_color = structure(seq_along(unique(annotation)), names = unique(annotation)),
	reorder_column = TRUE, ha = NULL, chromosome = paste0("chr", 1:22), by_chr = FALSE, max_cov = 100,
	background = NULL, p = 0.001, meth_range = c(0, 1)) {

	annotation_color = annotation_color[intersect(names(annotation_color), unique(annotation))]
	
	###############################################
	# distribution of global methylation
	if(inherits(background, "list")) {
		meth_list = NULL
		cov_list = NULL

		if(length(background) != length(sample_id)) {
			stop("Since you specified `background` as a list, the length should be same as `sample_id`.")
		}

		for(chr in chromosome) {

			methylation_hooks$set(chr)

			meth_gr = methylation_hooks$GRanges()
			ind_list = lapply(seq_along(sample_id), function(i) {
				mtch = as.matrix(findOverlaps(meth_gr, background[[i]]))
				ind = unique(mtch[, 1])
				nr = length(ind)
				ind = ind[sample(c(FALSE, TRUE), nr, replace = TRUE, prob = c(1-p, p))]

				message(qq("random sampled @{length(ind)} sites from @{nr} sites on @{chr} in @{sample_id[i]} (with p = @{p})"))
			})

			current_meth_list = lapply(seq_along(ind_list), function(i) {
				m = methylation_hooks$meth(row_index = ind_list[[i]], col_index = sample_id[i])[, 1]
				cov = methylation_hooks$coverage(row_index = ind_list[[i]], col_index = sample_id[i])[, 1]
				m[cov == 0] = NA
				m
			})
			current_cov_list = lapply(seq_along(ind_list), function(i) {
				cov = methylation_hooks$coverage(row_index = ind_list[[i]], col_index = sample_id[i])[, 1]
				cov[cov == 0] = NA
				cov[cov > max_cov] = NA
				log10(cov)
			})

			message(qq("on average there are @{round(mean(sapply(current_meth_list, function(x) sum(is.na(x)))))} CpG with 0 coverage.\n"))

			
			if(by_chr) {
				try(od <- .mat_dist(current_meth_list, reorder_column = reorder_column, anno = annotation, col = annotation_color, ha = ha, title = qq("methylation:@{chr}"), range = meth_range, ylab = "methylation"))
				try(.mat_dist(current_cov_list, reorder_column = FALSE, od = od, anno = annotation, col = annotation_color, ha = ha, title = qq("coverage:@{chr}"), range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})")))
			}

			meth_list = lapply(seq_along(meth_list), function(i) {
				c(meth_list[[i]], current_meth_list[[i]])
			})
			cov_list = lapply(seq_along(cov_list), function(i) {
				c(cov_list[[i]], current_cov_list[[i]])
			})
		}

		if(!by_chr) {
			meth_list = lapply(meth_list, function(meth) {
				n = length(meth)
				if(n > 100000) {
					meth[sample(n, 100000)]
				}
			})
			cov_list = lapply(cov_list, function(cov) {
				n = length(cov)
				if(n > 100000) {
					cov[sample(n, 100000)]
				}
			})
			
			od = .mat_dist(meth_list, reorder_column = reorder_column, anno = annotation, col = annotation_color, ha = ha, title = "methylation", range = meth_range, ylab = "methylation")
			.mat_dist(cov_list, reorder_column = FALSE, od = od, anno = annotation, col = annotation_color, ha = ha, title = "coverage", range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})"))
			
			
			return(invisible(od))
		}
	} else {
		meth_mat = NULL
		cov_mat = NULL
		for(chr in chromosome) {

			methylation_hooks$set(chr)

			meth_gr = methylation_hooks$GRanges()
			if(!is.null(background)) {
				mtch = as.matrix(findOverlaps(meth_gr, background))
				ind = unique(mtch[, 1])
			} else {
				ind = seq_len(length(meth_gr))
			}
			
			nr = length(ind)
			ind = ind[sample(c(FALSE, TRUE), nr, replace = TRUE, prob = c(1-p, p))]
			
			message(qq("random sampled @{length(ind)} sites from @{nr} sites on @{chr} (with p = @{p})"))
			mm = methylation_hooks$meth(row_index = ind, col_index = sample_id)
			cm = methylation_hooks$coverage(row_index = ind, col_index = sample_id)
			mm[cm == 0] = NA
			cm[cm == 0] = NA
			cm[cm > max_cov] = NA
			
			message(qq("on average there are @{round(mean(apply(mm, 2, function(x) sum(is.na(x)))))} CpG with 0 coverage.\n"))

			meth_mat = rbind(meth_mat, mm)
			cov_mat = rbind(cov_mat, cm)

			if(by_chr) {
				
				try(od <- .mat_dist(mm, reorder_column = reorder_column, anno = annotation, col = annotation_color, ha = ha, title = qq("methylation:@{chr}"), range = meth_range, ylab = "methylation"))
				try(.mat_dist(log10(cm), reorder_column = FALSE, od = od, anno = annotation, col = annotation_color, ha = ha, title = qq("coverage:@{chr}"), range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})")))
			}
		}

		if(!by_chr) {
			nr = nrow(meth_mat)
			if(nr > 100000) {
				meth_mat = meth_mat[sample(nr, 100000), ]
				cov_mat = cov_mat[sample(nr, 100000), ]
			}
			
			od = .mat_dist(meth_mat, reorder_column = reorder_column, anno = annotation, col = annotation_color, ha = ha, title = "methylation", range = meth_range, ylab = "methylation")
			.mat_dist(log10(cov_mat), reorder_column = FALSE, od = od, anno = annotation, col = annotation_color, ha = ha, title = "coverage", range = c(0, log10(max_cov)), ylab = qq("log10(CpG coverage, 1~@{max_cov})"))
			
			
			return(invisible(od))
		}
	}
}

get_mean_methylation_in_genomic_features = function(sample_id, gf_list, average = TRUE, p = 0.001,
	chromosome = paste0("chr", 1:22),
	filter_fun = function(s) length(s) >= 10 && any(diff(s) < 50)) {
	
	# initialize the mat_list
	# mat_list has the same name as `gf_list`
	
	gr_list = rep(list(GRanges(seqnames = character(0), ranges = IRanges(start = integer(0), end = integer(0)))), 
		          length(gf_list))
	names(gr_list) = names(gf_list)

	for(chr in chromosome) {
		methylation_hooks$set(chr)
		meth_mat = methylation_hooks$meth(col_index = sample_id)
		meth_gr = methylation_hooks$GRanges()
		meth_site = methylation_hooks$site()
		
		for(i in seq_along(gf_list)) {
			qqcat("overlapping to @{names(gf_list)[i]} on @{chr}\n")
			mtch = as.matrix(findOverlaps(gf_list[[i]], meth_gr))  ## <- test memory usage
			if(average){
				l = tapply(mtch[,2], mtch[,1], function(i) {
						x = meth_site[i]
						filter_fun(x)
					})
				mean_meth = tapply(mtch[,2], mtch[,1], function(i) colMeans(meth_mat[i, , drop = FALSE], na.rm = TRUE))
				mean_meth = mean_meth[l]

				ncpg = tapply(mtch[,2], mtch[,1], length)

				mean_meth_mat = matrix(unlist(mean_meth), nrow = length(mean_meth), byrow = TRUE)
				rownames(mean_meth_mat) = names(mean_meth); colnames(mean_meth_mat) = sample_id
				ind = as.integer(names(mean_meth))
				gr = gf_list[[i]][ind]
				mcols(gr) = cbind(as.data.frame(mean_meth_mat), ncpg = ncpg[names(mean_meth)])
				l = apply(mean_meth_mat, 1, function(x) sum(is.na(x))/length(x)) < 0.1
				gr = gr[l]
			} else {
				l = unique(mtch[, 2])
				gr = meth_gr[l]
				mcols(gr) = as.data.frame(meth_mat[l, , drop = FALSE])
				if(length(gr) > 10000) gr = gr[sample(c(TRUE, FALSE), length(gr), prob = c(p, 1-p), replace = TRUE)]
				gc(verbose = FALSE)
			}
			gr_list[[i]] = c(gr_list[[i]], gr)
		}
	}

	attr(gr_list, "generated_by") = "get_mean_methylation_in_genomic_features"

	return(gr_list)
}
