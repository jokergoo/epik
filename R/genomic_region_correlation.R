###############################################################
# correlation between two sets of genomic regions
# or enrichment of one type of genomic region in the other
###############################################################


# == title
# Correlation between two sets of genomic regions
#
# == param
# -gr_list_1 a list of `GenomicRanges::GRanges` objects, should be a named list, e.g. low methylated regions in different samples.
# -gr_list_2 a list of `GenomicRanges::GRanges` objects, should be a named list, e.g. a list of different genomic features.
# -background a `GenomicRanges::GRanges` object. The correlation is only looked in the background regions.
# -chromosome a vector of chromosome names
# -species species, used for random shuffling genomic regions
# -nperm number of random shufflings. If it is set to 0 or 1, no random shuffling will be performed.
# -mc.cores number of cores for parallel calculation
# -stat_fun method to calculate correlations. There are some pre-defined functions:
#           `genomic_corr_reldist`, `genomic_corr_absdist` measure how two sets of genomic regions are close; `genomic_corr_jaccard`,
#           `genomic_corr_intersect` measures how two sets of genomic regions are overlapped.
#           The self-defined function should accept at least two arguments which are two GRanges object.
#           The third argument is ``...`` which is passed from the main function. The function
#           should only return a numeric value.
# -... pass to ``stat_fun``
# -bedtools_binary random shuffling is perfomed by ``bedtools``. If ``bedtools`` is not in ``PATH``, the path of ``bedtools`` can be set here.
# -tmpdir dir for temporary files
#
# == details
# The correlation between two sets of genomic regions basically means how much the first type of genomic regions
# are overlapped or close to the second type of genomic regions.
#
# The significance of the correlation is calculated by random shuffling the regions. 
# In random shuffling, regions in ``gr_list_1`` will be shuffled. So if you want to shuffle ``gr_list_2``,
# just switch the first two arguments.
#
# Pleast note random shuffling is done by bedtools, so bedtools should be installed and exists in ``PATH``
# and should support ``-i -g -incl`` options.
#
# == value
# A list containing following elements:
#
# -stat statistic value
# -fold_change stat/E(stat), stat divided by expected value which is generated from random shuffling
# -p.value p-value for over correlated. So, 1 - p.value is the significance for being less correlated
# -stat_random_mean mean value of stat in random shuffling
# -stat_random_sd standard deviation in random shuffling
#
# If ``perm`` is set to 0 or 1, ``fold_change``, ``p.value``, ``stat_random_mean`` and ``stat_random_sd`` are all ``NULL``.
#
# == seealso
# `genomic_corr_reldist`, `genomic_corr_jaccard`, `genomic_corr_absdist`, `genomic_corr_intersect`, 
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
# gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
# genomic_regions_correlation(gr1, gr2, nperm = 0)
# genomic_regions_correlation(list(gr1 = gr1), list(gr2 = gr2), nperm = 0)
#
genomic_regions_correlation = function(gr_list_1, gr_list_2, background = NULL,
	chromosome = paste0("chr", c(1:22, "X", "Y")), species = "hg19",
	nperm = 0, mc.cores = 1, stat_fun = genomic_corr_jaccard, ..., 
	bedtools_binary = Sys.which("bedtools"), tmpdir = tempdir()) {
	
	## check input
	if(inherits(gr_list_1, "GRanges")) {
		gr_name_1 = deparse(substitute(gr_list_1))
		gr_list_1 = list(gr_list_1)
		names(gr_list_1) = gr_name_1
	}
	if(inherits(gr_list_2, "GRanges")) {
		gr_name_2 = deparse(substitute(gr_list_2))
		gr_list_2 = list(gr_list_2)
		names(gr_list_2) = gr_name_2
	}
	
	if(is.null(names(gr_list_1))) {
		stop("`gr_list_1` should have names.\n")
	}
	if(is.null(names(gr_list_2))) {
		stop("`gr_list_2` should have names.\n")
	}

	if(nperm > 1) {
		if(!file.exists(bedtools_binary) && nperm > 1) {
			stop("Cannot find binary file for bedtools.")
		}
	}

	qqcat("set strand to * and merge potential overlapped regions in gr_list_1.\n")
	gr_list_1 = lapply(gr_list_1, function(gr) {
		strand(gr) = "*"
		reduce(sort(gr))
	})
	qqcat("set strand to * and merge potential overlapped regions in gr_list_2.\n")
	gr_list_2 = lapply(gr_list_2, function(gr) {
		strand(gr) = "*"
		reduce(sort(gr))
	})
	
	# limit in chromosomes
	qqcat("subset regions in selected chromosomes.\n")
	gr_list_1 = lapply(gr_list_1, function(gr) gr[ seqnames(gr) %in% chromosome])
	gr_list_2 = lapply(gr_list_2, function(gr) gr[ seqnames(gr) %in% chromosome])

	gr_name_1 = names(gr_list_1)
	gr_name_2 = names(gr_list_2)
	
	# limit in background
	if(!is.null(background)) {
		strand(background) = "*"
		background = background[ seqnames(background) %in% chromosome ]
		background = reduce(background)

		qqcat("overlaping `gr_list_1` to background\n")
		gr_list_1 = lapply(gr_list_1, function(gr) {
			intersect(gr, background)
		})
		qqcat("overlaping `gr_list_2` to background\n")
		gr_list_2 = lapply(gr_list_2, function(gr) {
			intersect(gr, background)
		})

		# since `background` will be send to bedtoos, here we need a data frame
		background_df = as.data.frame(background)[1:3]
	}

	# prepare values that will be returned
	foldChange = matrix(0, nrow = length(gr_list_2), ncol = length(gr_list_1))
	rownames(foldChange) = names(gr_list_2)
	colnames(foldChange) = names(gr_list_1)
	p = foldChange
	stat = foldChange
	stat_random_mean = stat
	stat_random_sd = stat

	if(nperm > 1) {
		chr_len_df = getChromInfoFromUCSC(species)
		chr_len_df = chr_len_df[chr_len_df[[1]] %in% chromosome, , drop = FALSE]  # needed for bedtools shuffle
		# cache chr_len_df and background files
		chr_len_df_tmp = tempfile(tmpdir = tmpdir)
		write.table(chr_len_df, file = chr_len_df_tmp, sep = "\t", row.names = FALSE, col.names= FALSE, quote = FALSE)
	}

	if(!is.null(background)) {
		background_df_tmp = tempfile(tmpdir = tmpdir)
		write.table(background_df, file = background_df_tmp, sep = "\t", row.names = FALSE, col.names= FALSE, quote = FALSE)
	}

	# loop in gr_list_1
	for(i in seq_along(gr_list_1)) {

		stat_random = matrix(0, nrow = length(gr_list_2), ncol = nperm)
		
		# stat
		for(j in seq_along(gr_list_2)) {
			qqcat("calculating correlation between @{gr_name_1[i]} and @{gr_name_2[[j]]}\n")
			stat[j, i] = suppressWarnings(do.call("stat_fun", list(gr_list_1[[i]], gr_list_2[[j]], ...)))
		}

		if(nperm > 1) {
			# random shuffle gr_list_1
			# cache gr_list_1
			gr_list_1_df = as.data.frame(gr_list_1[[i]])[1:3]
			gr_list_1_df_tmp = tempfile()
			write.table(gr_list_1_df, file = gr_list_1_df_tmp, sep = "\t", row.names = FALSE, col.names= FALSE, quote = FALSE)

			res = mclapply(seq_len(nperm), mc.cores = mc.cores, function(k) {
				
				if(is.null(background)) {
					gr_random = systemdf(qq("'@{bedtools_binary}' shuffle -i '@{gr_list_1_df_tmp}' -g '@{chr_len_df_tmp}'"))
				} else {
					gr_random = systemdf(qq("'@{bedtools_binary}' shuffle -i '@{gr_list_1_df_tmp}' -g '@{chr_len_df_tmp}' -incl '@{background_df_tmp}'"))
				}

				gr_random = GRanges(seqnames = gr_random[[1]], ranges = IRanges(gr_random[[2]], gr_random[[3]]))

				# x contains stat for every gr_list_2
				x = numeric(length(gr_list_2))
				for(j in seq_along(gr_list_2)) {
					
					qqcat("calculating correlation between random_@{gr_name_1[i]} and @{gr_name_2[[j]]}, @{k}/@{nperm}\n")

					x[j] = suppressWarnings(do.call("stat_fun", list(gr_random, gr_list_2[[j]], ...)))
				}

				return(x)
			})

			# `res` is a list, convert to a matrix
			for(k in seq_along(res)) {
				stat_random[, k] = res[[k]]
			}

			file.remove(gr_list_1_df_tmp)

			stat_random_mean[, i] = rowMeans(stat_random)
			stat_random_sd[, i] = rowSds(stat_random)

			foldChange[, i] = stat[, i]/stat_random_mean[, i]
			p[, i] = sapply(seq_len(length(gr_list_2)), function(i) sum(stat_random[i, ] > stat[i])/nperm)
		} else {
			stat_random_mean[, i]  = NA
			stat_random_sd[, i] = NA
			foldChange[, i] = NA
			p[, i] = NA
		}
	}

	if(nperm > 1) {
		file.remove(chr_len_df_tmp)
	}
	if(!is.null(background)) {
		file.remove(background_df_tmp)
	}

	if(nperm <= 1) {
		foldChange = NULL
		p = NULL
		stat_random_mean = NULL
		stat_random_sd = NULL
	}

	res = list(stat = stat,
		       fold_change = foldChange, 
		       p.value = p,
		       stat_random_mean = stat_random_mean,
		       stat_random_sd = stat_random_sd)
	
	return(res)
}

# == title
# Relative distance between two sets of genomic regions
#
# == param
# -gr1 genomic region 1, a `GenomicRanges::GRanges` object
# -gr2 genomic region 2, a `GenomicRanges::GRanges` object
#
# == details
# For regions in ``gr1`` and ``gr2``, they are all degenerated as single points
# which are the middle points of regions. For each middle point in ``gr1``, it looks 
# for two nearest points in ``gr2`` on its left and right. The statistic is defined as the ratio of the distance
# to the nearest neighbour point to the distance of two neighbour points. If ``gr1`` and ``gr2`` are not correlated at all,
# It is expected that the ratio follows a uniform distribution. So final statisitic are the KS-statistic
# between the real distribution of rations to the uniform distribution.
#
# == reference
# Favoriv A, et al. Exploring massive, genome scale datasets with the GenometriCorr package. PLoS Comput Biol. 2012 May; 8(5):e1002529
# 
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_reldist(gr1, gr2)
#
genomic_corr_reldist = function(gr1, gr2) {
	# GRanges for mid-points
	gr1_mid = ceiling( (start(gr1) + end(gr1))/2 )
	gr2_mid = ceiling( (start(gr2) + end(gr2))/2 )
	
	# we don't need strand information here
	gr1_mid_gr = GRanges(seqnames = seqnames(gr1),
		                   ranges = IRanges(start = gr1_mid,
		                   	                end = gr1_mid))
	gr2_mid_gr = GRanges(seqnames = seqnames(gr2),
		                   ranges = IRanges(start = gr2_mid,
		                   	                end = gr2_mid))
	
	gr2_mid_gr = sort(gr2_mid_gr)
	
	# look for nearest position
	mtch = nearest(gr1_mid_gr, gr2_mid_gr) # index in gr2
	l = !is.na(mtch)        # if there is no site in gr2 on some chromosome, the value would be NA
	gr1_mid_gr = gr1_mid_gr[l]   # remove these regions in gr1
	m2 = gr2_mid_gr[ mtch[l] ]
	mtch = mtch[l]
	
	# look for another point in gr2
	ind = ifelse(start(gr1_mid_gr) > start(m2), mtch + 1, mtch - 1)  # ind contains index in gr2
	l = ind >= 1 & ind <= length(gr2)
	m3 = gr2_mid_gr[ ind[l] ]
	gr1_mid_gr = gr1_mid_gr[l]
	m2 = m2[l]

	suppressWarnings(d <- distance(gr1_mid_gr, m2) / distance(m2, m3))
	d = d[!is.na(d) & !is.infinite(d)]  # this is not necessary, just double check
	
	stat = 0
	try(stat <- (integrate(ecdf(d), lower = 0, upper = 0.5)$value - 0.25) / 0.25, silent = TRUE)
	
	return(stat)
}

# == title
# Jaccard coefficient between two sets of genomic regions
#
# == param
# -gr1 genomic region 1, a `GenomicRanges::GRanges` object
# -gr2 genomic region 2, a `GenomicRanges::GRanges` object
# -background background regions that should be only looked in, a `GenomicRanges::GRanges` object
#
# == details
# Jaccard coefficient is defined as the total length of intersection divided by total
# length of union of two sets of genomic regions.
#
# You can set the background when calculating Jaccard coefficient. For example,
# if the interest is the Jaccard coefficient between CpG sites in ``gr1`` and in ``gr2``
# ``background`` can be set with a `GenomicRanges::GRanges` object which contains positions of CpG sites.
#
# Be careful with the ``strand`` in your `GenomicRanges::GRanges` object!
#
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_jaccard(gr1, gr2)
#
genomic_corr_jaccard = function(gr1, gr2, background = NULL) {
	if(is.null(background)) {
		res = sum(as.numeric(width(intersect(gr1, gr2)))) / sum(as.numeric(width(union(gr1, gr2))))
	} else {
		gr1 = intersect(gr1, gr2)
		gr1 = intersect(gr1, background)

		gr2 = union(gr1, gr2)
		gr2 = intersect(gr2, background)
		res = sum(as.numeric(width(gr1))) / sum(as.numeric(width(gr2)))
	}
	return(res)
}

# == title
# Absolute distance between two sets of genomic regions
#
# == param
# -gr1 genomic region 1, a `GenomicRanges::GRanges` object
# -gr2 genomic region 2, a `GenomicRanges::GRanges` object
# -method function in which input is a vector of distance and output is a scalar
# -... pass to ``method``
#
# == details
# For regions in ``gr1`` and ``gr2``, they are all degenerated as single points
# which are the middle points of corresponding regions. For each middle point in ``gr1``, it looks 
# for the nearest point in ``gr2``. Assuming the distance vector is ``d``, the final statistic is ``method(d)``.
#
# == reference
# Favoriv A, et al. Exploring massive, genome scale datasets with the GenometriCorr package. PLoS Comput Biol. 2012 May; 8(5):e1002529
#
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_absdist(gr1, gr2)
#
genomic_corr_absdist = function(gr1, gr2, method = mean, ...) {
	# GRanges for mid-points
	gr1_mid = ceiling( (start(gr1) + end(gr1))/2 )
	gr2_mid = ceiling( (start(gr2) + end(gr2))/2 )
	
	gr1_mid_gr = GRanges(seqnames = seqnames(gr1),
		                   ranges = IRanges(start = gr1_mid,
		                   	                end = gr1_mid))
	gr2_mid_gr = GRanges(seqnames = seqnames(gr2),
		                   ranges = IRanges(start = gr2_mid,
		                   	                end = gr2_mid))
	
	# look for nearest position
	suppressWarnings(mtch <- distanceToNearest(gr1_mid_gr, gr2_mid_gr))
	dst = mtch@elementMetadata[, 1]
	dst = dst[!is.na(dst)]
	stat = method(as.numeric(dst), ...)
	
	return(stat)
}

# == title
# Intersections between two sets of genomic regions
#
# == param
# -gr1 genomic region 1, a `GenomicRanges::GRanges` object
# -gr2 genomic region 2, a `GenomicRanges::GRanges` object
# -method how to calculate the intersection statistic, see "details"
# -... pass to `GenomicRanges::countOverlaps` or `percentOverlaps`
#
# == details
# There are three metrics for the intersection statistic:
#
# -number It calculates number of regions in ``gr1`` that overlap with ``gr2``.
#       Please note this value is not equal to the number of intersections betweenn two sets of regions,
#       because one region in ``gr1`` may overlap with more than one
#       regions in ``gr2``.
# -percent It calculates for each region in ``gr1``, how much it is covered by regions in ``gr2``.
# -length sum of length of the intersection of the two sets of regions.
#
# With methods of"number" and "percent", ``genomic_corr_intersect(gr1, gr2)`` is always not identical
# to ``genomic_corr_intersect(gr2, gr1)``.
#
# == value
# A single correlation value.
#
# == seealso
# `genomic_regions_correlation`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr1 = GRanges(seqnames = "chr1", ranges = IRanges(c(1, 5), c(3, 8)))
# gr2 = GRanges(seqnames = "chr1", ranges = IRanges(c(2, 6), c(4, 8)))
# genomic_corr_intersect(gr1, gr2, method = "number")
# genomic_corr_intersect(gr1, gr2, method = "percent")
# genomic_corr_intersect(gr1, gr2, method = "length")
genomic_corr_intersect = function(gr1, gr2, method = c("number", "percent", "length"), ...) {
	method = match.arg(method)[1]

	if(method == "number") {
		x = countOverlaps(gr1, gr2, ...)
		res = sum(x > 0)
	} else if(method == "percent") {
		x = percentOverlaps(gr1, gr2, ...)
		w = width(gr1)
		res = sum(x*w)/sum(as.numeric(w))
	} else if(method == "length") {
		x = intersect(gr1, gr2)
		res = sum(as.numeric(width(x)))
	}
	return(res)
}

