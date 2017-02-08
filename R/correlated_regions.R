
# == title
# Correlated regions in a specified region
#
# == param
# -site position of CpG sites in this region, should be sorted
# -meth methylation matrix corresponding to ``site``
# -expr expression for the associated gene
# -chr chromosome name, used to construct the `GenomicRanges::GRanges` object
# -cov CpG coverage matrix. CpG coverage is important when ``meth`` is the raw methylation which means
#      CpG sites with extremely low coverage will be removed when calculating correlations
# -cov_cutoff cutoff for the CpG coverage
# -min_dp minimal number of non-NA values for calculating correlations. When ``meth`` is the raw methylation,
#          values for which CpG coverage is too low will be replaced with ``NA``, We only use non-NA values
#          to calculate correlations.
# -cor_method method for calcualting correlations
# -window_size how many CpG sites in a window
# -window_step step of the sliding window, measured in number of CpG sites
# -subgroup subgroup information
# -max_width maximum width of a window
#
# == details
# ``cov`` and ``cov_cutoff`` should be set when the methylation is unsmoothed, because
# for the unsmoothed data, the methylation rate is not unreliable when the CpG coverage is low.
#
# == value
# a `GenomicRanges::GRanges` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
correlated_regions_by_window = function(site, meth, expr, chr, cov = NULL, cov_cutoff = 3, min_dp = 4,
	cor_method = "spearman", window_size = 5, window_step = window_size, subgroup = NULL, max_width = 10000) {

	if(ncol(meth) != length(expr)) {
		stop("number of columsn of `meth` should be same as length of `expr`.\n")
	}

	index = seq(1, length(site) - window_size, by = window_step)  # just forget last few sites
	
	i = seq_len(length(index))
	ir = IRanges(site[index[i]], site[index[i]+window_size-1])
	l = width(ir) <= max_width
	i = i[l]
	ir = ir[l]

	if(length(i) == 0) {
		message(qq("no window left by max_width <= @{max_width} bp"))
		return(GRanges())
	}

	if(!is.null(cov)) {
		meth[cov < cov_cutoff] = NA
	}

	m = lapply(index[i], function(x) {
		ind = x+0:(window_size-1)
		meth_m = meth[ind, , drop = FALSE]
		sapply(seq_len(ncol(meth_m)), function(i) {
			xm = meth_m[, i]
			mean(xm, na.rm = TRUE)
		})
	})
	m = do.call('rbind', m)
	colnames(m) = paste0("mean_meth_", colnames(meth))
	corr = apply(m, 1, function(x) {
		l = !is.na(x)
		if(sum(l) < min_dp) return(NA)
		cor(x[l], expr[l], method = cor_method)
	})
	corr_p = suppressWarnings(apply(m, 1, function(x) {
		l = !is.na(x)
		if(sum(l) < min_dp) return(NA)
		cor.test(x[l], expr[l], method = cor_method)$p.value
	}))

	if(!is.null(subgroup)) {
		if(length(unique(subgroup)) == 1) subgroup = NULL
	}

	if(!is.null(subgroup)) {
		subgroup = as.vector(subgroup)
		meth_anova = apply(m, 1, function(x) {
			l = !is.na(x)
			data = data.frame(value = x[l], class = subgroup[l], stringsAsFactors = FALSE)
			if(length(unique(data$class)) < 2) return(NA)
			if(any(table(data$class) < 2)) return(NA)
			oneway.test(value ~ class, data = data)$p.value
		})
		meth_diameter = apply(m, 1, function(x) {
			l = !is.na(x)
			if(any(table(subgroup[l]) < 2)) return(NA)
			diameter(as.vector(tapply(x[l], subgroup[l], mean)))
		})
		if(length(unique(subgroup)) == 2) {
			unique_subgroup = sort(unique(subgroup))
			meth_diff = apply(m, 1, function(x) {
				l = !is.na(x)
				if(any(table(subgroup[l]) < 2)) return(NA)
				y = tapply(x[l], subgroup[l], mean)
				y[unique_subgroup[1]] - y[unique_subgroup[2]]
			})
		}
	}
	if(nrow(m) == 1) {
		meth_IQR = IQR(m, na.rm = TRUE)
	} else {
		meth_IQR = rowIQRs(m, na.rm = TRUE)
	}
	gr = GRanges(seqnames = rep(chr, length(ir)), ranges = ir)

	if(is.null(subgroup)) {
		df = DataFrame(ncpg = window_size,
		    m, # mean methylation
		    corr = corr,
		    corr_p = corr_p,
		    meth_IQR = meth_IQR)	
	} else {
		df = DataFrame(ncpg = window_size,
		    m, # mean methylation
		    corr = corr,
		    corr_p = corr_p,
		    meth_IQR = meth_IQR,
		    meth_anova = meth_anova,
		    meth_diameter = meth_diameter,
		    meth_diff = meth_diff)
	}
	mcols(gr) = df

	return(gr)
}

# == title
# Correlation between methylation and expression
#
# == param
# -sample_id a vector of sample ids
# -expr expression matrix
# -txdb a `GenomicFeatures::TxDb-class` object.
# -chr a single chromosome name
# -extend extension of gene model, both upstream and downstream
# -cov_filter if ``coverage`` hook is set in `methylation_hooks`, this option can be set to filter out CpG sites with low coverage across samples.
#          the value for this option is a function for which the argument is a vector of coverage values for current CpG in all samples.
# -cor_method method for calcualting correlations
# -subgroup subgroup information
# -window_size how many CpG sites in a window
# -window_step step of the sliding window, measured in number of CpG sites
# -max_width maximum width of a window
# -raw_meth whether use raw methylation value (values from ``raw`` hook set in `methylation_hooks`)
# -cov_cutoff cutoff for CpG coverage
# -min_dp minimal number of non-NA values for calculating correlations. When ``meth`` is the raw methylation,
#          values for which CpG coverage is too low will be replaced with ``NA``, We only use non-NA values
#          to calculate correlations.
# -col color for subgroups. This setting will be saved in the returned object and will be used in downstream analysis.
#      If not set, random colors are assigned.
# -species species. This setting will be used in downstream analysis
#
# == details
# The detection for correlated regions is gene-centric. For every gene, the process are as follows:
#
# - extend to both upstream and downstream by ``extend``
# - from the most upstream, use a sliding window which contains ``windows_size`` CpG sites, moving step of ``window_step`` CpG sites;
# - filter each window by CpG coverage (by ``cov_filter`` and ``cov_cutoff``) if ``raw_meth`` is ``TRUE``
# - calculate correlation between methylation and gene expression for this window
# - calculate other statistics
#
# Following meta columns are attached to the `GenomicRanges::GRanges` objects:
#
# -ncpg number of CpG sites
# -mean_meth_* mean methylation in each window in every sample.
# -corr correlation between methylation and expression
# -corr_p p-value for the correlation test
# -meth_IQR IQR of mean methylation if ``subgroup`` is not set
# -meth_anova p-value from oneway ANOVA test if ``subgroup`` is set
# -meth_diameter range between maximum mean and minimal mean in all subgroups if ``subgroup`` is set
# -meth_diff when there are two subgroups, the mean methylation in subgroup 1 substracting mean methylation in subgroup 2.
# -gene_id gene id
# -gene_tss_dist distance to tss of genes
# -tx_tss_dist if genes have multiple transcripts, this is the distance to the nearest transcript
# -nearest_txx_tss transcript id of the nearest transcript
#
# This function keeps all the information for all CpG windows. Uses can use `cr_add_fdr_column` to add fdr columns to the object,
# filter significant correlated regions by p-value, fdr and meth_diff columns, or use `cr_reduce` to reduce the significant regions.
#
# == value
# A `GenomicRanges::GRanges` object which contains correlations and associated statistics for every CpG windows.
# 
# The settings for finding correlated regions are stored as metadata of the `GenomicRanges::GRanges` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
correlated_regions = function(sample_id, expr, txdb, chr, extend = 50000,
	cov_filter = function(x) sum(x > 0, na.rm = TRUE) > length(x)/2,
	cor_method = "spearman", subgroup = NULL, window_size = 5, window_step = window_size, 
	max_width = 10000, raw_meth = FALSE, cov_cutoff = 3, min_dp = 4, col = NULL, species = "hg19") {

	message(qq("extracting gene model (extend = @{extend}, chr = @{chr})..."))
	gene = genes(txdb)
	tx_list = transcriptsBy(txdb, by = "gene")

	gene = gene[seqnames(gene) == chr]

	g = intersect(rownames(expr), names(gene))
	expr = expr[g, , drop = FALSE]
	gene = gene[g]
	tx_list = tx_list[g]

	genemodel = gene
	start(genemodel) = start(genemodel) - extend
	end(genemodel) = end(genemodel) + extend
	s = start(genemodel)
	start(genemodel) = ifelse(s > 0, s, 1)

	all_gi = rownames(expr)
	n_gene = length(all_gi)

	methylation_hooks$set_chr(chr, verbose = FALSE)
	meth_gr = methylation_hooks$gr
	site = start(meth_gr)

	# sample_id = c(sample_id, methylation_hooks$sample_id)
	sample_id = intersect(sample_id, colnames(expr))

	if(!is.null(subgroup)) {
		lt = format_sample_id_and_subgroup(sample_id, subgroup)
		sample_id = lt$sample_id
		subgroup = lt$subgroup
	}
	message(qq("@{length(sample_id)} samples are used"))
	
	if(raw_meth) {
		meth = methylation_hooks$raw[, sample_id]
	} else {
		meth = methylation_hooks$meth[, sample_id]
	}
	cov = methylation_hooks$cov[, sample_id]

	expr = expr[, sample_id, drop = FALSE]

	if(!is.null(cov_filter)) {
		
		l = apply(cov, 1, cov_filter)
		if(any(is.na(l))) {
			stop("`cov_filter` generates `NA`, check it.")
		}
		site = site[l]
		meth = meth[l, , drop = FALSE]
		cov = cov[l, , drop = FALSE]
	}

	message(qq("on chromosome @{chr}, @{length(site)} CpG sites, @{n_gene} genes"))

	if(!raw_meth) cov_cutoff = 0
	
	res = GRanges()
	for(i in seq_len(n_gene)) {
			
		# current gene name
		gi = all_gi[i]

		# expression of current gene
		e = expr[gi, ]

		# if gene has low expression in many samples
		if(all(e == 0) || sd(e) == 0)  {
			message(qq("[@{chr}:@{gi}, @{i}/@{n_gene}] @{gi} has zero expression in all samples, skip"))
			next
		}

		start = start(genemodel[gi])
		end = end(genemodel[gi])
		gm_site_index = extract_sites(start, end, site, TRUE, 0)
		gm_site = site[gm_site_index]
		gm_meth = meth[gm_site_index, , drop = FALSE]
		gm_cov = cov[gm_site_index, , drop = FALSE]

		if(length(gm_site) < 10) {
			message(qq("[@{chr}:@{gi}, @{i}/@{n_gene}] @{gi} has too few cpg sites, skip"))
			next
		}

		message(qq("[@{chr}:@{gi}, @{i}/@{n_gene}] @{length(gm_site)} CpG sites ..."))
		gr = correlated_regions_by_window(gm_site, gm_meth, e, gm_cov, cov_cutoff = cov_cutoff, chr = chr,
			subgroup = subgroup, cor_method = cor_method, window_size = window_size, window_step = window_step,
			min_dp = min_dp, max_width = max_width)
		gr$gene_id = rep(gi, length(gr))

		## distance to gene tss
		tss = promoters(gene[gi], upstream = 1, downstream = 0)
		gene_tss_dist = as.data.frame(distanceToNearest(gr, tss))[, 3]
		if(as.vector(strand(gene[gi])) == "+") {
			gene_tss_dist = ifelse(end(gr) < start(tss), -gene_tss_dist, gene_tss_dist)
		} else if(as.vector(strand(gene[gi])) == "-") {
			gene_tss_dist = ifelse(start(gr) > end(tss), -gene_tss_dist, gene_tss_dist)
		}
		gr$gene_tss_dist = gene_tss_dist

		## distance to tx tss
		tx = tx_list[[gi]]
		tx = tx[tx$tx_name != gi]
		tx_width = width(tx)
		tss = promoters(tx, upstream = 1, downstream = 0)
		dist = as.data.frame(distanceToNearest(gr, tss, select = "all"))
		dist = data.frame(queryHits = as.vector(tapply(dist[[1]], dist[[1]], unique)),
			              subjectHits = as.vector(tapply(dist[[2]], dist[[1]], function(x) x[which.max(tx_width[x])[1]])),
			              distance = as.vector(tapply(dist[[3]], dist[[1]], unique)))
		tx_tss_dist = dist[, 3]
		if(as.vector(strand(gene[gi])) == "+") {
			tx_tss_dist = ifelse(end(gr) < start(tss[dist[,2]]), -tx_tss_dist, tx_tss_dist)
		} else if(as.vector(strand(gene[gi])) == "-") {
			tx_tss_dist = ifelse(start(gr) > end(tss[dist[,2]]), -tx_tss_dist, tx_tss_dist)
		}
		gr$tx_tss_dist = tx_tss_dist
		gr$nearest_tx_tss = tss[dist[,2]]$tx_name

		if(length(gr)) {
			res = c(res, gr)
		}
	}

	param = list()
	param$subgroup = subgroup
	param$cor_method = cor_method
	param$extend = extend
	param$window_size = window_size
	param$window_step = window_step
	param$max_width = max_width
	param$sample_id = sample_id
	param$cov_filter = cov_filter
	param$cov_cutoff = cov_cutoff
	param$raw_meth = raw_meth
	param$min_dp = min_dp
	if(is.null(col)) {
		if(is.null(subgroup)) {
			col = rand_color(length(sample_id))
		} else {
			n_subgroup = length(unique(subgroup))
			col = structure(rand_color(n_subgroup), names = unique(subgroup))
		}
	}
	param$col = col
	param$species = species

	metadata(res) = list(cr_param = param)
	return(res)
}


# == title
# Calcualte FDRs for CRs
#
# == param
# -cr original correlated regions from `correlated_regions`
# -fdr_method method to calculate FDR
#
# == details
# Since correlated region detection is per-chromosome, after merging correlated regions from all chromosomes, FDR
# can be calcualted based on ``corr_p`` and/or ``meth_anno`` column.
#
# Please note, FDRs are calculated for negative CRs and positive CRs separatedly.
#
# == value
# Correlated regions with two/one columns (``corr_fdr``, ``meth_anova_fdr``)
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
cr_add_fdr_column = function(cr, fdr_method = "BH") {

	# separate neg and pos
	l_neg = cr$corr < 0
	l_pos = cr$corr > 0

	# add corr_fdr column
	cr$corr_fdr = 1
	cr$corr_fdr[l_neg] = p.adjust(cr$corr_p[l_neg], fdr_method)
	cr$corr_fdr[l_pos] = p.adjust(cr$corr_p[l_pos], fdr_method)

	if(!is.null(cr$meth_anova)) {
		cr$meth_anova_fdr = 1
		cr$meth_anova_fdr[l_neg] = p.adjust(cr$meth_anova[l_neg], fdr_method)
		cr$meth_anova_fdr[l_pos] = p.adjust(cr$meth_anova[l_pos], fdr_method)
	}

	cr_param = metadata(cr)$cr_param
	cr_param$fdr_method = fdr_method
	metadata(cr)$cr_param = cr_param
	cr
}

# gr is sorted
reduce_cr_by_gene = function(gr, gap = 1) {
	if(length(gr) == 0) return(GRanges())

	window_size = metadata(gr)$cr_param$window_size
	window_step = metadata(gr)$cr_param$window_step

	n = length(gr)
	if(all(start(gr)[-1] > end(gr)[-n]+1)) {
		gr2 = gr
		mcols(gr2) = NULL
		# gr2$mean_corr = gr$corr
		gr2$merged_windows = rep(1, n)
		gr2$ncpg = window_size
	} else {
		gr2 = reduce(gr, min.gap = gap, with.revmap = TRUE)
		revmap = mcols(gr2)$revmap
		mcols(gr2) = NULL
		# gr2$mean_corr = sapply(revmap, function(ind) mean(gr[ind]$corr, na.rm = TRUE))
		gr2$merged_windows = sapply(revmap, length)
		gr2$ncpg = gr2$merged_windows*window_size - (gr2$merged_windows-1)*(window_size - window_step)
	}
	message(qq("  @{deparse(substitute(gr))} has been reduced from @{length(gr)} to @{length(gr2)}"))
	return(gr2)
}

# == title
# Merge correlated regions
#
# == param
# -cr correlated regions from `correlated_regions`. In most cases, it is correlated regions with significant correlations.
# -txdb the transcriptome annotation which is same as the one used in `correlated_regions`
# -expr the expression matrix which is same as the one used in `correlated_regions`. If it is set
#        the correlation will be re-calculated for the merged regions.
# -gap gap for the merging, pass to `reduce2`
# -mc.cores cores for parallel computing. It is paralleled by gene
#
# == details
# As there are overlaps between two neighbouring correlated regions with the default settings, it is possible to merge them into
# large regions. The mering is gene-wise, and all statistics (e.g. mean methylation, correlation) will be
# re-calculated.
#
# == value
# A `GenomicRanges::GRanges` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
cr_reduce = function(cr, txdb, expr = NULL, gap = bp(1), mc.cores = 1) {

	cr_param = metadata(cr)$cr_param
	cr = sort(cr)

	message("extracting gene and tx models.")
	gene = genes(txdb)
	tx_list = transcriptsBy(txdb, by = "gene")

	cr_list = split(cr, cr$gene_id)
	i = 0
	cr_reduced_list = mclapply(names(cr_list), function(gi) {

		message(qq("reducing cr on @{gi}, @{i <<- i+1}/@{length(cr_list)}..."))
		cr = cr_list[[gi]]
		neg_cr = cr[cr$corr < 0]
		pos_cr = cr[cr$corr > 0]
		gr = c(reduce_cr_by_gene(neg_cr, gap = gap),
			   reduce_cr_by_gene(pos_cr, gap = gap))
		gr = sort(gr)
		gr$gene_id = rep(gi, length(gr))

		## distance to gene tss
		tss = promoters(gene[gi], upstream = 1, downstream = 0)
		gene_tss_dist = as.data.frame(distanceToNearest(gr, tss))[, 3]
		if(as.vector(strand(gene[gi])) == "+") {
			gene_tss_dist = ifelse(end(gr) < start(tss), -gene_tss_dist, gene_tss_dist)
		} else if(as.vector(strand(gene[gi])) == "-") {
			gene_tss_dist = ifelse(start(gr) > end(tss), -gene_tss_dist, gene_tss_dist)
		}
		gr$gene_tss_dist = gene_tss_dist

		## distance to tx tss
		## distance to tx tss
		tx = tx_list[[gi]]
		tx = tx[tx$tx_name != gi]
		tx_width = width(tx)
		tss = promoters(tx, upstream = 1, downstream = 0)
		dist = as.data.frame(distanceToNearest(gr, tss, select = "all"))
		dist = data.frame(queryHits = as.vector(tapply(dist[[1]], dist[[1]], unique)),
			              subjectHits = as.vector(tapply(dist[[2]], dist[[1]], function(x) x[which.max(tx_width[x])[1]])),
			              distance = as.vector(tapply(dist[[3]], dist[[1]], unique)))
		tx_tss_dist = dist[, 3]
		if(as.vector(strand(gene[gi])) == "+") {
			tx_tss_dist = ifelse(end(gr) < start(tss[dist[,2]]), -tx_tss_dist, tx_tss_dist)
		} else if(as.vector(strand(gene[gi])) == "-") {
			tx_tss_dist = ifelse(start(gr) > end(tss[dist[,2]]), -tx_tss_dist, tx_tss_dist)
		}
		gr$tx_tss_dist = tx_tss_dist
		gr$nearest_tx_tss = tss[dist[,2]]$tx_name

		gr
	}, mc.cores = mc.cores)

	# add mean methylation and meth_diameter, meth_IQR columns
	sample_id = cr_param$sample_id
	subgroup = cr_param$subgroup
	cov_filter = cr_param$cov_filter
	cov_cutoff = cr_param$cov_cutoff
	raw_meth = cr_param$raw_meth
	min_dp = cr_param$min_dp
	cor_method = cr_param$cor_method

	prev_chr = ""

	for(i in seq_along(cr_list)) {

		gi = names(cr_list)[i]
		message(qq("re-calculating mean methylation for reduced cr for @{gi}"))

		cr = cr_list[[i]]
		cr = c(reduce(cr[cr$corr > 0]), reduce(cr[cr$corr < 0]))
		cr_reduced = cr_reduced_list[[i]]

		mtch_reduced = as.matrix(findOverlaps(cr, cr_reduced))

		chr = as.vector(seqnames(cr))[1]

		if(chr != prev_chr) {
			methylation_hooks$set_chr(chr, verbose = FALSE)
			meth_gr = methylation_hooks$gr
			if(raw_meth) {
				meth = methylation_hooks$raw[, sample_id]
			} else {
				meth = methylation_hooks$meth[, sample_id]
			}
			cov = methylation_hooks$cov[, sample_id]

			meth[cov < cov_cutoff] = NA

			if(!is.null(cov_filter)) {
				
				l = apply(cov, 1, cov_filter)
				if(any(is.na(l))) {
					stop("`cov_filter` generates `NA`, check it.")
				}
				meth_gr = meth_gr[l]
				meth = meth[l, , drop = FALSE]
				cov = cov[l, , drop = FALSE]
			}
			prev_chr = chr
		}

		# mean methylation for each reduced cr
		mtch = as.matrix(findOverlaps(cr, meth_gr))
		mtch_cr_list = split(mtch[, 2], mtch[, 1])

		m = tapply(mtch_reduced[, 1], mtch_reduced[, 2], function(ind) {
			ind = unlist(mtch_cr_list[as.character(ind)])
			colMeans(meth[ind, , drop = FALSE], na.rm = TRUE)
		})
		m = do.call("rbind", m)

		colnames(m) = paste0("mean_meth_", colnames(meth))
		
		if(!is.null(expr)) {
			e = expr[gi, sample_id]
			corr = apply(m, 1, function(x) {
				l = !is.na(x)
				if(sum(l) < min_dp) return(NA)
				cor(x[l], e[l], method = cor_method)
			})
			corr_p = suppressWarnings(apply(m, 1, function(x) {
				l = !is.na(x)
				if(sum(l) < min_dp) return(NA)
				cor.test(x[l], e[l], method = cor_method)$p.value
			}))
		}
		if(!is.null(subgroup)) {
			if(length(unique(subgroup)) == 1) subgroup = NULL
		}

		if(!is.null(subgroup)) {
			subgroup = as.vector(subgroup)
			meth_anova = apply(m, 1, function(x) {
				l = !is.na(x)
				data = data.frame(value = x[l], class = subgroup[l], stringsAsFactors = FALSE)
				if(length(unique(data$class)) < 2) return(NA)
				if(any(table(data$class) < 2)) return(NA)
				oneway.test(value ~ class, data = data)$p.value
			})
			meth_diameter = apply(m, 1, function(x) {
				l = !is.na(x)
				if(any(table(subgroup[l]) < 2)) return(NA)
				diameter(as.vector(tapply(x[l], subgroup[l], mean)))
			})
			if(length(unique(subgroup)) == 2) {
				unique_subgroup = sort(unique(subgroup))
				meth_diff = apply(m, 1, function(x) {
					l = !is.na(x)
					tb = table(subgroup[l])
					if(any(tb < 2)) return(NA)
					y = as.vector(tapply(x[l], subgroup[l], mean))
					names(y) = unique(subgroup[l])
					y[unique_subgroup[1]] - y[unique_subgroup[2]]
				})
			}
		}
		if(nrow(m) == 1) {
			meth_IQR = IQR(m, na.rm = TRUE)
		} else {
			meth_IQR = rowIQRs(m, na.rm = TRUE)
		}
		
		if(is.null(expr)) {
			if(is.null(subgroup)) {
				df = DataFrame(m, # mean methylation
				    meth_IQR = meth_IQR)	
			} else {
				df = DataFrame(m, # mean methylation
				    meth_IQR = meth_IQR,
				    meth_anova = meth_anova,
				    meth_diameter = meth_diameter,
				    meth_diff = meth_diff)
			}
		} else {
			if(is.null(subgroup)) {
				df = DataFrame(m, # mean methylation
					corr = corr,
		   			corr_p = corr_p,
				    meth_IQR = meth_IQR)	
			} else {
				df = DataFrame(m, # mean methylation
					corr = corr,
		    		corr_p = corr_p,
				    meth_IQR = meth_IQR,
				    meth_anova = meth_anova,
				    meth_diameter = meth_diameter,
				    meth_diff = meth_diff)
			}
		}
		mcols(cr_reduced_list[[i]]) = cbind(mcols(cr_reduced_list[[i]]), df)
	}

	cr2 = do.call("c", cr_reduced_list)

	metadata(cr2)$cr_param = cr_param
	return(cr2) 
}


# == title
# Add subtype specificity columns in cr
#
# == param
# -cr correlated regions
# -cutoff cutoff for p-values of ANOVA test
# -suffix suffix of column names
#
# == details
# If ``subgroup`` is set in `correlated_regions`, this function can assign subtype specificity to each subtype.
#
# We use following digits to represent subtype specificity: 1 is defined as the methylation is higher than all other subtypes and the difference is significant.
# -1 is defined as the methylation is lower than all other subtypes and the difference is significant.
# All the others are defined as 0.
#
# == value
# A `GenomicRanges::GRanges` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
cr_add_subtype_specificity = function(cr, cutoff = 0.05, suffix = "_ss") {
	
	cr_param = metadata(cr)$cr_param
	subgroup = cr_param$subgroup
	subgroup_level = unique(subgroup)
	n_subgroup = length(subgroup_level)
	sample_id = cr_param$sample_id

	if(n_subgroup <= 1) {
		warning("Number of subgroups < 2.")
		return(cr)
	}

	subtype_ss = matrix(nrow = length(cr), ncol = n_subgroup)
	colnames(subtype_ss) = subgroup_level

	meth_mat = mcols(cr)
	meth_mat = meth_mat[, grep("^mean_meth_", colnames(meth_mat))]
	meth_mat = as.matrix(meth_mat)
	counter = set_counter(length(cr))
	for(i in seq_len(length(cr))) {
		x = meth_mat[i, paste0("mean_meth_", sample_id)]
		# pairwise t-test, t-value and p-value
		mat_t = matrix(nrow = n_subgroup, ncol = n_subgroup)
		rownames(mat_t) = subgroup_level
		colnames(mat_t) = subgroup_level
		mat_p = mat_t

		for(i1 in 2:n_subgroup) {
			for(i2 in 1:(i1-1)) {
				x1 = x[subgroup == subgroup_level[i1]]
				x2 = x[subgroup == subgroup_level[i2]]
				test = t.test(x1, x2)
				mat_t[i1, i2] = test$statistic
				mat_p[i1, i2] = test$p.value
				mat_t[i2, i1] = -mat_t[i1, i2]
				mat_p[i2, i1] = mat_p[i1, i2]	
			}
		}

		ss = apply(mat_t, 1, function(x) {
			x = x[!is.na(x)]
			if(all(x > 0)){
				return(1)
			} else if(all(x < 0)) {
				return(-1)
			} else {
				return(0)
			}
		})

		l = apply(mat_p, 1, function(x) {
			x = x[!is.na(x)]
			all(x < cutoff)
		})
		ss[!l] = 0
		subtype_ss[i, ] = ss

		counter()
	}

	colnames(subtype_ss) = paste0(colnames(subtype_ss), suffix)

	mcols(cr) = cbind(as.data.frame(mcols(cr)), subtype_ss)
	return(cr)
}
