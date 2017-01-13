
# title
# correlated regions in a specified regio
#
# == param
# -site position of CpG sites, should be sorted
# -meth methylation matrix corresponding to ``site``
# -expr expression for the associated gene
# -chr chromosome, used to construct the `GenomicRanges::GRanges` object
# -cov CpG coverage matrix
# -cov_cutoff cutoff for coverage
# -min_dp minimal number of non-NA values for calculating correlations
# -cor_method method for calcualting correlations
# -window_size how many CpG sites in a window
# -subgroup subgroup information
# -max_width maximum width of a window
#
# == details
# ``cov`` and ``cov_cutoff`` should be set when the methylation is unsmoothed, because
# for the unsmoothed data, the methylation rate is not unreliable when the CpG coverage is low.
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
		    meth_diameter = meth_diameter)
	}
	mcols(gr) = df

	return(gr)
}

# == title
# Correlation between methylation and expression
#
# == param
# -sample_id a vector of sample id
# -expr expression matrix in which columns correspond to sample ids
# -txdb a ``GenomicFeatures::GRanges`` object. Gene names should be same type as row names in ``expr``
# -chr a single chromosome
# -extend extension of gene model, both upstream and downstream
# -cov_filter if ``coverage`` hook is set in `methylation_hooks`, this option can be set to filter out CpG sites with low coverage across samples.
#          the value for this option is a function for which the argument is a vector of coverage values for current CpG in all samples.
# -cor_method method to calculate correlation
# -factor classes of samples
# -window_size number of CpGs in a window
# -max_width maximum width of a window
# -raw_meth whether use raw methylation value (values from ``raw`` hook set in `methylation_hooks`)
# -cov_cutoff cutoff for coverage
# -min_dp minimal non-NA values for calculating correlations
# -col color for classes
#
# == details
# The detection for correlated regions is gene-centric. For every gene, the process are as follows:
#
# - extend to both upstream and downstream
# - from the most upstream, use a sliding window which contains ``windows_size`` CpG sites
# - filter each window by CpG coverage (by ``cov_filter`` and ``cov_cutoff``)
# - calculate correlation between methylation and gene expression for this window
#
# Following meth columns are attached to the `GenomicRanges::GRanges` objects:
#
# -n number of CpG sites
# -mean_meth_* mean methylation in each window in every sample.
# -corr correlation
# -corr_p p-value for the correlation test
# -meth_IQR IQR of mean methylation if ``factor`` is not set
# -meth_anova p-value from oneway ANOVA test if ``factor`` is set
# -meth_diameter range between maximum mean and minimal mean in all subgroups if ``factor`` is set
# -gene_id gene id
# -gene_tss_dist distance to tss of genes
# -tx_tss_dist if genes have multiple transcripts, this is the distance to the nearest transcript
# -nearest_txx_tss transcript id of the nearest transcript
#
# This function keeps all the information for all CpG windows. Users can get `filter_correlated_regions` to get correlated regions
# with significant correlations and use `reduce_cr` to merge neighbouring windows.
#
# Since information for all CpG windows are kept, the size of the object is always very huge, thus, it is reasonable
# to analyze each chromosome separately and save each object as a single file. Some downstream functions expect a formatted
# path of the cr file.
#
# == value
# A `GenomicRanges::GRanges` object which contains associated statistics for every CpG windows.
#
# == seealso
# `filter_correlated_regions`, `reduce_cr`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
correlated_regions = function(sample_id, expr, txdb, chr, extend = 50000,
	cov_filter = function(x) sum(x > 0, na.rm = TRUE) > length(x)/2,
	cor_method = "spearman", subgroup = NULL, window_size = 5, window_step = window_size, 
	max_width = 10000, raw_meth = FALSE, cov_cutoff = 3, min_dp = 4, col = NULL) {

	qqcat("extracting gene model (extend = @{extend}, chr = @{chr})...\n")
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
	qqcat("@{length(sample_id)} samples are used\n")
	
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

	qqcat("on chromosome @{chr}, @{length(site)} CpG sites, @{n_gene} genes\n")

	if(!raw_meth) cov_cutoff = 0
	
	res = GRanges()
	for(i in seq_len(n_gene)) {
			
		# current gene name
		gi = all_gi[i]

		# expression of current gene
		e = expr[gi, ]

		# if gene has low expression in many samples
		if(all(e == 0) || sd(e) == 0)  {
			qqcat("[@{chr}:@{gi}, @{i}/@{n_gene}] @{gi} has zero expression in all samples, skip\n")
			next
		}

		start = start(genemodel[gi])
		end = end(genemodel[gi])
		gm_site_index = extract_sites(start, end, site, TRUE, 0)
		gm_site = site[gm_site_index]
		gm_meth = meth[gm_site_index, , drop = FALSE]
		gm_cov = cov[gm_site_index, , drop = FALSE]

		if(length(gm_site) < 10) {
			qqcat("[@{chr}:@{gi}, @{i}/@{n_gene}] @{gi} has too few cpg sites, skip\n")
			next
		}

		qqcat("[@{chr}:@{gi}, @{i}/@{n_gene}] @{length(gm_site)} CpG sites ...\n")
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

		res = c(res, gr)
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

	metadata(res) = list(cr_param = param)
	return(res)
}

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
		gr2$mean_corr = gr$corr
		gr2$merged_windows = rep(1, n)
		gr2$ncpg = window_size
	} else {
		gr2 = reduce(gr, min.gap = gap, with.revmap = TRUE)
		revmap = mcols(gr2)$revmap
		mcols(gr2) = NULL
		gr2$mean_corr = sapply(revmap, function(ind) mean(gr[ind]$corr, na.rm = TRUE))
		gr2$merged_windows = sapply(revmap, length)
		gr2$ncpg = gr2$merged_windows*window_size - (gr2$merged_windows-1)*(window_size - window_step)
	}
	qqcat("  @{deparse(substitute(gr))} has been reduced from @{length(gr)} to @{length(gr2)}\n")
	return(gr2)
}

reduce_cr = function(cr, txdb, gap = 1, mc.cores = 1) {

	cr_param = metadata(cr)$cr_param

	qqcat("extracting gene and tx models.\n")
	gene = genes(txdb)
	tx_list = transcriptsBy(txdb, by = "gene")

	cr_list = split(cr, cr$gene_id)
	i = 0
	res = mclapply(names(cr_list), function(gi) {

		qqcat("reducing cr on @{gi}, @{i <<- i+1}/@{length(cr_list)}...\n")
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

	cr2 = do.call("c", res)
	metadata(cr2)$cr_param = cr_param
	return(cr2) 
}


