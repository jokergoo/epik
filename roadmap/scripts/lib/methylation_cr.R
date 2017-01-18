diameter = function (x) {
    r = range(x)
    r[2] - r[1]
}

correlated_regions_per_gene = function(site, meth, cov, expr, chr, cov_cutoff = 3, min_dp = 4,
	cor_method = "spearman", window_size = 5, window_step = window_size, 
	factor = NULL, max_width = 100000) {

	if(ncol(meth) != length(expr)) {
		stop("number of columsn of `meth` should be same as length of `expr`.\n")
	}

	# index = seq(1, length(site), by = window_size)
	
	# i = seq_len(length(index) - 1)
	# ir = IRanges(site[index[i]], site[index[i+1]-1])

	# index for the starting site
	index = seq(1, length(site) - window_size, by = window_step)  # just forget last few sites
	
	i = seq_len(length(index))
	ir = IRanges(site[index[i]], site[index[i]+window_size-1])
	l = width(ir) <= max_width
	i = i[l]
	ir = ir[l]

	m = lapply(index[i], function(x) {
		ind = x+0:(window_size-1)
		meth_m = meth[ind, , drop = FALSE]
		cov_m = cov[ind, , drop = FALSE]
		sapply(seq_len(ncol(meth_m)), function(i) {
			xm = meth_m[, i]
			ym = cov_m[, i]
			mean(xm[ym >= cov_cutoff], na.rm = TRUE)
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

	if(!is.null(factor)) {
		if(length(unique(factor)) == 1) factor = NULL
	}

	if(!is.null(factor)) {
		factor = as.vector(factor)
		meth_anova = apply(m, 1, function(x) {
			l = !is.na(x)
			data = data.frame(value = x[l], class = factor[l], stringsAsFactors = FALSE)
			if(length(unique(data$class)) < 2) return(NA)
			if(any(table(data$class) < 2)) return(NA)
			oneway.test(value ~ class, data = data)$p.value
		})
		meth_diameter = apply(m, 1, function(x) {
			l = !is.na(x)
			if(any(table(factor[l]) < 2)) return(NA)
			diameter(as.vector(tapply(x[l], factor[l], mean)))
		})
	}
	if(nrow(m) == 1) {
		meth_IQR = IQR(m, na.rm = TRUE)
	} else {
		meth_IQR = rowIQRs(m, na.rm = TRUE)
	}
	gr = GRanges(seqnames = rep(chr, length(ir)),
		    ranges = ir)
	if(is.null(factor)) {
		df = DataFrame(n = window_size,
		    m, # mean methylation
		    corr = corr,
		    corr_p = corr_p,
		    meth_IQR = meth_IQR)	
	} else {
		df = DataFrame(n = window_size,
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

correlated_regions = function(sample_id, expr, txdb, chr, extend = 50000,
	cov_filter = function(x) sum(x > 0, na.rm = TRUE) > length(x)/2,
	cor_method = "spearman", factor = NULL, window_size = 5, window_step = window_size, max_width = 100000,
	raw_meth = FALSE, cov_cutoff = 3, min_dp = 4, col = NULL) {

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

	expr = expr[, sample_id, drop = FALSE]

	all_gi = rownames(expr)
	n_gene = length(all_gi)

	methylation_hooks$set(chr)
	site = methylation_hooks$site()
	if(raw_meth) {
		meth = methylation_hooks$raw(col_index = sample_id)
	} else {
		meth = methylation_hooks$meth(col_index = sample_id)
	}
	cov = methylation_hooks$coverage(col_index = sample_id)

	if(!is.null(cov_filter)) {
		
		l = apply(cov, 1, cov_filter)
		if(any(is.na(l))) {
			stop("`cov_filter` generates `NA`, check it.")
		}
		site = site[l]
		meth = meth[l, , drop = FALSE]
		cov = cov[l, , drop = FALSE]
	}

	if(!raw_meth) cov_cutoff = 0
	
	res = GRanges()
	for(i in seq_len(n_gene)) {
			
		# current gene name
		gi = all_gi[i]

		# expression of current gene
		e = expr[gi, ]

		qq.options(cat_prefix = qq("[@{chr}:@{gi}, @{i}/@{n_gene}]"))

		# if gene has low expression in many samples
		if(all(e == 0) || sd(e) == 0)  {
			qqcat("@{gi} has zero expression in all samples, skip\n")
			next
		}

		start = start(genemodel[gi])
		end = end(genemodel[gi])
		gm_site_index = extract_sites(start, end, site, TRUE, 0)
		gm_site = site[gm_site_index]
		gm_meth = meth[gm_site_index, sample_id, drop = FALSE]
		gm_cov = cov[gm_site_index, sample_id, drop = FALSE]

		if(length(gm_site) < 10) {
			qqcat("@{gi} has too few cpg sites, skip\n")
			next
		}

		qqcat("...\n")
		gr = correlated_regions_per_gene(gm_site, gm_meth, gm_cov, e, cov_cutoff = cov_cutoff, chr = chr,
			factor = factor, cor_method = cor_method, window_size = window_size, window_step = window_step, min_dp = min_dp,
			max_width = max_width)
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
		dist = data.frame(queryHits = tapply(dist[[1]], dist[[1]], unique),
			              subjectHits = tapply(dist[[2]], dist[[1]], function(x) x[which.max(tx_width[x])[1]]),
			              distance = tapply(dist[[3]], dist[[1]], unique))
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
	
	qq.options(cat_prefix = NULL)

	attr(res, "factor") = factor
	attr(res, "col") = col
	attr(res, "cor_method") = cor_method
	attr(res, "extend") = extend
	attr(res, "window_size") = window_size
	attr(res, "window_step") = window_step
	attr(res, "max_width") = max_width
	attr(res, "sample_id") = sample_id
	attr(res, "cor_method") = cor_method
	attr(res, "cov_filter") = cov_filter
	attr(res, "raw_meth") = raw_meth
	attr(res, "cov_filter") = cov_filter
	attr(res, "min_dp") = min_dp

	return(res)
}

readRDS_or_readRData = function(file) {
	if(grepl("\\.rds$", file, ignore.case = TRUE)) {
		cr = readRDS(file)
	} else if(grepl("\\.rdata", file, ignore.case = TRUE)) {
		var_name = load(file)
		eval(parse(text = paste0("cr = ", var_name)))
	}
	cr
}

# == title
# Get correlated regions with significant correlations
# 
# == param
# -chromosome a vector of chromosome names
# -template template path to find cr files
# -cutoff cutoff of adjusted correlation p-values
# -adj_method method for calculating adjusted p-values
# -meth_diameter_cutoff cutoff for methylation diameters
# -meth_IQR_cutoff cutoff for IQR, if there is no subtype information, IQR is used to remove less variable methylation
# -anova_cutoff cutoff for adjust ANOVA p-values
#
# == details
# As explained in `correlated_regions`, original cr is huge and is always saved as a separated file for single chromosome.
# Here ``template`` defined how to get the cr files. E.g. if ``template``is defined as "path_to/@{chr}_cr.rds", the funciton
# will replace "@{chr}" to every chromosome and read the data.
#
# Two additional columns are attached:
#
# -corr_fdr FDR for correlation p-values
# -meth_anova_fdr FDR for anova test, will be added only if ``factor`` is set in `correlated_regions`.
#
# == value
# A `GenomicRanges::GRanges` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
filter_correlated_regions = function(chromosome = paste0("chr", 1:22), template, 
	cutoff = 0.05, adj_method = "BH", meth_diameter_cutoff = 0.25, meth_IQR_cutoff = 0.25,
	anova_cutoff = 0.05, type = "all") {

	if(length(cutoff) == 1) cutoff = rep(cutoff, 2)

	cat("calculate fdr...\n")
	corr_p = NULL
	meth_anova = NULL
	meth_diameter = NULL
	meth_IQR = NULL
	chr_name = NULL
	corr = NULL
	for(chr in chromosome) {
		qqcat("reading cr for @{chr}\n")
		cr = readRDS_or_readRData(qq(template))
		if(type == "neg") {
			l = cr$corr < 0
			l[is.na(l)] = FALSE
			cr = cr[l]
		} else if(type == "pos") {
			l = cr$corr > 0
			l[is.na(l)] = FALSE
			cr = cr[l]
		}

		has_anova = FALSE
		if("meth_anova" %in% colnames(mcols(cr))) {
			has_anova = TRUE
		}

		if(has_anova) {
			cr = cr[(!is.na(cr$corr)) & (!is.na(cr$corr_p)) & (!is.na(cr$meth_anova)) & (!is.na(cr$meth_diameter))]
		} else {
			cr = cr[(!is.na(cr$corr)) & (!is.na(cr$corr_p))]
		}
	    
	    corr_p = c(corr_p, cr$corr_p)
	    corr = c(corr, cr$corr)
	   	if(has_anova) {
	   		meth_anova = c(meth_anova, cr$meth_anova)
	   		meth_diameter = c(meth_diameter, cr$meth_diameter)
	   	} else {
	   		meth_mat = as.matrix(mcols(cr)[, grep("^mean_meth", colnames(mcols(cr)))])
	   		meth_IQR = c(meth_IQR, rowIQRs(meth_mat, na.rm = TRUE))
	   		# meth_IQR = c(meth_IQR, cr$meth_IQR)
	   	}
	    chr_name = c(chr_name, rep(chr, length(cr)))
	}

	if(has_anova) {
		anova_fdr = p.adjust(meth_anova, method = adj_method)
		l = anova_fdr <= cutoff[1] & meth_diameter >= meth_diameter_cutoff
		qqcat("filter out @{sum(!l)}/@{length(l)} by differential methylation.\n")
		corr_fdr = rep(Inf, length(corr_p))
		corr_fdr[l] = p.adjust(corr_p[l], method = adj_method)
		l = l & ifelse(corr > 0, corr_fdr <= cutoff[1], corr_fdr <= cutoff[2])
	} else {
		corr_fdr = p.adjust(corr_p, method = adj_method)
		l = ifelse(corr > 0, corr_fdr <= cutoff[1], corr_fdr <= cutoff[2]) & meth_IQR >= meth_IQR_cutoff & !is.na(meth_IQR)
	}

	l = l & !is.na(corr)

	cat("filter by fdr...\n")
	cr2 = GRanges()
	for(chr in chromosome) {
		qqcat("reading cr for @{chr}\n")
		cr = readRDS_or_readRData(qq(template))
		if(type == "neg") {
			l = cr$corr < 0
			l[is.na(l)] = FALSE
			cr = cr[l]
		} else if(type == "pos") {
			l = cr$corr > 0
			l[is.na(l)] = FALSE
			cr = cr[l]
		}

		if(has_anova) {
			cr = cr[(!is.na(cr$corr)) & (!is.na(cr$corr_p)) & (!is.na(cr$meth_anova)) & (!is.na(cr$meth_diameter))]
		} else {
			cr = cr[(!is.na(cr$corr)) & (!is.na(cr$corr_p))]
		}

		lo = chr_name == chr
		cr$corr_fdr = corr_fdr[lo]
		if(has_anova) {
			cr$meth_anova_fdr = anova_fdr[lo]
		}
		if(sum(l[lo])) {
			cr2 = suppressWarnings(c(cr2, cr[l[lo]]))
		}
	}

	attr(cr2, "factor") = attr(cr, "factor")
	attr(cr2, "col") = attr(cr, "col")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "extend") = attr(cr, "extend")
	attr(cr2, "window_size") = attr(cr, "window_size")
	attr(cr2, "window_step") = attr(cr, "window_step")
	attr(cr2, "max_width") = attr(cr, "max_width")
	attr(cr2, "sample_id") = attr(cr, "sample_id")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "cov_filter") = attr(cr, "cov_filter")
	attr(cr2, "cov_cutoff") = attr(cr, "cov_cutoff")
	attr(cr2, "raw_meth") = attr(cr, "raw_meth")
	attr(cr2, "min_dp") = attr(cr, "min_dp")

	cr2
}

# gr is sorted
reduce_cr_by_gene = function(gr, gap = 1, window_size, window_step) {
	if(length(gr) == 0) return(GRanges())

	n = length(gr)
	if(all(start(gr)[-1] > end(gr)[-n]+1)) {
		gr2 = gr
		mcols(gr2) = NULL
		gr2$mean_corr = gr$corr
		gr2$merged_windows = rep(1, n)
		gr2$n = window_size
	} else {
		gr2 = reduce(gr, min.gap = gap, with.revmap = TRUE)
		revmap = mcols(gr2)$revmap
		mcols(gr2) = NULL
		gr2$mean_corr = sapply(revmap, function(ind) mean(gr[ind]$corr, na.rm = TRUE))
		gr2$merged_windows = sapply(revmap, length)
		gr2$n = gr2$merged_windows*window_size - (gr2$merged_windows-1)*(window_size - window_step)
	}
	qqcat("  reduced from @{length(gr)} to @{length(gr2)}\n")
	return(gr2)
}

reduce_cr = function(cr, txdb, gap = 1, mc.cores = 1, window_size = 6, window_step = 3) {

	sample_id = attr(cr, "sample_id")
	
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
		gr = c(reduce_cr_by_gene(neg_cr, gap = gap, window_size = window_size, window_step = window_step),
			   reduce_cr_by_gene(pos_cr, gap = gap, window_size = window_size, window_step = window_step))
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
		tx = tx_list[[gi]]
		tx = tx[tx$tx_name != gi]
		tss = promoters(tx, upstream = 1, downstream = 0)
		dist = as.data.frame(distanceToNearest(gr, tss))
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

	cr2 = copy_cr_attribute(cr, cr2)

	cr2
}

copy_cr_attribute = function(cr, cr2) {
	attr(cr2, "factor") = attr(cr, "factor")
	attr(cr2, "col") = attr(cr, "col")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "extend") = attr(cr, "extend")
	attr(cr2, "window_size") = attr(cr, "window_size")
	attr(cr2, "window_step") = attr(cr, "window_step")
	attr(cr2, "max_width") = attr(cr, "max_width")
	attr(cr2, "sample_id") = attr(cr, "sample_id")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "cov_filter") = attr(cr, "cov_filter")
	attr(cr2, "cov_cutoff") = attr(cr, "cov_cutoff")
	attr(cr2, "raw_meth") = attr(cr, "raw_meth")
	attr(cr2, "min_dp") = attr(cr, "min_dp")
	cr2
}

cr_qc = function(neg_cr, pos_cr, txdb, chromosome = paste0("chr", 1:22), region = c("all", "intergenic", "tss", "gene")) {

	gm = genes(txdb)
	region = match.arg(region)[1]

	n_neg = matrix(0, nrow = 8, ncol = 5)
	IQR_cutoff = c(0, 0.1, 0.2, 0.3, 0.4)
	colnames(n_neg) = paste0("IQR > ", IQR_cutoff)
	cor_cutoff = 0:7/10
	rownames(n_neg) = paste(">", cor_cutoff)
	n_pos = n_neg
	w_neg = n_neg
	w_pos = n_neg
	for(i in seq_along(cor_cutoff)) {
		for(j in seq_along(IQR_cutoff)) {
			qqcat("cor_cutoff = @{cor_cutoff[i]}, IQR_cutoff = @{IQR_cutoff[j]}\n")

			neg_cr2 = neg_cr[abs(neg_cr$corr) > cor_cutoff[i] & neg_cr$meth_IQR > IQR_cutoff[j]]
			pos_cr2 = pos_cr[abs(pos_cr$corr) > cor_cutoff[i] & pos_cr$meth_IQR > IQR_cutoff[j]]

			neg_cr2_list = split(neg_cr2, neg_cr2$gene_id)
			gi = unique(neg_cr2$gene_id)
			gl = width(gm[gi])
			names(gl) = gi
			for(k in seq_along(neg_cr2_list)) {
				neg_cr3 = neg_cr2_list[[k]]
				if(region == "tss") {
					neg_cr3 = neg_cr3[neg_cr3$gene_tss_dist > -1000 & neg_cr3$gene_tss_dist < 2000]
				} else if(region == "gene") {
					neg_cr3 = neg_cr3[neg_cr3$gene_tss_dist > 2000 & neg_cr3$gene_tss_dist < gl[gi[k]]]
				} else if(region == "intergenic") {
					neg_cr3 = neg_cr3[neg_cr3$gene_tss_dist < -1000 | neg_cr3$gene_tss_dist > gl[gi[k]]]
				}
				# neg_cr3 = reduce_cr_by_gene(neg_cr3)
				n_neg[i, j] = n_neg[i, j] + length(neg_cr3)
				w_neg[i, j] = w_neg[i, j] + sum(width(neg_cr3))
				qqcat("  @{k}/@{length(neg_cr2_list)} neg_cr gene...\n")
			}

			pos_cr2_list = split(pos_cr2, pos_cr2$gene_id)
			gi = unique(pos_cr2$gene_id)
			for(k in seq_along(pos_cr2_list)) {
				pos_cr3 = pos_cr2_list[[k]]
				if(region == "tss") {
					pos_cr3 = pos_cr3[pos_cr3$gene_tss_dist > -1000 & pos_cr3$gene_tss_dist < 2000]
				} else if(region == "gene") {
					pos_cr3 = pos_cr3[pos_cr3$gene_tss_dist > 2000 & pos_cr3$gene_tss_dist < gl[gi[k]]]
				} else if(region == "intergenic") {
					pos_cr3 = pos_cr3[pos_cr3$gene_tss_dist < -1000 | pos_cr3$gene_tss_dist > gl[gi[k]]]
				}
				# pos_cr3 = reduce_cr_by_gene(pos_cr3)
				n_pos[i, j] = n_pos[i, j] + length(pos_cr3)
				w_pos[i, j] = w_pos[i, j] + sum(width(pos_cr3))
				qqcat("  @{k}/@{length(pos_cr2_list)} pos_cr gene...\n")
			}
		}
	}

	r1 = n_neg/n_pos
	r2 = w_neg/w_pos
	par(mfrow = c(1, 2))
	matplot(r1, xlab = "abs_correlation", type = "b", pch = 1, ylab = "#neg/#pos", col = 1:5, axes = FALSE, main = qq("region = @{region}"))
	legend("topleft", lty = 1, col = 1:5, legend = colnames(r1))
	axis(side = 1, at = 1:8, labels = rownames(r1), cex = 0.8)
	axis(side = 2)
	box()
	matplot(r2, xlab = "abs_correlation", type = "b", pch = 1, ylab = "width(neg)/width(pos)", col = 1:5, axes = FALSE, main = qq("region = @{region}"))
	axis(side = 1, at = 1:8, labels = rownames(r1), cex = 0.8)
	axis(side = 2)
	box()
}



cr_overlap_to_genomic_features = function(gf_list, template, cutoff = 0, species = NULL, 
	chromosome = paste0("chr", 1:22)) {

	cr = GRanges()
	for(chr in chromosome) {
		qqcat("loading @{chr}...\n")
		cr_tmp = readRDS_or_readRData(qq(template))
		l = !is.na(cr_tmp$corr)
		cr = c(cr, cr_tmp[l])
	}

	jaccard = matrix(0, ncol = length(gf_list), nrow = length(cutoff))
	colnames(jaccard) = names(gf_list)
	rownames(jaccard) = paste(ifelse(cutoff > 0, ">", ifelse(cutoff < 0, "<", "=")), cutoff)

	gf_list = lapply(gf_list, function(gr) {
		gr = reduce(gr)
	    gr[seqnames(gr) %in% chromosome]
	})
	
	
	for(j in seq_along(cutoff)) {
		if(cutoff[j] > 0) {
	    	cr2 = cr[cr$corr > cutoff[j]]
	    } else if(cutoff[j] < 0) {
	    	cr2 = cr[cr$corr < cutoff[j]]
	    } else {
	    	cr2 = cr
	    }
	    cr_reduced = reduce(cr2)

		for(i in seq_along(gf_list)) {
		    gr = gf_list[[i]]
		    mtch = as.matrix(findOverlaps(cr_reduced, gr))
		    qqcat("[@{names(gf_list)[i]}] [@{cutoff[j]}] @{length(mtch)} matches\n")
		    jaccard[j, i] = sum(as.numeric(width(pintersect(cr_reduced[mtch[,1]], gr[mtch[,2]]))))/sum(as.numeric(width(punion(cr_reduced[mtch[,1]], gr[mtch[,2]]))))
		}
	}
	return(jaccard)
}

