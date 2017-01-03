
cr_enriched_at_tss = function(neg_cr, pos_cr, txdb, main = NULL) {
	
	gene = genes(txdb)
	tx = transcripts(txdb)
	tx_list = transcriptsBy(txdb, "gene")

	cr = c(neg_cr, pos_cr)

	gene = gene[names(gene) %in% cr$gene_id]
	tx = tx[tx$tx_name %in% cr$nearest_tx_tss]
	names(tx) = tx$tx_name
	tx$gene_id = structure(cr$gene_id, names = cr$nearest_tx_tss)[names(tx)]

	gene_tss = promoters(gene, upstream = 1, downstream = 0)
	tx_tss = promoters(tx, upstream = 1, downstream = 0)

	# align to gene tss
	mat_pos = normalizeToMatrix(pos_cr, gene_tss, mapping_column = "gene_id",
        extend = c(10000, 20000), mean_mode = "absolute", w = 50)
	mat_neg = normalizeToMatrix(neg_cr, gene_tss, mapping_column = "gene_id",
        extend = c(10000, 20000), mean_mode = "absolute", w = 50)
	x_pos = colSums(mat_pos)
	x_neg = colSums(mat_neg)


	# tx_list = tx_list[names(gene)]
	# tx_list = lapply(names(tx_list), function(gi) {
	# 	gr = tx_list[[gi]]
	# 	gr = gr[gr$tx_name != gi]
	# 	gr$gene_id = rep(gi, length(gr))
	# 	gr
	# })

	# tx2 = do.call("c", tx_list)

	# mat_tx = normalizeToMatrix(tx2, gene_tss, mapping_column = "gene_id",
	# 	extend = c(10000, 20000), mean_mode = "w0", w = 50)
	# mat_nearest_tx = normalizeToMatrix(tx, gene_tss, mapping_column = "gene_id",
	# 	extend = c(10000, 20000), mean_mode = "w0", w = 50)
	
	par(mar = c(4, 4, 4, 4))
	plot(x_pos, ylim = c(0, max(c(x_pos, x_neg))), type = "n", col = "red", axes = FALSE, xlab = 'pos relateive to gene tss', ylab = "coverage", main = main)
	if(length(pos_cr)) lines(x_pos, col = "red")
	if(length(neg_cr)) lines(x_neg, col = "green")
	abline(v = length(x_pos)/3, lty = 2, col = "grey")
	axis(side = 2)
	axis(side = 1, at = seq(1, length(x_pos), length = 7), labels = c(-10000, -5000, 0, 5000, 10000, 15000, 20000))

	legend("topright", lty = 1, col = c("green", "red"), legend = c("neg_cr", "pos_cr"))
	par(new = FALSE)

	# align to tx tss
	mat_pos = normalizeToMatrix(pos_cr, tx_tss, mapping_column = "nearest_tx_tss",
        extend = c(10000, 20000), mean_mode = "absolute", w = 50)
	mat_neg = normalizeToMatrix(neg_cr, tx_tss, mapping_column = "nearest_tx_tss",
        extend = c(10000, 20000), mean_mode = "absolute", w = 50)
	x_pos = colSums(mat_pos)
	x_neg = colSums(mat_neg)

	par(mar = c(4, 4, 4, 4))
	plot(x_pos, ylim = c(0, max(c(x_pos, x_neg))), type = "n", col = "red", axes = FALSE, xlab = 'pos relateive to nearest tx tss', ylab = "coverage", main = main)
	if(length(pos_cr)) lines(x_pos, col = "red")
	if(length(neg_cr)) lines(x_neg, col = "green")
	abline(v = length(x_pos)/3, lty = 2, col = "grey")
	axis(side = 2)
	axis(side = 1, at = seq(1, length(x_pos), length = 7), labels = c(-10000, -5000, 0, 5000, 10000, 15000, 20000))

	legend("topright", lty = 1, col = c("green", "red"), legend = c("neg_cr", "pos_cr"))
}


cr_enriched_at_gene_body = function(neg_cr, pos_cr, txdb, extend = attr(neg_cr, "extend"), main = NULL) {
	
	gene = genes(txdb)

	cr = c(neg_cr, pos_cr)
	
	gene = gene[names(gene) %in% cr$gene_id]
	if(length(extend) == 1) extend = c(extend, extend)
	
	# align to gene tss
	mat_pos = normalizeToMatrix(pos_cr, gene, mapping_column = "gene_id",
        extend = extend, mean_mode = "absolute", w = 50, target_ratio = 0.6)
	mat_neg = normalizeToMatrix(neg_cr, gene, mapping_column = "gene_id",
        extend = extend, mean_mode = "absolute", w = 50, target_ratio = 0.6)
	x_pos = colSums(mat_pos)
	x_neg = colSums(mat_neg)
	
	par(mar = c(4, 4, 4, 4))
	plot(x_pos, ylim = c(0, max(c(x_pos, x_neg))), type = "n", col = "red", axes = FALSE, xlab = 'pos', ylab = "coverage", main = main)
	if(length(pos_cr)) lines(x_pos, col = "red")
	if(length(neg_cr)) lines(x_neg, col = "green")
	upstream_index = attr(mat_pos, "upstream_index")
	target_index = attr(mat_pos, "target_index")
	abline(v = length(upstream_index) + 0.5, lty = 2, col = "grey")
	axis(side = 2)
	axis(side = 1, at = c(1, length(upstream_index), length(upstream_index) + length(target_index), length(x_neg)), 
		labels = c(-extend[1], "TSS", "TES", extend[2]))

	legend("topright", lty = 1, col = c("green", "red"), legend = c("neg_cr", "pos_cr"))
	par(new = FALSE)
}


# mean signal across samples
enrich_with_histone_mark = function(target, mark, sample_id, target_ratio = 0.1, return_arr = FALSE, 
	value_column = "density", extend = 5000, mean_mode = "w0", w = 50, ...) {

	target_name = deparse(substitute(target))
	# 200 is number of windows (5000 + 5000)/50

	hm_list = get_peak_list(mark, sample_id = intersect(sample_id, chipseq_hooks$sample_id(mark)))

	sample = names(hm_list)
	flag = 0
	for(i in seq_along(sample)) {
		qqcat("@{sample[i]}: normalize histone modifications to @{target_name}.\n")
	    tm = normalizeToMatrix(hm_list[[i]], target, target_ratio = target_ratio, 
	    	value_column = value_column, extend = extend, mean_mode = mean_mode, w = w, trim = c(0, 0.01), ...)
	    if(!flag) {
	    	arr = array(dim = c(length(target), dim(tm)[2], length(hm_list)))
	    	dimnames(arr) = list(rownames(tm), colnames(tm), names(hm_list))
	    	flag = 1
	    }
	    arr[, , i] = tm
	}

	mat = apply(arr, c(1, 2), mean, na.rm = TRUE)
	mat = copyAttr(tm, mat)

	if(return_arr) {
		return(list(arr = arr, list = mat))
	} else {
		return(mat)
	}
}

# mean methylation across samples
enrich_with_methylation = function(target, sample_id, target_ratio = 0.1, mode = rowMeans,
	extend = 5000, w = 50, mean_mode = "absolute", empty_value = NA, smooth = TRUE, ...) {

	target_name = deparse(substitute(target))

	# extrace methylation value which in [-5k, 5k] from TSS, and calculate mean methylation in each subgroup
	target_extend = GRanges(seqnames = seqnames(target), ranges = IRanges(start(target)-extend, end(target)+extend))
	# process raw methylation data
	meth_gr = GRanges()
	for(chr in sort(unique(as.vector(seqnames(target))))) {
	    methylation_hooks$set(chr)
	    gr = methylation_hooks$GRanges()
	    mtch = as.matrix(findOverlaps(gr, target_extend))
	    ind = unique(mtch[, 1])
	    mm = mode(methylation_hooks$meth(row_index = ind, col_index = sample_id), na.rm = TRUE)
	    gr = gr[ind]
	    mcols(gr) = mm
	    meth_gr = c(meth_gr, gr)
	}
	rm(gr)
	gc(verbose = FALSE)
	colnames(mcols(meth_gr)) = "mean_meth"
	qqcat("normalize methylation signals to @{target_name}...\n")
	mat = normalizeToMatrix(meth_gr, target, value_column = "mean_meth", target_ratio = target_ratio, extend = extend, 
		w = w, mean_mode = mean_mode, empty_value = empty_value, smooth = smooth, ...)
	return(mat)
}


set_counter = function(n, fmt = "%s") {

	n = as.integer(n)
	i = 1

	f = function() {
		if(interactive()) {
			pct = round(i/n*100, 1)
			str = paste0(i, "/", n, " (", pct, "%)")
			str = sprintf(fmt, str)

			cat(paste(rep("\b", 200), collapse=""))
			cat(str)
			if(i == n) cat("\n")

			i = i + 1
			assign("i", i, envir = parent.env(environment()))
			return(invisible(i))
		}
	}
}

dist_by_closeness2 = function(...) {
	as.dist(dist_by_closeness(...))
}


anno_enriched_by_sign = function(gp = gpar(col = "red"), pos_line = TRUE, pos_line_gp = gpar(lty = 2),
	yaxis = TRUE, ylim = NULL, value = c("mean", "sum", "abs_mean", "abs_sum"), yaxis_side = "right", 
	yaxis_gp = gpar(fontsize = 8), pos_col = "red", neg_col = "green") {

	# in case of lazy loading
	gp = gp
	pos_line = pos_line
	pos_line_gp = pos_line_gp
	yaxis = yaxis
	ylim = ylim
	yaxis_side = yaxis_side
	yaxis_gp = yaxis_gp

	value = match.arg(value)[1]
	function(index) {

		ht = get("object", envir = parent.frame(n = 5))
		mat = ht@matrix
		mat_pos = mat
		mat_pos[mat_pos < 0] = 0
		mat_neg = mat
		mat_neg[mat_neg > 0] = 0
		mat_pos = abs(mat_pos)
		mat_neg = abs(mat_neg)

		upstream_index = attr(mat, "upstream_index")
		downstream_index = attr(mat, "downstream_index")
		target_index = attr(mat, "target_index")

		n1 = length(upstream_index)
		n2 = length(target_index)
		n3 = length(downstream_index)
		n = n1 + n2 + n3
		
		y_pos = sapply(ht@row_order_list, function(i) {
			colMeans(mat_pos[i, , drop = FALSE], na.rm = TRUE)
		})
		y_neg = sapply(ht@row_order_list, function(i) {
			colMeans(mat_neg[i, , drop = FALSE], na.rm = TRUE)
		})

		if(is.null(ylim)) {
			ylim = range(c(y_neg, y_pos), na.rm = TRUE)
			ylim[2] = ylim[2] + (ylim[2] - ylim[1]) * 0.05
		}
		
		gp = EnrichedHeatmap:::recycle_gp(gp, length(ht@row_order_list))

		pushViewport(viewport(xscale = c(0, n), yscale = ylim))
		grid.rect(gp = gpar(col = "black", fill = NA))
		for(i in seq_len(ncol(y_pos))) {
			gp2 = EnrichedHeatmap:::subset_gp(gp, i); gp2$col = pos_col;class(gp2) = "gpar"
			grid.lines(seq_len(n)-0.5, y_pos[,i], default.units = "native", gp = gp2)
			gp2 = EnrichedHeatmap:::subset_gp(gp, i); gp2$col = neg_col;class(gp2) = "gpar"
			grid.lines(seq_len(n)-0.5, y_neg[,i], default.units = "native", gp = gp2)
		}
		if(pos_line) {
		    if(n1 && n2 && n3) {
                grid.lines(rep((n1-0.5)/n, 2), c(0, 1), gp = pos_line_gp)
                grid.lines(rep((n1+n2-0.5)/n, 2), c(0, 1), gp = pos_line_gp)
            } else if(n1 && !n2 && n3) {
                grid.lines(rep((n1-0.5)/n, 2), c(0, 1), gp = pos_line_gp)
            } else if(!n1 && n2 && n3) {
                grid.lines(rep((n1+n2-0.5)/n, 2), c(0, 1), gp = pos_line_gp)
            } else if(n1 && n2 && !n3) {
                grid.lines(rep((n1-0.5)/n, 2), c(0, 1), gp = pos_line_gp)
            }
		}
		if(yaxis) {
			if(yaxis_side == "right") {
				grid.yaxis(main = FALSE, gp = yaxis_gp)
			} else {
				grid.yaxis(gp = yaxis_gp)
			}
		}
	    upViewport()
	}
}

