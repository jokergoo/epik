
# == title
# Normalize histone modification signals to target
#
# == param
# -target target regions
# -mark histone mark name
# -sample_id a vector of sample ids
# -mode how to summarize histone modification signals among samples, by defualt is the cross-sample mean signal
# -return_arr whether also return the three dimension array itself
# -... pass to `EnrichedHeatmap::normalizeToMatrix`
#
# == details
# For each sample, the signal is normalized as a matrix, which results an array in which the third
# dimension corresponds to samples. The final normalized matrix which shows e.g. mean signal matrix is calcualted
# by ``apply(array, c(1, 2), mode)``.
#
# == values
# If ``return_arr`` is set to ``FALSE``, the funtion returns a matrix which can be directly sent to 
# `EnrichedHeatmap::EnrichedHeatmap`. If ``return_arr`` is ``TRUE``, the returned value is a list in which
# the first element is the original array that each slice in the third dimension is the normalize matrix
# in each sample.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
enrich_with_histone_mark = function(target, mark, sample_id, mode = mean, return_arr = FALSE, ...) {

	target_name = deparse(substitute(target))
	# 200 is number of windows (5000 + 5000)/50

	all_target_chromosome = unique(as.vector(seqnames(target)))
	hm_list = get_peak_list(mark, sample_id = intersect(sample_id, chipseq_hooks$sample_id(mark)), chr = all_target_chromosome)

	sample = names(hm_list)
	flag = 0
	for(i in seq_along(sample)) {
		message(qq("@{sample[i]}: normalize @{mark} signal to @{target_name}."))
	    tm = normalizeToMatrix(hm_list[[i]], target, trim = c(0, 0.01), ...)
	    if(!flag) {
	    	arr = array(dim = c(length(target), dim(tm)[2], length(hm_list)))
	    	dimnames(arr) = list(rownames(tm), colnames(tm), names(hm_list))
	    	flag = 1
	    }
	    arr[, , i] = tm
	}

	mat = apply(arr, c(1, 2), mode, na.rm = TRUE)
	mat = copyAttr(tm, mat)

	if(return_arr) {
		return(list(arr = arr, list = mat))
	} else {
		return(mat)
	}
}

# == title
# Normalize methylation to target
#
# == param
# -target target regions
# -sample_id a vector of sample ids
# -mode how to summarize methylation among samples, by default it is the cross-sample mean methylation. 
#       Since methylation is represented as matrix, here we use ``row*``-family functions (e.g. `rowMeans`, `matrixStats::rowMedians`)
# -extend pass to `EnrichedHeatmap::normalizeToMatrix`
# -smooth pass to `EnrichedHeatmap::normalizeToMatrix`
# -... pass to `EnrichedHeatmap::normalizeToMatrix`
#
# == value
# A matrix which can be directly visualized by `EnrichedHeatmap::EnrichedHeatmap`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
enrich_with_methylation = function(target, sample_id, mode = rowMeans, extend = 5000, smooth = TRUE, ...) {

	target_name = deparse(substitute(target))

	if(length(extend) == 1) extend = rep(extend, 2)

	# extrace methylation value which in [-5k, 5k] from TSS, and calculate mean methylation in each subgroup
	target_extend = GRanges(seqnames = seqnames(target), ranges = IRanges(start(target)-extend[1], end(target)+extend[2]))
	# process raw methylation data
	meth_gr = GRanges()
	for(chr in sort(unique(as.vector(seqnames(target))))) {
	    methylation_hooks$set_chr(chr, verbose = FALSE)
	    gr = methylation_hooks$gr
	    mtch = as.matrix(findOverlaps(gr, target_extend))
	    ind = unique(mtch[, 1])
	    mm = mode(methylation_hooks$meth[ind, sample_id, drop = FALSE], na.rm = TRUE)
	    gr = gr[ind]
	    mcols(gr) = mm
	    meth_gr = c(meth_gr, gr)
	}
	colnames(mcols(meth_gr)) = "mean_meth"
	message(qq("normalize methylation signals to @{target_name}..."))
	mat = normalizeToMatrix(meth_gr, target, value_column = "mean_meth", extend = extend, smooth = smooth, ...)
	return(mat)
}


dist_by_closeness2 = function(...) {
	as.dist(dist_by_closeness(...))
}

