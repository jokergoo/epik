### function for methylation analysis


# == title
# Heatmap for differential methylated genomic features
#
# == param
# -gr a `GenomicRanges::GRanges` object returned from `get_mean_methylation_in_genomic_features`
# -subgroup subgroup information
# -ha column annotations, a `ComplexHeatmap::HeatmapAnnotation-class` object
# -genomic_features a single or a list of `GenomicRanges::GRanges` objects that are used to annotate ``gr``
# -meth_diff minimal range between mean value in subgroups
# -cutoff if subgroup information is provided, cutoff of p-value of the oneway ANOVA test
# -adj_method method to calculate adjusted p-values
# -cluster_columns how to cluster columns. "subgroup" means only to cluster samples in every subgroup.
# -... pass to `ComplexHeatmap::Heatmap`
#
# == details
# Only those regions having differential methylation are visualized. 
#
# == value
# A `GenomicRanges::GRanges` object which only contains regions with significant differential methylation.
# 
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
heatmap_diff_methylation_in_genomic_features = function(gr, subgroup, 
	ha = HeatmapAnnotation(subgroup = subgroup, show_annotation_name = TRUE),
	genomic_features = NULL, 
	meth_diff = 0, cutoff = 0.05, adj_method = "BH", 
	cluster_columns = c("subgroup", "all", "none"), ...) {
		
	cluster_columns = match.arg(cluster_columns)[1]
	
	mat = mcols(gr)
	mat = as.matrix(mat[, grepl("mean_meth_", colnames(mat)), drop = FALSE])
	nr0 = nrow(mat)

	sample_id = gsub("mean_meth_", "", colnames(mat))
	colnames(mat) = sample_id

	l = gr$ncpg == 0
	if(sum(l) > 0) {
		message(qq("remove @{sum(l)}/@{length(l)} rows with ncpg == 0"))
		gr = gr[!l]
		mat = mat[!l, , drop = FALSE]
	}

	if(meth_diff > 0) {
		x = tapply(seq_len(ncol(mat)), subgroup, function(ind) {
				rowMeans(mat[, ind, drop = FALSE], na.rm = TRUE)
			})
		dim(x) = NULL
		x = as.matrix(data.frame(x))
		v = apply(x, 1, function(y) max(y, na.rm = TRUE) - min(y, na.rm = TRUE))
		l = v > meth_diff

		gr$meth_diff = v

		message(qq("remove @{sum(!l)}/@{length(l)} rows by filtering `meth_diff > @{meth_diff}`"))
		gr = gr[l]
		mat = mat[l, , drop = FALSE]
	}
	
	if(cutoff < 1) {
		# p = numeric(nrow(mat))
		# for(j in seq_len(nrow(mat))) {
		# 	group = factor(subgroup)
		# 	data = mat[j, ]
		# 	data = data + rnorm(length(data), 0, 0.01)
		# 	df = data.frame(data = data, group = group)
		# 	p[j] = oneway.test(data ~ group, data = df)$p.value
		# }
		data = mat
		for(j in seq_len(ncol(mat))) {
			data[, j] = data[, j] + rnorm(nrow(mat), 0, 0.01)
		}
		p = rowFtests(mat, factor(subgroup))[, "p.value"]
		adjp = p.adjust(p, method = adj_method)
		l = adjp < cutoff

		gr$p_value = p
		gr$adjp = adjp

		gr = gr[l]
		mat = mat[l, , drop = FALSE]
			
		message(qq("remove @{sum(!l)}/@{length(l)} rows after filtering by oneway-ANOVA test (@{adj_method} < @{cutoff})."))
	}
	
	ogr = gr
	# don't need that much rows
	nr = nrow(mat)
	if(nr > 5000) {
		l = sort(sample(nr, 5000))
		gr = gr[l]
		mat = mat[l, , drop = FALSE]
	}

	### draw heatmap ###

	### cluster columns ###########
	if(cluster_columns == "subgroup") {
		type = unique(subgroup)
		mm = NULL
		for(t in type) {
			m1 = mat[, subgroup == t, drop = FALSE]
			hc = hclust(dist(t(m1)))
			m1 = m1[, hc$order, drop = FALSE]
			cn = c(colnames(mm), colnames(m1))
			mm = cbind(mm, m1)
			colnames(mm) = cn
		}
		mat = mm
		cluster_columns = FALSE
	} else if(cluster_columns == "all") {
		cluster_columns = TRUE
	} else {
		cluster_columns = FALSE
	}

	ht_list = Heatmap(mat, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
		top_annotation = ha,
		cluster_columns = cluster_columns, show_row_dend = FALSE,
		use_raster = nrow(mat) > 2500, raster_quality = 2,
		heatmap_legend_param = list(title = "methylation"), ...)

	ncpg = gr$ncpg
	ncpg[ncpg > quantile(ncpg, 0.95)] = quantile(ncpg, 0.95)
	ht_list = ht_list + rowAnnotation(ncpg = row_anno_points(ncpg, axis = TRUE, axis_side = "top", size = unit(0.5, "mm"), gp = gpar(col = "#00000040")), width = unit(1, "cm"),
		show_annotation_name = TRUE)
	
	# make a copy of `gr`
	gr2 = gr
	mcols(gr2) = NULL

	# len = width(gr2)
	# len[len > quantile(len, 0.95)] = quantile(len, 0.95)
	# if(length(unique(len)) > 1) {
	# 	ht_list = ht_list + rowAnnotation(length = row_anno_points(len, axis = TRUE, axis_side = "top", size = unit(0.5, "mm"), gp = gpar(col = "#00000040")), width = unit(1, "cm"),
	# 		show_annotation_name = TRUE)
	# }
	
	if(!is.null(genomic_features)) {
		if(inherits(genomic_features, "GRanges")) {
			gf_name = deparse(substitute(genomic_features))
			genomic_features = list(genomic_features)
			names(genomic_features) = gf_name

		}
		gr_combine = annotate_to_genomic_features(gr2, genomic_features, prefix = "")
		ht_list = ht_list + Heatmap(as.matrix(mcols(gr_combine)), name = "anno", col = colorRamp2(c(0, 1), c("white", "orange")),
			cluster_columns = FALSE, use_raster = length(gr2) > 2500, raster_quality = 2,
			width = unit(4, "mm")*length(genomic_features))
	}
	
	draw(ht_list)
	
	return(invisible(ogr))
}


# == title
# Calculate mean methylation in a list of genomic features
#
# == param
# -sample_id  a vector of sample IDs
# -genomic_features a list or a single `GenomicRanges::GRanges` objects
# -chromosome a vector of chromosome names
#
# == value
# A list of or a single `GenomicRanges::GRanges` objects (according to ``genomic_features`` you specified) in which mean methylation matrix and number of CpG in each region
# are attached. The variable can be sent to `heatmap_diff_methylation_in_genomic_features` to visualize.
#
# Note it should be kept in mind that it doesn't make any sense to calculate mean methylation in long regions where
# there are hetergenuous methylation patterns.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
get_mean_methylation_in_genomic_features = function(sample_id, genomic_features, chromosome = paste0("chr", 1:22)) {
	
	is_single_gf = FALSE
	if(inherits(genomic_features, "GRanges")) {
		gf_name = deparse(substitute(genomic_features))
		genomic_features = list(genomic_features)
		names(genomic_features) = gf_name
		is_single_gf = TRUE
	}

	genomic_features = lapply(genomic_features, function(gf) {
		gf = gf[seqnames(gf) %in% chromosome]
		# mcols(gf) = NULL
		gf
	})

	mean_meth_list = lapply(genomic_features, function(gf) {
		m = matrix(NA, nrow = length(gf), ncol = length(sample_id))
		colnames(m) = paste0("mean_meth_", sample_id)
		m
	})
	ncpg_list = lapply(genomic_features, function(gf) {
		rep(0, length(gf))
	})

	for(chr in chromosome) {
		methylation_hooks$set_chr(chr, verbose = FALSE)
		meth_mat = methylation_hooks$meth[, sample_id, drop = FALSE]
		meth_gr = methylation_hooks$gr
		
		for(i in seq_along(genomic_features)) {
			message(qq("overlapping to @{names(genomic_features)[i]} on @{chr}"))
			mtch = as.matrix(findOverlaps(genomic_features[[i]], meth_gr))
			mean_meth = tapply(mtch[,2], mtch[,1], function(i) colMeans(meth_mat[i, , drop = FALSE], na.rm = TRUE))
			ncpg = tapply(mtch[,2], mtch[,1], length)
			ind = as.integer(names(mean_meth))
			if(is.list(mean_meth)) {
				mean_meth_list[[i]][ind, ] = do.call("rbind", mean_meth)
			} else {
				mean_meth_list[[i]][ind, ] = mean_meth
			}
			ncpg_list[[i]][ind] = ncpg
		}
	}

	for(i in seq_along(genomic_features)) {
		mcols(genomic_features[[i]]) = cbind(as.data.frame(mcols(genomic_features[[i]])), mean_meth_list[[i]], ncpg = ncpg_list[[i]])
	}

	if(is_single_gf) {
		genomic_features = genomic_features[[1]]
	}

	return(genomic_features)
}

get_mean_histone_density_in_genomic_features = function(mark, sample_id, genomic_features, chromosome = paste0("chr", 1:22)) {
	
	is_single_gf = FALSE
	if(inherits(genomic_features, "GRanges")) {
		gf_name = deparse(substitute(genomic_features))
		genomic_features = list(genomic_features)
		names(genomic_features) = gf_name
		is_single_gf = TRUE
	}

	genomic_features = lapply(genomic_features, function(gf) {
		gf = gf[seqnames(gf) %in% chromosome]
		# mcols(gf) = NULL
		gf
	})

	sample_id = intersect(chipseq_hooks$sample_id(mark), sample_id)

	mean_dend_list = lapply(genomic_features, function(gf) {
		m = matrix(NA, nrow = length(gf), ncol = length(sample_id))
		colnames(m) = qq("mean_@{mark}_@{sample_id}", collapse = FALSE)
		m
	})

	for(sid in sample_id) {
		peak = chipseq_hooks$peak(mark, sid)
		for(i in seq_along(genomic_features)) {
			message(qq("overlapping to @{names(genomic_features)[i]} on @{sid}"))
			mean_dend_list[[i]][, qq("mean_@{mark}_@{sid}")] = overlap_mean_signals(genomic_features[[i]], peak, peak$density, "w0", 0)
		}
	}
	
	for(i in seq_along(genomic_features)) {
		mcols(genomic_features[[i]]) = cbind(as.data.frame(mcols(genomic_features[[i]])), mean_dend_list[[i]])
	}

	if(is_single_gf) {
		genomic_features = genomic_features[[1]]
	}

	return(genomic_features)
}

overlap_mean_signals = function(gr1, gr2, value, mean_mode = c("w0", "absolute", "weighted"), empty_value = 0) {
	mean_mode = match.arg(mean_mode)[1]
	mtch = as.matrix(findOverlaps(gr1, gr2))
	if(is.atomic(value) && length(value) == 1) {
		value = mcols(gr2)[, value]
	}
	if(!is.null(ncol(value))) {
		value = as.data.frame(value)
	}
	v = HilbertCurve:::average_in_window(gr1, gr2, mtch, value, mean_mode, empty_value)
	v2 = matrix(empty_value, nc = ncol(v), nr = length(gr1))
	v2[unique(mtch[, 1]), ] = v
	colnames(v2) = colnames(v)
	return(v2)
}
