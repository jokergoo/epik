normalize_epigenomic_signals = function(cr, target, marks, expr, include_correlation_matrix = TRUE,
	extend = 5000, target_ratio = 0.1) {

	cr_param = metadata(cr)$cr_param
	sample_id = cr_param$sample_id
	subgroup = cr_param$subgroup
	subgroup_level = unique(subgroup)
	n_subgroup = length(subgroup_level)

	############ enriched to methylations ##################
	message("normalizing methylation to @{deparse(substitute(target))}")

	if(include_correlation_matrix) {
		mat_meth_corr = normalizeToMatrix(cr, target, mapping_column = "gene_id", value_column = "corr",
			extend = extend, mean_mode = "absolute", w = 50, target_ratio = target_ratio, trim = 0, empty_value = 0)
	} else {
		mat_meth_corr = NA
	}

	meth_mat_mean = enrich_with_methylation(target, sample_id, target_ratio = target_ratio, extend = extend)
	meth_mat_mean[attr(meth_mat_mean, "failed_rows"), ] = 0.5

	if(n_subgroup <= 1) {
		meth_mat_diff = enrich_with_methylation(target, sample_id, target_ratio = target_ratio, extend = extend, mode = rowIQRs)
		meth_mat_diff[attr(meth_mat_diff, "failed_rows"), ] = 0
	} else if(n_subgroup == 2) {
		meth_mat_mean_1 = enrich_with_methylation(target, sample_id[subgroup == subgroup_level[1]], target_ratio = target_ratio, extend = extend)
		failed_rows = attr(meth_mat_mean_1, "failed_rows")
		qqcat("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets.\n")
		meth_mat_mean_1[failed_rows, ] = 0.5

		meth_mat_mean_2 = enrich_with_methylation(target, sample_id[subgroup == subgroup_level[2]], target_ratio = target_ratio, extend = extend)
		failed_rows = attr(meth_mat_mean_2, "failed_rows")
		qqcat("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets.\n")
		meth_mat_mean_2[failed_rows, ] = 0.5
		meth_mat_diff = meth_mat_1 - meth_mat_mean_2
	} else {
		meth_mat_mean_list = lapply(subgroup_level, function(le) {
			meth_mat_mean = enrich_with_methylation(target, sample_id[subgroup == le, target_ratio = target_ratio, extend = extend)
			failed_rows = attr(meth_mat_mean, "failed_rows")
			qqcat("There are @{length(failed_rows)} failed rows when normalizing methylation to the targets.\n")
			meth_mat_mean[failed_rows, ] = 0.5
			meth_mat_mean
		})
		meth_array_mean = array(dim = c(dim(meth_mat_mean), n_subgroup))
		for(i in seq_along(meth_mat_mean_list)) {
			meth_array_mean[, , i] = meth_mat_mean_list[[i]]
		}
		meth_mat_diff = apply(meth_array_mean, c(1, 2), max)  - apply(meth_array_mean, c(1, 2), min)
		meth_mat_diff = copyAttr(meth_mat_mean, meth_mat_diff)
	}


	################# enrich to histome modifications #################
	hist_mat_corr_list = list()
	hist_mat_mean_list = list()
	hist_mat_diff_list = list()
	
	for(k in seq_along(marks)) {
		
		message(qq("normalizing @{marks[i]} signals to @{deparse(substitute(target))}"))

		hm_sample_id = intersect(sample_id, chipseq_hooks$sample_id(marks[k]))
		hm_subgroup = subgroup[sample_id %in% hm_sample_id]
		hm_subgroup_level = unique(hm_subgroup)
		n_hm_subgroup = length(hm_subgroup_level)

		# applied to each sample, each mark
		lt = enrich_with_histone_mark(target, sample_id = hm_sample_id, mark = marks[k], return_arr = TRUE, 
			target_ratio = target_ratio, extend = extend)
		hist_mat_mean_list[[k]] = lt[[2]]
		arr = lt[[1]]

		# only calculate the correlation when there are enough samples
		if(length(hm_sample_id) >= 5 && include_correlation_matrix) {
			# detect regions that histone MARKS correlate to expression
			expr2 = expr[target$gene_id, intersect(colnames(expr), hm_sample_id)]
			hist_mat_corr = matrix(nrow = nrow(expr2), ncol = ncol(mat))

			counter = set_counter(nrow(hist_mat_corr), fmt = "  calculate correlation for %s rows.")
			for(i in seq_len(nrow(hist_mat_corr))) {
				counter()
			    for(j in seq_len(ncol(hist_mat_corr))) {
			        x = cor(arr[i, j, ], expr2[i, ], method = "spearman")
			        hist_mat_corr[i, j] = x
			    }
			}
			cat("\n")
			hist_mat_corr[is.na(hist_mat_corr)] = 0
			hist_mat_corr = copyAttr(meth_mat_mean, hist_mat_corr)
			cor_mat_list[[k]] = hist_mat_corr
		} else {
			hist_mat_corr_list[[k]] = NA
		}

		if(n_hm_subgroup <= 1) {
			hist_mat_diff_list[[k]] = apply(arr, c(1, 2), IQR, na.rm = TRUE)
			hist_mat_diff_list[[k]] = copyAttr(meth_mat_mean, hist_mat_diff_list[[k]])
		} else if(n_hm_subgroup == 2) {
			hm_sample_id_subgroup1 = hm_sample_id[hm_subgroup == hm_subgroup_level[1]]
			hm_sample_id_subgroup2 = hm_sample_id[hm_subgroup == hm_subgroup_level[2]]
			h1 = apply(arr[, , hm_sample_id_subgroup1], c(1, 2), mean, na.rm = TRUE)
			h2 = apply(arr[, , hm_sample_id_subgroup2], c(1, 2), mean, na.rm = TRUE)
			hist_mat_diff_list[[k]] = h1 - h2
			hist_mat_diff_list[[k]] = copyAttr(meth_mat_mean, hist_mat_diff_list[[k]])
		} else {
			h_list = lapply(hm_subgroup_level, function(le) {
				apply(arr[, , hm_sample_id[hm_subgroup == le]], c(1, 2), mean, na.rm = TRUE)
			})
			hm_array_mean = array(dim = c(dim(meth_mat_mean), n_hm_subgroup))
			for(i in seq_along(h_list)) {
				hm_array_mean[, , i] = h_list[[i]]
			}
			hist_mat_diff_list[[k]] = apply(hm_array_mean, c(1, 2), max)  - apply(hm_array_mean, c(1, 2), min)
			hist_mat_diff_list[[k]] = copyAttr(meth_mat_mean, hist_mat_diff_list[[k]])
		}
	}
	names(hist_mat_corr_list) = marks
	names(hist_mat_mean_list) = marks
	names(hist_mat_diff_list) = marks

	return(list(meth_mat_corr = meth_mat_corr,
		        meth_mat_mean = meth_mat_mean,
		        meth_mat_diff = meth_mat_diff,
		        hist_mat_corr_list = hist_mat_corr_list,
		        hist_mat_mean_list = hist_mat_mean_list,
		        hist_mat_diff_list = hist_mat_diff_list))
}

