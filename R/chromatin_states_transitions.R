# == title
# Generate transition matrix from chromHMM results
#
# == param
# -gr_list_1 a list of `GenomicRanges::GRanges` objects which contain chromatin states in group 1.
#            The first column in meta columns should be the states. Be careful when importing bed files to 
#            `GenomicRanges::GRanges` objects (start positions in bed files are 0-based while 1-based in ``GRanges`` objects).
# -gr_list_2 a list of `GenomicRanges::GRanges` objects which contains chromatin states in group 2.
# -window window size which was used to do chromHMM states segmentation If it is not specified, the greatest common divisor
#         of the width of all regions is used.
# -min_1  If there are multiple samples in the group, it is possible that a segment has more than one states asigned to it.
#         If the recurrency of each state is relatively low, it means there is no one dominant state for this segment and it should 
#         be removed. This argument controls the minimal value for the recurrency of states in a given segment.
# -min_2 same as ``min_1``, but for samples in group 2.
# -meth_diff If methylation dataset is provided, the segments for which the methylation difference between two groups is less than
#             this value are removed.
# -chromosome subset of chromosomes
#
# == detail
# For a segment in the genome, the chromatin state may be different in different subgroups. This is called chromatin state transistion.
# This function visualize such kind of genome-wide transitions.
#
# The whole genome is segmentated by ``window`` and states with highest occurence among samples are assigned to segments.
#
# To make the function run successfully, number of segments (after binned by ``window``) in all samples 
# should be all the same and there should be no gaps between segments. If the segmentation data is directly imported from
# chromHMM results, you dont need to worry.
#
# == value
# A transition matrix in which values represent total width of segments that transite from one state to the other in the two groups. Rows correspond
# to group 1 and columns correspond to group 2.
#
# If methylation dataset is provided, the mean methylation for each state in each group is attached, which will be used to calculate
# mean methylation difference in `chromatin_states_transition_chord_diagram`.
#
# == seealso
# The matrix can be sent to `chromatin_states_transition_chord_diagram` to visualize.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# set.seed(123)
# gr_list_1 = lapply(1:5, function(i) {
# 	pos = sort(c(0, sample(1:9999, 99), 10000))*200
# 	GRanges(seqnames = "chr1", ranges = IRanges(pos[-101] + 1, pos[-1]), 
# 	    states = paste0("state_", sample(1:9, 100, replace = TRUE)))
# })
# gr_list_2 = lapply(1:5, function(i) {
# 	pos = sort(c(0, sample(1:9999, 99), 10000))*200
# 	GRanges(seqnames = "chr1", ranges = IRanges(pos[-101] + 1, pos[-1]), 
# 	    states = paste0("state_", sample(1:9, 100, replace = TRUE)))
# })
# mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2)
make_transition_matrix_from_chromHMM = function(gr_list_1, gr_list_2, window = NULL, 
	min_1 = floor(length(gr_list_1)/2), min_2 = floor(length(gr_list_2)/2), 
	meth_diff = 0, chromosome = paste0("chr", 1:22)) {

	if(inherits(gr_list_1, "GRanges")) {
		gr_list_1 = list(gr_list_1)
		min_1 = 0
	}
	if(inherits(gr_list_2, "GRanges")) {
		gr_list_2 = list(gr_list_2)
		min_2 = 0
	}

	if(!is.null(chromosome)) {
		gr_list_1 = lapply(gr_list_1, function(gr) gr[seqnames(gr) %in% chromosome])
		gr_list_2 = lapply(gr_list_2, function(gr) gr[seqnames(gr) %in% chromosome])
	}

	if(is.null(window)) {
		window = numbers::mGCD(c(sapply(gr_list_1, function(gr) numbers::mGCD(unique(width(gr)))),
		                sapply(gr_list_2, function(gr) numbers::mGCD(unique(width(gr))))))
		message(qq("window is set to @{window}"))
		if(window == 1) {
			message("when converting bed files to GRanges objects, be careful with the 0-based and 1-based coordinates.")
		}
	}

	# if(is.null(all_states)) {
		all_states = unique(c(unlist(lapply(gr_list_1, function(gr) {unique(as.character(mcols(gr)[, 1]))})),
			                unlist(lapply(gr_list_2, function(gr) {unique(as.character(mcols(gr)[, 1]))}))))
		all_states = sort(all_states)
		message(qq("@{length(all_states)} states in total"))
	# }

	message("extracting states")
	m1 = as.data.frame(lapply(gr_list_1, function(gr) {
		k = round(width(gr)/window)
		s = as.character(mcols(gr)[, 1])
		as.integer(factor(rep(s, times = k), levels = all_states))
	}))
	m2 = as.data.frame(lapply(gr_list_2, function(gr) {
		k = round(width(gr)/window)
		s = as.character(mcols(gr)[, 1])
		as.integer(factor(rep(s, times = k), levels = all_states))
	}))

	m1 = as.matrix(m1)
	m2 = as.matrix(m2)

	gr = makeWindows(gr_list_1[[1]], w = window)

	message("counting for each state")
	t1 = rowTabulates(m1)
	t2 = rowTabulates(m2)
	l = rowMaxs(t1) >= min_1 & rowMaxs(t2) >= min_2
	t1 = t1[l, , drop = FALSE]
	t2 = t2[l, , drop = FALSE]
	gr = gr[l]

	message("determine states")
	states1 = rowWhichMax(t1)
	states2 = rowWhichMax(t2)

	meth_hooks_defined = !is.null(methylation_hooks$get_by_chr)

	if(meth_hooks_defined) {
		## get sum_meth and n_cpg in each window
		sum_meth1 = rep(NA, length(states1))
		sum_meth2 = rep(NA, length(states1))
		n_cpg = numeric(length(states1))
		for(chr in unique(seqnames(gr))) {
			l_chr = as.vector(seqnames(gr)) == chr
			gr_subset = gr[l_chr]
			
			methylation_hooks$set_chr(chr, verbose = FALSE)
			meth_gr = methylation_hooks$gr
			meth_mat = methylation_hooks$meth
			
			mtch = as.matrix(findOverlaps(gr_subset, meth_gr))

			message("calculating mean methylation for transistions from `gr_list_1`")
			x = tapply(mtch[, 2], mtch[, 1], function(ind) {
				sum(rowMeans(meth_mat[ind, intersect(names(gr_list_1), colnames(meth_mat)), drop = FALSE]))
			})
			sum_meth1[which(l_chr)[as.numeric(names(x))]] = x
			
			message("calculating mean methylation for transistions from `gr_list_2`")
			x = tapply(mtch[, 2], mtch[, 1], function(ind) {
				sum(rowMeans(meth_mat[ind, intersect(names(gr_list_2), colnames(meth_mat)), drop = FALSE]))
			})
			sum_meth2[which(l_chr)[as.numeric(names(x))]] = x
			n_cpg[which(l_chr)[as.numeric(names(x))]] = tapply(mtch[, 2], mtch[, 1], length)
		}

		l = abs(sum_meth1/n_cpg - sum_meth2/n_cpg) >= meth_diff
		l[is.na(l)] = FALSE
		states1 = states1[l]
		states2 = states2[l]
		sum_meth1 = sum_meth1[l]
		sum_meth2 = sum_meth2[l]
		n_cpg = n_cpg[l]
	}

	message("generate transition matrix")
	mat = as.matrix(table(states1, states2))
	rownames(mat) = all_states[as.numeric(rownames(mat))]
	colnames(mat) = all_states[as.numeric(colnames(mat))]
	mat = mat * as.numeric(window)
	class(mat) = "matrix"
	names(dimnames(mat)) = NULL

	mat2 = matrix(0, nrow = length(all_states), ncol = length(all_states))
	rownames(mat2) = all_states
	colnames(mat2) = all_states
	mat2[rownames(mat), colnames(mat)] = mat

	meth_mean_1 = NULL
	meth_mean_2 = NULL
	if(meth_hooks_defined) {
		meth_mean_1 = matrix(NA, nrow = nrow(mat2), ncol = ncol(mat2))
		dimnames(meth_mean_1) = dimnames(mat2)
		meth_mean_2 = matrix(NA, nrow = nrow(mat2), ncol = ncol(mat2))
		dimnames(meth_mean_2) = dimnames(mat2)
		for(s1 in rownames(mat2)) {
			for(s2 in colnames(mat2)) {
				l1 = all_states[states1] == s1 & all_states[states2] == s2 & !is.na(sum_meth1)
				if(sum(l1)) meth_mean_1[s1, s2] = mean(sum_meth1[l1]/n_cpg[l1])
				l2 = all_states[states1] == s1 & all_states[states2] == s2 & !is.na(sum_meth2)
				if(sum(l2)) meth_mean_2[s1, s2] = mean(sum_meth2[l2]/n_cpg[l2])
			}
		}
	}

	attr(mat2, "meth_mean_1") = meth_mean_1
	attr(mat2, "meth_mean_2") = meth_mean_2

	class(mat2) = c("chromatin_states_transition_matrix", "matrix")

	return(mat2)
}

# == title
# Print chromatin_states_transition_matrix class object
#
# == param
# -x a ``chromatin_states_transition_matrix`` class object
# -... additional arguments
#
# == value
# no value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
print.chromatin_states_transition_matrix = function(x, ...) {
	qqcat("A chromatin states transition matrix defined by @{nrow(x)} states:\n")
	print(rownames(x))

	if(!is.null(attr(x, "meth_mean_1"))) {
		cat("Mean methylation in each transition provided\n")
	}
}

# == title
# Subset chromatin_states_transition_matrix class object
#
# == param
# -x a ``chromatin_states_transition_matrix`` class object
# -i index of rows
# -j index of columns
# -drop whether degenerate the matrix
#
# == value
# a ``chromatin_states_transition_matrix`` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
"[.chromatin_states_transition_matrix" = function(x, i, j, drop = FALSE) {
	
	meth_mean_1 = attr(x, "meth_mean_1")
	meth_mean_2 = attr(x, "meth_mean_2")
	class(x) = "matrix"
	n = nrow(x)

	if(missing(i) && missing(j)) {
		i = 1:n
		j = 1:n
	} else if(missing(i)) {
		i = 1:n
	} else if(missing(j)) {
		j = 1:n
	}
		
	if(nargs() == 2) {
		return(x[i])
	}

	x = x[i, j, drop = FALSE]
	if(!is.null(meth_mean_1)) {
		meth_mean_1 = meth_mean_1[i, j, drop = FALSE]
		meth_mean_2 = meth_mean_2[i, j, drop = FALSE]
	}
	attr(x, "meth_mean_1") = meth_mean_1
	attr(x, "meth_mean_2") = meth_mean_2
	class(x) = c("chromatin_states_transition_matrix", "matrix")
	return(x)
}

# == title
# Transpose the chromatin transition matrix
#
# == param
# -x a ``chromatin_states_transition_matrix`` class object
#
# == value
# a ``chromatin_states_transition_matrix`` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
t.chromatin_states_transition_matrix = function(x) {
	meth_mean_1 = attr(x, "meth_mean_1")
	meth_mean_2 = attr(x, "meth_mean_2")
	class(x) = "matrix"

	if(is.null(meth_mean_1)) {
		x = t(x)
	} else{
		x = t(x)
		meth_mat_1 = t(meth_mat_1)
		meth_mat_2 = t(meth_mat_2)
	}
	attr(x, "meth_mean_1") = meth_mean_1
	attr(x, "meth_mean_2") = meth_mean_2
	class(x) = c("chromatin_states_transition_matrix", "matrix")
	return(x)
}

# == title
# Simply return names of chromatin states
#
# == param
# -x a ``chromatin_states_transition_matrix`` object returned from `make_transition_matrix_from_chromHMM`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
state_names = function(x) {
	return(rownames(x))
}



# == title
# Change chromatin state names
#
# == param
# -x a ``chromatin_states_transition_matrix`` object returned from `make_transition_matrix_from_chromHMM`
# -value new chromatin state names
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
"state_names<-" = function(x, value) {
	meth_mean_1 = attr(x, "meth_mean_1")
	meth_mean_2 = attr(x, "meth_mean_2")
	class(x) = "matrix"
	n = nrow(x)

	rownames(x) = value
	colnames(x) = value
	if(!is.null(meth_mean_1)) {
		rownames(meth_mean_1) = value
		colnames(meth_mean_1) = value
		rownames(meth_mean_2) = value
		colnames(meth_mean_2) = value
	}
	attr(x, "meth_mean_1") = meth_mean_1
	attr(x, "meth_mean_2") = meth_mean_2
	class(x) = c("chromatin_states_transition_matrix", "matrix")
	return(x)
}

# == title
# Chord diagram for visualizing chromatin states transitions
#
# == param
# -mat the transition matrix. It should be a square matrix in which row names and column names should be all the same.
#      If it is not, the function will try to re-format it.
# -group_names name for the two groups under comparison. You also add it afterwards by using `graphics::text`
# -max_mat if there are several transition matrix to be compared, set it to the matrix with maximum absolute and it will make
#          scales of all matrix the same and comparable.
# -remove_unchanged_transition whether to remove transitions that states are not changed (set the values in diagonal to 0)
# -state_col color for states. It should be a vector of which names correspond to states.
# -legend_position positions of legends. Possible values are "bottomleft", "bottomright", "topright" and "topleft".
#             If the value is specified as vector with length larger than two, the legend will be split into several parts.
#              Set the value to ``NULL`` to suppress legends.
# -... pass to `circlize::chordDiagram`
#
# == details
# Rows of ``mat`` locate at the bottom of the circle by default. You can transpose the matrix to move rows to the top of the circle.
#
# The chord diagram visualizes how much chromatin states change. In the diagram, width of each link represents the total
# width of segments in a certain chromatin state in group 1 that transite to other chromatin state in group 2. The width of 
# each grid represents total width of segments in a certain chromatin in group 1 that transite to all states in group 2.
#
# If methylation dataset is provided when making the transistion matrix by using `make_transition_matrix_from_chromHMM`,
# there will be extra tracks on the outside of the circlie to represenst the mean methylation difference in two groups.
#
# Chord diagram is implemented in base graphic system, which means, you can add titles or other graphics by base graphic 
# functions (e.g. `graphics::title`, `graphics::text`, ...)
#
# If you want to adjust order of states in the chord diagram, directly change row and column order of the matrix.
# 
# == value
# No value is returned.
#
# == seealso
# `make_transition_matrix_from_chromHMM` which generates transition matrix directly from chromHMM results.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# set.seed(123)
# gr_list_1 = lapply(1:5, function(i) {
# 	pos = sort(c(0, sample(1:9999, 99), 10000))*200
# 	GRanges(seqnames = "chr1", ranges = IRanges(pos[-101] + 1, pos[-1]), 
# 	    states = paste0("state_", sample(1:9, 100, replace = TRUE)))
# })
# gr_list_2 = lapply(1:5, function(i) {
# 	pos = sort(c(0, sample(1:9999, 99), 10000))*200
# 	GRanges(seqnames = "chr1", ranges = IRanges(pos[-101] + 1, pos[-1]), 
# 	    states = paste0("state_", sample(1:9, 100, replace = TRUE)))
# })
# mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2)
# chromatin_states_transition_chord_diagram(mat, legend_position = "bottomleft")
chromatin_states_transition_chord_diagram = function(mat, group_names = NULL, max_mat = mat, 
	remove_unchanged_transition = TRUE, state_col = NULL, legend_position = NULL, ...) {

	op = par(no.readonly = TRUE)
	on.exit(par(op))
	par(xpd = NA)

	meth_mean_1 = attr(mat, "meth_mean_1")
	meth_mean_2 = attr(mat, "meth_mean_2")

	contain_methylation = !is.null(meth_mean_1)

	if(contain_methylation) {
		# only those with big methylation difference
		meth_mean_1 = attr(mat, "meth_mean_1")
		meth_mean_2 = attr(mat, "meth_mean_2")
	}

	mat[is.na(mat)] = 0

	if(length(intersect(rownames(mat), colnames(mat))) == 0) {
		rownames(mat) = paste0("E", seq_len(nrow(mat)))
		colnames(mat) = paste0("E", seq_len(ncol(mat)))
		if(contain_methylation) {
			dimnames(meth_mean_1) = dimnames(mat)
			dimnames(meth_mean_2) = dimnames(mat)
		}
	}

	all_states = union(rownames(mat), colnames(mat))
	n_states = length(all_states)
	mat2 = matrix(0, nrow = n_states, ncol = n_states)
	rownames(mat2) = all_states
	colnames(mat2) = all_states
	mat2[rownames(mat), colnames(mat)] = mat
	if(contain_methylation) {
		meth_mean_foo_1 = matrix(NA, nrow = n_states, ncol = n_states)
		dimnames(meth_mean_foo_1) = dimnames(mat2)
		meth_mean_foo_2 = meth_mean_foo_1
		meth_mean_foo_1[rownames(meth_mean_1), colnames(meth_mean_1)] = meth_mean_1
		meth_mean_foo_2[rownames(meth_mean_2), colnames(meth_mean_2)] = meth_mean_2
	}

	mat = mat2
	if(contain_methylation) {
		meth_mean_1 = meth_mean_foo_1
		meth_mean_2 = meth_mean_foo_2
	}

	if(is.null(state_col)) {
		col_fun = colorRamp2(0:10, brewer.pal(11, "Spectral"))
		x = (seq_along(all_states)-1)*10/(length(all_states)-1)
		state_col = col_fun(x)
		names(state_col) = all_states
	} else {
		if(is.null(names(state_col))) {
			names(state_col) = all_states
		}
	}

	if(remove_unchanged_transition) {
		for(i in all_states) {
			mat[i, i] = 0
		}
	}

	rownames(mat) = paste0("R_", seq_len(n_states))
	colnames(mat) = paste0("C_", seq_len(n_states))
	if(contain_methylation) {
		dimnames(meth_mean_1) = dimnames(mat)
		dimnames(meth_mean_2) = dimnames(mat)
	}

	state_col2 = c(state_col, state_col)
	names(state_col2) = c(rownames(mat), colnames(mat))

	colmat = rep(state_col2[rownames(mat)], n_states)
	colmat = rgb(t(col2rgb(colmat)), maxColorValue = 255)
	
	# make thin links very light
	qati = quantile(mat, 0.7)
	colmat[mat > qati] = paste0(colmat[mat > qati], "A0")
	colmat[mat <= qati] = paste0(colmat[mat <= qati], "20")
	dim(colmat) = dim(mat)

	de = 360 - (360 - 20 - 30) * sum(mat)/sum(max_mat) - 30
	circos.par(start.degree = -de/4, gap.degree = c(rep(1, n_states-1), de/2, rep(1, n_states-1), de/2),
		points.overflow.warning = FALSE)

	cdm_res = chordDiagram(mat, col = colmat, grid.col = state_col2,
		directional = TRUE, annotationTrack = "grid", preAllocateTracks = list(track.height = ifelse(contain_methylation, 0.1, 0.01)), ...)

	for(sn in get.all.sector.index()) {
		if(abs(get.cell.meta.data("cell.start.degree", sector.index = sn) - get.cell.meta.data("cell.end.degree", sector.index = sn)) > 3) {
			xcenter = get.cell.meta.data("xcenter", sector.index = sn, track.index = 2)
			ycenter = get.cell.meta.data("ycenter", sector.index = sn, track.index = 2)
			i_state = as.numeric(gsub("(C|R)_", "", sn))
			circos.text(xcenter, ycenter, i_state, col = "white", font = 2, cex = 0.7, 
				sector.index = sn, track.index = 2, adj = c(0.5, 0.5), niceFacing = TRUE)
			circos.axis(sector.index = sn, track.index = 2, major.tick.percentage = 0.2, labels.away.percentage = 0.2, labels.cex = 0.5)
		}
	}

	if(contain_methylation) {
		oljoin = par("ljoin")
	    par(ljoin = "mitre")
		abs_max = quantile(abs(c(meth_mean_1, meth_mean_2) - 0.5), 0.95, na.rm = TRUE)
		col_fun = colorRamp2(c(0.5 - abs_max, 0.5, 0.5 + abs_max), c("blue", "white", "red"))
		col_fun2 = colorRamp2(c(-abs_max, 0, abs_max), c("green", "white", "orange"))

		ylim = get.cell.meta.data("ylim", sector.index = rownames(mat)[1], track.index = 1)
		y1 = ylim[1] + (ylim[2] - ylim[1])*0.4
		y2 = ylim[2]
		for(i in seq_len(nrow(cdm_res))) {
			if(cdm_res$value[i] > 0) {
				circos.rect(cdm_res[i, "x1"], y1, cdm_res[i, "x1"] - abs(cdm_res[i, "value"]), y1 + (y2-y1)*0.45, 
					col = col_fun(meth_mean_1[cdm_res$rn[i], cdm_res$cn[i]]), border = col_fun(meth_mean_1[cdm_res$rn[i], cdm_res$cn[i]]),
					sector.index = cdm_res$rn[i], track.index = 1)

				circos.rect(cdm_res[i, "x1"], y1 + (y2-y1)*0.55, cdm_res[i, "x1"] - abs(cdm_res[i, "value"]), y2, 
					col = col_fun2(meth_mean_2[cdm_res$rn[i], cdm_res$cn[i]] - meth_mean_1[cdm_res$rn[i], cdm_res$cn[i]]), 
					border = col_fun2(meth_mean_2[cdm_res$rn[i], cdm_res$cn[i]] - meth_mean_1[cdm_res$rn[i], cdm_res$cn[i]]),
					sector.index = cdm_res$rn[i], track.index = 1)

				circos.rect(cdm_res[i, "x1"], -0.5, cdm_res[i, "x1"] - abs(cdm_res[i, "value"]), -0.7, 
					col = state_col2[cdm_res$cn[i]], border = state_col2[cdm_res$cn[i]],
					sector.index = cdm_res$rn[i], track.index = 2)

				circos.rect(cdm_res[i, "x2"], y1, cdm_res[i, "x2"] - abs(cdm_res[i, "value"]), y1 + (y2-y1)*0.45, 
					col = col_fun(meth_mean_2[cdm_res$rn[i], cdm_res$cn[i]]), 
					border = col_fun(meth_mean_2[cdm_res$rn[i], cdm_res$cn[i]]),
					sector.index = cdm_res$cn[i], track.index = 1)

				circos.rect(cdm_res[i, "x2"], y1 + (y2-y1)*0.55, cdm_res[i, "x2"] - abs(cdm_res[i, "value"]), y2, 
					col = col_fun2(meth_mean_1[cdm_res$rn[i], cdm_res$cn[i]] - meth_mean_2[cdm_res$rn[i], cdm_res$cn[i]]), 
					border = col_fun2(meth_mean_1[cdm_res$rn[i], cdm_res$cn[i]] - meth_mean_2[cdm_res$rn[i], cdm_res$cn[i]]),
					sector.index = cdm_res$cn[i], track.index = 1)
			}
		}
		par(ljoin = oljoin)
	}

	if(length(legend_position) == 1) { 
		add_chord_diagram_legend(legend_position, 1:n_states, all_states, state_col[all_states])
	} else if(length(legend_position)== 2) {
		ib = ceiling(n_states/2)
		ind = 1:ib
		add_chord_diagram_legend(legend_position[1], ind, all_states[ind], state_col[all_states][ind])
		ind = (ib+1):n_states
		add_chord_diagram_legend(legend_position[2], ind, all_states[ind], state_col[all_states][ind])
	}
	if(!is.null(group_names)) {
		text(0, 1.1, group_names[2], adj = c(0.5, 0.5))
		text(0, -1.1, group_names[1], adj = c(0.5, 0.5))
	}
	text(-1, 1.1, paste0(sum(mat), " bp"), adj = c(0, 0.5))

	if(contain_methylation) {
		vps = baseViewports()
		pushViewport(vps$inner, vps$figure, vps$plot)

		at = round(c(0.5 - abs_max, 0.5, 0.5 + abs_max), digits = 1)
		lgd1 = Legend(at = at, labels = at, direction = "horizontal", col_fun = col_fun, title_position = "topcenter",
			title = "Mean methylation", legend_width = unit(3, "cm"))

		at = round(c(-abs_max, 0, abs_max), digits = 1)
		lgd2 = Legend(at = at, labels = at, direction = "horizontal", col_fun = col_fun2, title_position = "topcenter",
			title = "Mean difference", legend_width = unit(3, "cm"))

		pushViewport(viewport(x = unit(0.5, "npc") - unit(2.5, "cm"), y = unit(-3, "mm"), width = grobWidth(lgd1), 
			height = grobHeight(lgd1), just = c("right", "bottom")))
		grid.draw(lgd1)
		upViewport()
		pushViewport(viewport(x = unit(0.5, "npc") + unit(2.5, "cm"), y = unit(-3, "mm"), 
			width = grobWidth(lgd2), height = grobHeight(lgd2), just = c("left", "bottom")))
		grid.draw(lgd2)
		upViewport()
		upViewport()
	}

	circos.clear()
}


add_chord_diagram_legend = function(position = c("bottomleft", "bottomright", "topleft", "topright"), 
	index = seq_along(labels), labels, col) {
	
	position = match.arg(position)[1]
	if(length(index) == 0) {
		return(NULL)
	}

	coor = par("usr")
	n = length(labels)
	text_height = strheight("a")
	labels_max_width = max(strwidth(labels))
	legend_width = text_height*(1+0.5) + labels_max_width
	if(position == "bottomleft") {
		x1 = rep(coor[1], n)
		x2 = x1 + text_height
		y1 = coor[3] + (rev(seq_len(n))-1)*1.5*text_height
		y2 = y1 + text_height
		rect(x1, y1, x2, y2, col = col, border = col)
		text((x1+x2)/2, (y1+y2)/2, index, cex = 0.6, font = 2, col = "white")
		text(x2 + 0.5*text_height, (y1+y2)/2, labels, adj = c(0, 0.5), cex = 0.8)
	} else if(position == "bottomright") {
		x1 = rep(coor[2] - labels_max_width, n)
		x2 = x1 + text_height
		y1 = coor[3] + (rev(seq_len(n))-1)*1.5*text_height
		y2 = y1 + text_height
		rect(x1, y1, x2, y2, col = col, border = col)
		text((x1+x2)/2, (y1+y2)/2, index, cex = 0.6, font = 2, col = "white")
		text(x2 + 0.5*text_height, (y1+y2)/2, labels, adj = c(0, 0.5), cex = 0.8)
	} else if(position == "topleft") {
		x1 = rep(coor[1], n)
		x2 = x1 + text_height
		y1 = coor[4] - (seq_len(n)-1)*1.5*text_height
		y2 = y1 - text_height
		rect(x1, y1, x2, y2, col = col, border = col)
		text((x1+x2)/2, (y1+y2)/2, index, cex = 0.6, font = 2, col = "white")
		text(x2 + 0.5*text_height, (y1+y2)/2, labels, adj = c(0, 0.5), cex = 0.8)
	} else if(position == "topright") {
		x1 = rep(coor[2] - labels_max_width, n)
		x2 = x1 + text_height
		y1 = coor[4] - (seq_len(n)-1)*1.5*text_height
		y2 = y1 - text_height
		rect(x1, y1, x2, y2, col = col, border = col)
		text((x1+x2)/2, (y1+y2)/2, index, cex = 0.6, font = 2, col = "white")
		text(x2 + 0.5*text_height, (y1+y2)/2, labels, adj = c(0, 0.5), cex = 0.8)
	}
}

