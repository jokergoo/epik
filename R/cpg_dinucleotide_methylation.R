
two_rows_are_similar = function(x1, x2, meth_diff = 0.05, p = 0.8) {
	sum(abs(x1 - x2) < meth_diff)/length(x1) > p
}

two_row_means = function(x1, x2, w1, w2) {
	y = (x1*w1 + x2*w2)/(w1 + w2)
	y[is.na(y)] = 0
	y[is.infinite(y)] = 0
	y
}


# == title
# Merge methylation for CpG dinucleotide
#
# == param
# -pos positions of CpG sites, must be sorted
# -meth methylation matrix
# -cov CpG coverage matrix
# -meth_diff normally methylation for the two Cs in a CpG dinucleotide is very similar.
#    This cutoff helps to find whether two neighbouring Cs belong to a CpG dinucleotide.
# -p for a CpG dinucleotide, the minimal proportion of samples that the two Cs have methylation difference
#     less than ``meth_diff`
#
# == details
# normally methylation for the two Cs in a CpG dinucleotide is very similar. This function
# helps to reduce the redundency of the methylation dataset.
#
# For two Cs in a CpG dinucleotide, the merged methylation value is calculated by weighting
# the CpG coverage of the two Cs. Also the merged coverage is also calculated by weighting
# the coverage itself.
#
# == value
# a list containg a ``pos`` vector, a ``meth`` matrix and a ``cov`` matrix
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
# 
cpg_dinucleotide_methylation = function(pos, meth, cov, meth_diff = 0.05, p = 0.8) {

	if(length(pos) != nrow(meth)) {
		stop("length of `pos` should be same as the number of rows of `meth`.")
	}

	if(!identical(dim(meth), dim(cov))) {
		stop("dimension of `meth` should be the same as `cov`.")
	}

	meth = as.matrix(meth)
	cov = as.matrix(cov)

	meth[is.na(meth)] = 0
	cov[is.na(cov)] = 0
	cov[cov < 0] = 0
	cov[meth < 0] = 0

	# check the methylation is on CpG level or C level
	# if on C level, calcualte the mean methylation of two Cs weighted by coverage
	j = 0
	meth2 = matrix(nrow = nrow(meth), ncol = ncol(meth))
	cov2 = matrix(nrow = nrow(cov), ncol = ncol(cov))
	pos2 = numeric(length(pos))
	
	next_used = FALSE
	for(i in seq_len(nrow(meth))) {
		if(next_used) {
			next_used = FALSE
			next
		}

		if(i < length(pos)) {
			if(pos[i] >= pos[i + 1]) {
				stop("`pos` should be sorted.")
			}
		}

		if(pos[i] + 1 == pos[i + 1] && i < length(pos)) {
			# check whether the methylation values are same in all samples
			if(two_rows_are_similar(meth[i, ], meth[i+1, ], meth_diff, p)) {
				j = j + 1
				meth2[j, ] = two_row_means(meth[i, ], meth[i+1, ], cov[i, ], cov[i+1, ])
				cov2[j, ] = round(two_row_means(cov[i, ], cov[i+1, ], cov[i, ], cov[i+1, ]))
				pos2[j] = pos[i]
				next_used = TRUE
				qqcat("Merge CpG dinucleotide at position [@{pos[i]}, @{pos[i+1]}]\n")
			} else {
				j = j + 1
				meth2[j, ] = meth[i, ]
				cov2[j, ] = round(cov[i, ])
				pos2[j] = pos[i]
				qqcat("find an alone C in a continous CpG*** at position @{pos[i]}\n")
			}
		} else { 
			j = j + 1
			meth2[j, ] = meth[i, ]
			cov2[j, ] = round(cov[i, ])
			pos2[j] = pos[i]
			qqcat("find an alone C at position @{pos[i]}\n")
		}
	}
	l = apply(meth2, 1, function(x) all(is.na(x)))
	meth2 = meth2[!l, , drop = FALSE]
	cov2 = cov2[!l, , drop = FALSE]
	pos2 = pos2[!l]

	colnames(meth2) = colnames(meth)
	colnames(cov2) = colnames(cov)

	return(list(pos = pos, meth = meth2, cov = cov2))
}

