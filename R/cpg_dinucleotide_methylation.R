
weighted_mean = function(x1, x2, w1, w2) {
	y = ifelse(is.na(x1), x2, ifelse(is.na(x2), x1, (x1*w1 + x2*w2)/(w1 + w2)))
	y[is.infinite(y)] = NA
	y
}


# == title
# Merge methylation for CpG dinucleotide
#
# == param
# -pos positions of CpG sites, must be sorted
# -meth methylation matrix. Note the value should be the number of methylated CpGs at each CpG site
# -cov CpG coverage matrix
#
# == details
# Normally methylation for the two Cs in a CpG dinucleotide is very similar. This function
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
cpg_dinucleotide_methylation = function(pos, meth, cov) {

	if(length(pos) != nrow(meth)) {
		stop("length of `pos` should be same as the number of rows of `meth`.")
	}

	if(!identical(dim(meth), dim(cov))) {
		stop("dimension of `meth` should be the same as `cov`.")
	}

	if(all(meth <= 1)) {
		stop("`meth` should be the number of methylated CpGs.")
	}

	meth = as.matrix(meth)
	cov = as.matrix(cov)

	cov[is.na(cov)] = 0
	cov[cov < 0] = 0
	meth[is.na(meth)] = 0
	meth[meth < 0] = 0

	# check the methylation is on CpG level or C level
	# if on C level, calcualte the mean methylation of two Cs weighted by coverage
	j = 0
	meth2 = matrix(0, nrow = nrow(meth), ncol = ncol(meth))
	cov2 = matrix(0, nrow = nrow(cov), ncol = ncol(cov))
	pos2 = numeric(length(pos))

	n = length(pos)

	l_continous = rep(FALSE, n)
	l_continous[2:n] = l_continous[2:n] | pos[2:n] - pos[1:(n-1)] == 1
	l_continous[1:(n-1)] = l_continous[1:(n-1)] | pos[1:(n-1)] - pos[2:n] == -1

	l_single = !l_continous
	
	for(i in seq_len(n-1)) {
		if(!l_continous[i+1]) {
			l_continous[i] = FALSE
		}
		if(l_continous[i]) {
			l_continous[i+1] = FALSE
		}
	}
	if(l_continous[n]) {
		l_continous[n] = FALSE
	}

	pos2[l_single] = pos[l_single]
	meth2[l_single, ] = meth[l_single, ]
	cov2[l_single, ] = cov[l_single, ]

	ind = which(l_continous)
	pos2[ind] = pos[ind]
	meth2[ind, ] = weighted_mean(x1 = meth[ind, ], x2 = meth[ind+1, ],
		                         w1 = cov[ind, ], w2 = cov[ind+1, ])
	cov2[ind, ] = weighted_mean(x1 = cov[ind, ], x2 = cov[ind+1, ],
		                         w1 = cov[ind, ], w2 = cov[ind+1, ])
	
	l = pos2 > 0
	meth2 = meth2[l, , drop = FALSE]
	cov2 = cov2[l, , drop = FALSE]
	pos2 = pos2[l]

	cov2[is.na(cov2)] = 0
	meth2[is.na(meth2)] = 0

	cov2 = round(cov2)
	meth2 = round(meth2)

	colnames(meth2) = colnames(meth)
	colnames(cov2) = colnames(cov)

	return(list(pos = pos2, meth = meth2, cov = cov2))
}

