pos = c(1, 2, 6, 7, 8, 12, 15)

meth = cbind(c(0.9, 0.9, 0.4, 0.5, 0.5, 0.7, 0.8),
	         c(0.8, 0.8, 0.5, 0.7, 0.7, 0.2, 0.4))

cov = cbind(sample(20, 7), sample(20, 7))


cpg_dinucleotide_methylation(pos, meth, cov)

