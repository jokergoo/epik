
# use HMM to capture regions with high difference
library(STAN)
library(depmixS4)
# gr should be sorted
segment_by_hmm = function(gr, x, nStates = 5, method = c("depmixS4", "STAN")) {
	# maybe all chromosomes can be merged as a single huge chromosome
	od = order(gr)
	gr = gr[od]
	x = x[od]
	q = quantile(abs(x), 0.99)
	x[x > q] = q
	x[x < -q] = -q

	method = match.arg(method)[1]

	if(method == "STAN") {
		mat_list = list(matrix(x, ncol = 1))
		hmm = initHMM(mat_list, nStates = nStates, "Gaussian")
		hmm_fitted = fitHMM(mat_list, hmm, maxIters = 10, verbose = TRUE)
		viterbi = getViterbi(hmm_fitted, mat_list)
		states = viterbi[[1]]
	} else {
		df = data.frame(foo = x)
		mod = depmix(response = foo ~ 1, data = df, nstates = nStates)
		fm = fit(mod)
		states = posterior(fm)[, 1]
	}

	mean = tapply(x, states, mean)
	s = structure(rank(mean), names = names(mean))
	states2 = s[as.character(states)]
	gr$states = states2

	return(gr)
}


average_in_window = function(gr1, gr2, v, empty_value = 0, chromosome = unique(as.vector(seqnames(gr1)))) {
	x = rep(empty_value, length(gr1))
	chromosome = intersect(chromosome, intersect( unique(as.vector(seqnames(gr1))),  unique(as.vector(seqnames(gr2)))))
	for(chr in chromosome) {
		l = as.vector(seqnames(gr1) == chr)
		if(sum(l)) {
			window = ranges(gr1[seqnames(gr1) == chr])
			ir = ranges(gr2[seqnames(gr2) == chr])

			mtch = as.matrix(findOverlaps(window, ir))
			v2 = HilbertCurve:::average_in_window(window, ir, mtch, v, "w0", empty_value)
			x[as.vector(seqnames(gr1) == chr)][unique(mtch[, 1])] = v2
		}
	}
	return(x)
}
