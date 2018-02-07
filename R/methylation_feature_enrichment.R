
cg_percent_per_window = function(w, chromosome = paste0("chr", 1:22), genome = "hg19") {

	chromInfo = getChromInfoFromUCSC(genome)
	chromInfo = chromInfo[chromInfo$chrom %in% chromosome, ]
	chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
	chromGr_window = makeWindows(chromGr, w = w, short.keep = TRUE)

	chromGr_window_lt = as.list(split(chromGr_window, as.vector(seqnames(chromGr_window))))
	chromGr_window_lt = lapply(names(chromGr_window_lt), function(chr) {
		methylation_hooks$set_chr(chr)
		cpg = methylation_hooks$gr
		gr = chromGr_window_lt[[chr]]
		gr$cpg_per_1kbwindow = countOverlaps(gr, cpg)/window*w
		gr
	})

	chromGr_window = unlist(GRangesList(chromGr_window_lt))
	return(chromGr_window)
}

window_env = new.env()

generate_background_from_methylation_features = function(gr, test_chr = NULL, q = c(0.05, 0.95)) {
	if(inherits(gr, "list")) {
		gr = unlist(GRangesList(gr))
	}
	if(is.null(test_chr)) {
		w = sapply(split(gr, as.vector(seqnames(gr))), function(x) sum(as.numeric(width(x))))
		test_chr = names(w[which.max(w)])
	}
	methylation_hooks$set_chr(test_chr)
	cpg = methylation_hooks$gr
	gr = gr[seqnames(gr) == test_chr]
	w = round(quantile(width(gr), 0.25), digits = -3)
	w[w > 10000] = 10000
	w[w < 1000] = 1000

	qqcat("window size is set to @{w}\n")

	gr2 = makeWindows(gr, w = w, short.keep = TRUE)
	gr2 = gr2[width(gr2) > w/4]

	gr2$cpg_per_1kbwindow = countOverlaps(gr2, cpg)/width(gr2)*1000
	x = gr2$cpg_per_1kbwindow
	x = x[x > 0]
	qa = quantile(x, q)

	window = readRDS(qq("@{PROJECT_DIR}/rds/cg_percent_per_window_@{w}.rds"))
	l = window$cpg_per_1kbwindow > qa[1] & window$cpg_per_1kbwindow < qa[2]
	qqcat("@{sum(l)}/@{length(l)} were kept.\n")
	reduce(sort(window[l]))
}


