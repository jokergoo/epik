library(methods)
suppressPackageStartupMessages(library(GetoptLong))

chr = "chr21"
GetoptLong("chr=s", "chromosome")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

compare_meth = function(cr, chr, start, end, x = NULL, x2 = NULL, main = NULL) {

	sample_id = attr(cr, "sample_id")
	factor = attr(cr, "factor")
	col = attr(cr, "col")
	n_sample = length(sample_id)
	if(is.null(factor)) factor = 1:n_sample
	if(is.null(col)) col = rand_color(n_sample)

	methylation_hooks$set(chr)

	ind = extract_sites(start, end, methylation_hooks$site(), TRUE, 0)
	meth = methylation_hooks$meth(row_index = ind, col_index = sample_id)
	raw = methylation_hooks$raw(row_index = ind, col_index = sample_id)
	cov = methylation_hooks$coverage(row_index = ind, col_index = sample_id)
	site = methylation_hooks$site(index = ind)

	gr_input = GRanges(seqnames = chr, ranges = IRanges(start, end))
	mtch = as.matrix(findOverlaps(gr_input, cr))
	cr_intersect = cr[unique(mtch[, 2])]
	s = start(cr_intersect)
	e = end(cr_intersect)
	par(mfrow = c(5 + (!is.null(x)) + (!is.null(x2)), 1), mar = c(1, 4, 1, 1))

	# smoothed methylation
	matplot(site, meth, type = "l", col = col[factor], lty = 1,
	    ylab = "smoothed meth", xlab = NULL, main = main, ylim = c(0, 1))
	for(i in seq_along(s)) {
		rect(s[i], 0, e[i], 1, col = "#FF000010", border = NA)
	}
	# legend("bottomleft", lty = 1, col = col, legend = names(col))

	# if(!is.null(x)) plot((start(x)+end(x))/2, x$corr, xlim = range(site), ylab = "corr\n(from smoothed)", xlab = NULL, type = "l")

	# raw methylation
	matplot(site, raw, type = "l", col = col[factor], lty = 1,
	    ylab = "raw meth", xlab = NULL, ylim = c(0, 1))
	for(i in seq_along(s)) {
		rect(s[i], 0, e[i], 1, col = "#FF000010", border = NA)
	}
	# if(!is.null(x2)) plot((start(x2)+end(x2))/2, x2$corr, xlim = range(site), ylab = "corr\n(from raw)", xlab = NULL, type = "l")

	# raw methylation with CpG coverage > q25
	cov_cutoff_1 = round(quantile(cov, 0.25))
	plot(site, xlim = range(site), ylim = c(0, 1), type = "n",
	    ylab = qq("raw meth (cov >= @{cov_cutoff_1})(q25)"), xlab = NULL)
	for(i in seq_len(ncol(raw))) {
	    xx = raw[, i]
	    yy = cov[, i]
	    lines(site[yy >= cov_cutoff_1], xx[yy >= cov_cutoff_1], col = col[factor[i]])
	}

	# raw methylation with CpG coverage > q50
	cov_cutoff_2 = round(quantile(cov, 0.5))
	plot(site, xlim = range(site), ylim = c(0, 1), type = "n",
	    ylab = qq("raw meth (cov >= @{cov_cutoff_2})(q50)"), xlab = NULL)
	for(i in seq_len(ncol(raw))) {
	    xx = raw[, i]
	    yy = cov[, i]
	    lines(site[yy >= cov_cutoff_2], xx[yy >= cov_cutoff_2], col = col[factor[i]])
	}
	par(mar = c(4, 4, 1, 1))
	plot(site, rowMeans(cov), type = 'h', ylab = "CpG coverage", xlab = "CpG sites")

}

neg_cr_smooth = readRDS(qq("@{OUTPUT_DIR}/rds/all_neg_cr_w6s3_fdr_less_than_0.05_methdiff_larger_than_0.1.rds"))
pos_cr_smooth = readRDS(qq("@{OUTPUT_DIR}/rds/all_pos_cr_w6s3_fdr_less_than_0.05_methdiff_larger_than_0.1.rds"))

# merge CRs if they are too close to each other
gr1 = reduce(neg_cr_smooth, min.gapwidth = 1000)
gr2 = reduce(pos_cr_smooth, min.gapwidth = 1000)

# subset to the current chromosome
gr1 = gr1[seqnames(gr1) == chr]
gr2 = gr2[seqnames(gr2) == chr]

# for each region in `gr1` and `gr2`, compare methylation distribution of smoothed data and raw data
start = start(gr1)
end = end(gr1)
pdf(qq("@{OUTPUT_DIR}/plots/neg_cr_meth_compare_@{chr}.pdf"), width = 8, height = 10)
for(i in seq_len(length(gr1))) {
	compare_meth(neg_cr_smooth, chr, start[i], end[i], main = qq("neg_cr, @{chr}:@{start[i]}-@{end[i]}"))
}
dev.off()

start = start(gr2)
end = end(gr2)
pdf(qq("@{OUTPUT_DIR}/plots/pos_cr_meth_compare_@{chr}.pdf"), width = 8, height = 10)
for(i in seq_len(length(gr2))) {
	compare_meth(pos_cr_smooth, chr, start[i], end[i], main = qq("pos_cr, @{chr}:@{start[i]}-@{end[i]}"))
}
dev.off()


# for(chr in paste0("chr", 1:22)) {
#     cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/compare_smooth.R --chr @{chr}")    
#     cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=15G -N compare_smooth_@{chr}' '@{cmd}'")
#     system(cmd)
# }
