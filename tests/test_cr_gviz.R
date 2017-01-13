source("test_cr_chr21.R")

sig_cr = cr[cr$corr_p < 0.05]

chipseq_hooks$peak = function(mark, sid, chr) {
    qqcat("reading peaks: @{sid}, @{mark} at @{chr}\n")
    df = read.table(pipe(qq("zcat @{BASE_DIR}/data/narrow_peaks/@{sid}-@{mark}.narrowPeak.gz | grep '@{chr}'")), stringsAsFactors = FALSE)
    GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]] + 1, df[[3]]), density = df[[5]])
}

peak_list = lapply(MARKS[1:2], get_peak_list, chr = "chr21")
names(peak_list) = MARKS[1:2]

tx_list = transcriptsBy(TXDB, "gene")

for(gi in unique(sig_cr$gene_id)) {
	pdf(qq("gviz_@{gi}.pdf"), width = 8, height = 0.12*length(tx_list[[gi]]) + 8.85)
	cr_gviz(sig_cr, gi, EXPR, TXDB, gf_list = list(CGI = CGI), hm_list = peak_list)
	dev.off()
}
