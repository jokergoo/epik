source("/home/guz/project/development/epik/tests/test_cr_head.R")

library(Gviz)
source("/home/guz/project/development/epik/R/correlated_regions_gviz.R")
source("/home/guz/project/development/epik/R/correlated_regions.R")
library(circlize)
library(Rcpp)
sourceCpp("/home/guz/project/development/epik/src/extract_sites.cpp")

library(GenomicFeatures)
library(RColorBrewer)

cr = readRDS(qq("@{PROJECT_DIR}/rds/cr_chr21.rds"))

sig_cr = cr[cr$corr_p < 0.05]

peak_list = lapply(MARKS, get_peak_list, chr = "chr21")
names(peak_list) = MARKS

tx_list = transcriptsBy(TXDB, "gene")

for(gi in unique(sig_cr$gene_id)[1:6]) {
	pdf(qq("gviz_@{gi}.pdf"), width = 8, height = 0.12*length(tx_list[[gi]]) + 8.85)
	cr_gviz(sig_cr, gi, EXPR, TXDB, gf_list = list(CGI = CGI), hm_list = peak_list)
	dev.off()
}
