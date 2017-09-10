
library(memoise)
source("/home/guz/project/development/epik/tests/test_cr_head.R")
source("/home/guz/project/development/epik/R/code_snippet.R")
source("/home/guz/project/development/epik/R/correlated_regions_enrich.R")
source("/home/guz/project/development/epik/R/correlated_regions.R")
source("/home/guz/project/development/epik/R/correlated_regions_global_vis.R")
source("/home/guz/project/development/epik/R/correlated_regions_enrichedheatmap.R")
source("/home/guz/project/development/epik/R/common_utils.R")

library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(matrixStats)
library(EnrichedHeatmap)
library(Rcpp)
sourceCpp("/home/guz/project/development/epik/src/dist_by_closeness.cpp")

cr = readRDS(qq("@{PROJECT_DIR}/rds/cr_chr21.rds"))
cr = cr_add_fdr_column(cr)

g = genes(TXDB)
tss = promoters(g, upstream = 1, downstream = 0)
tss = tss[seqnames(tss) == "chr21"]

tss = tss[unique(cr$gene_id)]

normalize_epigenomic_signals(cr, tss, MARKS, EXPR)
normalize_epigenomic_signals(cr, tss, marks = NULL)
normalize_epigenomic_signals(cr, tss, MARKS, EXPR, include_correlation_matrix = FALSE)


cr_enriched_heatmap_at_cgi(cr, TXDB, EXPR, CGI, marks = MARKS)
cr_enriched_heatmap_at_cgi(cr, TXDB, EXPR, CGI, type = "pos")


cr_enriched_heatmap_at_tss(cr, TXDB, EXPR, CGI, marks = MARKS[1])

cr = cr_enriched_heatmap(cr, TXDB, EXPR)
cr_enriched_heatmap_at_gene(cr, TXDB, EXPR, CGI, marks = MARKS[1])

tfbs = read.table("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/encode_uniform_tfbs_merged_1kb.bed", sep = "\t", stringsAsFactors = FALSE)
tfbs = GRanges(seqnames = tfbs[[1]], ranges = IRanges(tfbs[[2]], tfbs[[3]]))
cr_enriched_heatmap_at_genomic_features(cr, TXDB, EXPR, tfbs, marks = MARKS[1])

