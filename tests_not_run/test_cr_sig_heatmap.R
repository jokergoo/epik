source("/home/guz/project/development/epik/tests/test_cr_head.R")
source("/home/guz/project/development/epik/R/correlated_regions_sig_heatmap.R")
source("/home/guz/project/development/epik/R/correlated_regions.R")
source("/home/guz/project/development/epik/R/genomic_region_merge.R")
source("/home/guz/project/development/epik/R/genomic_region_annotation.R")
source("/home/guz/project/development/epik/R/common_utils.R")

library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(matrixStats)

cr = readRDS(qq("@{PROJECT_DIR}/rds/cr_chr21.rds"))
cr = cr_add_fdr_column(cr)
sig_cr = cr[cr$corr_fdr < 0.1]
sig_cr = reduce_cr(sig_cr, TXDB, EXPR)
sig_cr_heatmap(sig_cr, TXDB, EXPR, gf_list = list(CGI = CGI, shore = CGI_SHORE))
