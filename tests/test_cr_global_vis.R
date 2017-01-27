source("/home/guz/project/development/epik/tests/test_cr_head.R")
source("/home/guz/project/development/epik/R/correlated_regions.R")
source("/home/guz/project/development/epik/R/correlated_regions_global_vis.R")
source("/home/guz/project/development/epik/R/correlated_regions_enrich.R")
source("/home/guz/project/development/epik/R/genomic_region_correlation.R")
source("/home/guz/project/development/epik/R/genomic_region_annotation.R")
source("/home/guz/project/development/epik/R/common_utils.R")

library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(gtrellis)
library(HilbertCurve)

library(Rcpp)
sourceCpp("/home/guz/project/development/epik/src/dist_by_closeness.cpp")

cr = readRDS(qq("@{PROJECT_DIR}/rds/cr_chr21.rds"))

cr_hilbert_curve(cr)
cr_hilbert_curve(cr, legend = NULL)
cr_hilbert_curve(cr, merge_chr = FALSE)
cr_hilbert_curve(cr, add_chr_name = FALSE)

cr = cr_add_fdr_column(cr)

cr = cr_enrichedheatmap(cr, TXDB, EXPR)

sig_cr_enrichedheatmap(cr, TXDB, fdr_cutoff = 0.1, meth_diff_cutoff = 0.1)
sig_cr_compare_cutoff(cr, TXDB)

cytoband_list = gtrellis_cr_genes(cr, TXDB, EXPR)
gtrellis_sig_cytoband(cr, TXDB, cytoband_list)

cr_genes_david(cr, david_user = "z.gu@dkfz.de")
