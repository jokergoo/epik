source("~/project/development/epik/tests/test_cr_head.R")
source("~/project/development/epik/R/correlated_regions_genic_stat.R")
source("~/project/development/epik/R/correlated_regions.R")
source("~/project/development/epik/R/genomic_region_merge.R")
source("~/project/development/epik/R/common_utils.R")

library(matrixStats)

cr = readRDS(qq("@{PROJECT_DIR}/rds/cr_chr21.rds"))

sig_cr = cr[cr$corr_p < 0.05]

cr2 = cr_reduce(sig_cr, TXDB, EXPR)

cr_genic_stat(cr2, TXDB)
