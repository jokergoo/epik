library(methods)
library(GetoptLong)

chr = "chr21"
GetoptLong("chr=s", "chr")

source("~/project/development/epik/R/correlated_regions.R")
source("~/project/development/epik/R/common_utils.R")
library(circlize)
library(matrixStats)
library(Rcpp)
sourceCpp("~/project/development/epik/src/extract_sites.cpp")

source("~/project/development/epik/roadmap/data_config.R")

cr = correlated_regions(SAMPLE_ID, EXPR, TXDB, chr = chr, window_size = 6, 
	window_step = 3, subgroup = SUBGROUP, col = SUBGROUP_COLOR)

saveRDS(cr, file = qq("@{PROJECT_DIR}/rds/cr_@{chr}.rds"))
