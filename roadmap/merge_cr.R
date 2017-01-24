library(methods)
library(GetoptLong)

library(GlobalOptions)
library(GenomicRanges)
source("~/project/development/epik/R/read_data_hooks.R")
library(GenomicFeatures)
source("~/project/development/epik/R/common_utils.R")
source("~/project/development/epik/R/correlated_regions.R")
library(circlize)

source("~/project/development/epik/roadmap/data_config.R")


cr = GRanges()
for(chr in CHROMOSOME) {
	cr = c(cr, readRDS(qq("@{PROJECT_DIR}/rds/cr_@{chr}.rds")))
}
cr = cr_add_fdr_column(cr)
saveRDS(cr, file = qq("@{PROJECT_DIR}/rds/cr_all.rds"))
