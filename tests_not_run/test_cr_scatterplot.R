source("/home/guz/project/development/epik/tests/test_cr_head.R")
source("/home/guz/project/development/epik/R/correlated_regions_scatterplot.R")

cr = readRDS(qq("@{PROJECT_DIR}/rds/cr_chr21.rds"))

cr_scatterplot(cr[1], EXPR)
