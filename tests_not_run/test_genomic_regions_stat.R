library(circlize)
library(GenomicRanges)
library(GlobalOptions)
library(ggplot2)

source("/home/guz/project/development/epik/R/read_data_hooks.R")
source("/home/guz/project/development/epik/R/genomic_region_stat.R")

set.seed(123)
gr_list = lapply(1:10, function(i) {
	df = generateRandomBed(1000)[1:sample(100:1000, 1), ]
	GRanges(df[[1]], ranges = IRanges(df[[2]], df[[3]]))
})
names(gr_list) = paste0("sample_", 1:10)
genomic_regions_basic_stat(gr_list)
genomic_regions_basic_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
    annotation_color = c("a" = "red", "b" = "blue"))
genomic_regions_basic_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
    annotation_color = c("a" = "red", "b" = "blue"), type = "number")
genomic_regions_basic_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
    annotation_color = c("a" = "red", "b" = "blue"), type = "median_width")
genomic_regions_basic_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
    annotation_color = c("a" = "red", "b" = "blue"), by_chr = TRUE)
