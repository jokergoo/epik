require(circlize)
set.seed(123)
gr_list = lapply(1:10, function(i) {
	df = generateRandomBed(1000)[1:sample(100:1000, 1), ]
	GRanges(df[[1]], ranges = IRanges(df[[2]], df[[3]]))
})
names(gr_list) = paste0("sample_", 1:10)
basic_genomic_regions_stat(gr_list)
basic_genomic_regions_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
    annotation_color = c("a" = "red", "b" = "blue"))
basic_genomic_regions_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
    annotation_color = c("a" = "red", "b" = "blue"), type = "number")
basic_genomic_regions_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
    annotation_color = c("a" = "red", "b" = "blue"), type = "median_width")
basic_genomic_regions_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
    annotation_color = c("a" = "red", "b" = "blue"), by_chr = TRUE)
