context("test gr_ss")

if(Sys.getenv("IS_PBS") != "") {

makeGRangesFromDataFrameWithFirstThreeColumns = function(df) {
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
}

files = dir("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/narrow_peaks/", pattern = "gz$")

peak_list = lapply(files[1:36], function(f) {
	df = read.table(paste0("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/narrow_peaks/", f))
	makeGRangesFromDataFrameWithFirstThreeColumns(df)
})

cr = common_regions(peak_list)
cr = common_regions(peak_list, min_width = 10000, min_coverage = 4)
cr = common_regions(peak_list, min_width = 1000)

factors = rep(c("T1", "T2", "T3"), each = 12)
res = subgroup_specific_genomic_regions(cr, factors = factors)
res = subgroup_specific_genomic_regions(cr, factors = factors, type = "001")
res = subgroup_specific_genomic_regions(cr, factors = factors, present = 0.8, absent = 0.2)
res = subgroup_specific_genomic_regions(cr, factors = factors, present = function(x) sum(x > 0)/length(x) > 0.5, absent = function(x) sum(x > 0)/length(x) < 0.3)

plot_subgroup_specificity_heatmap(res)

load("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/gr_list_1.RData")
genomic_features = lapply(gr_list_1[1:2], makeGRangesFromDataFrameWithFirstThreeColumns)

plot_subgroup_specificity_heatmap(res, genomic_features = genomic_features)

}
