library(GenomicRanges)
library(GetoptLong)
library(GenomicFeatures)
library(matrixStats)

source("/home/guz/project/development/epik/R/genomic_region_correlation.R")
source("/home/guz/project/development/epik/R/systemdf.R")
source("/home/guz/project/development/epik/R/genomic_region_annotation.R")

gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))
background = GRanges(seqnames = "chr1", ranges = IRanges(start = 7, end = 15))

test_that("test genomic correlation functions", {

	expect_that(genomic_corr_jaccard(gr1, gr2), equals(4/16))
	expect_that(genomic_corr_jaccard(gr1, gr2, background = background), equals(3/5))

	expect_that(genomic_corr_absdist(gr1, gr2), equals(2.5))

	expect_that(genomic_corr_intersect(gr1, gr2, method = "number"), equals(1))
	expect_that(genomic_corr_intersect(gr1, gr2, method = "percent"), equals(0.4))
	expect_that(genomic_corr_intersect(gr1, gr2, method = "length"), equals(4))
})


makeGRangesFromDataFrameWithFirstThreeColumns = function(df) {
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
}

PROJECT_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/epik_roadmap"

files = dir(qq("@{PROJECT_DIR}/data/narrow_peaks"), pattern = "gz$")
peak_list = lapply(sample(files, 4), function(f) {
	df = read.table(qq("@{PROJECT_DIR}/data/narrow_peaks/@{f}"))
	qqcat("reading @{f}\n")
	makeGRangesFromDataFrameWithFirstThreeColumns(df)
})
names(peak_list) = paste0("peak_", seq_along(peak_list))

txdb = loadDb(qq("@{PROJECT_DIR}/data/gen10.long.sqlite"))
genomic_features = list(
	gene = genes(txdb),
	exon = exons(txdb)
)
strand(genomic_features[[1]]) = "*"
strand(genomic_features[[2]]) = "*"

# general test

gr1 = peak_list[[1]]
gr2 = genomic_features[[1]]
genomic_corr_reldist(gr1, gr2)
genomic_corr_jaccard(gr1, gr2)
genomic_corr_absdist(gr1, gr2)
genomic_corr_absdist(gr1, gr2, method = median)
genomic_corr_absdist(gr1, gr2, trim = 0.1)
genomic_corr_intersect(gr1, gr2, method = "number")
genomic_corr_intersect(gr1, gr2, method = "percent")
genomic_corr_intersect(gr1, gr2, method = "length")

# test if some chromosomes donot exist in the gr2
gr1 = gr1[seqnames(gr1) %in% c("chr1", "chr2")]
gr2 = gr2[seqnames(gr2) %in% c("chr2", "chr3")]
genomic_corr_reldist(gr1, gr2)
genomic_corr_jaccard(gr1, gr2)
genomic_corr_absdist(gr1, gr2)
genomic_corr_intersect(gr1, gr2, method = "number")
genomic_corr_intersect(gr1, gr2, method = "percent")
genomic_corr_intersect(gr1, gr2, method = "length")

### test correlation main function
genomic_regions_correlation(peak_list, genomic_features, nperm = 4)
genomic_regions_correlation(peak_list, genomic_features, nperm = 4, chromosome = c("chr1", "chr2"))
genomic_regions_correlation(peak_list, genomic_features, nperm = 4, mc.cores = 2)
genomic_regions_correlation(peak_list, genomic_features, nperm = 4, stat_fun = genomic_corr_reldist)
genomic_regions_correlation(peak_list, genomic_features, nperm = 4, stat_fun = genomic_corr_absdist, trim = 0.1)

### with background
background = systemdf("bedtools random -l 1000 -n 10000 -g /icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/bed/hg19.len | sort -V -k1,1 -k2,2 | bedtools merge -i stdin")
background = makeGRangesFromDataFrameWithFirstThreeColumns(background)

# test with restrict set
genomic_corr_jaccard(gr1, gr2, background = background)
genomic_corr_intersect(gr1, gr2, background = background, method = "number")

genomic_regions_correlation(peak_list, genomic_features, nperm = 4, background = background)
