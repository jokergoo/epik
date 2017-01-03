context("test genomic annotations")

gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))

test_that("test percentOverlaps", {
	expect_that(percentOverlaps(gr1, gr2), equals(c(0, 4/7)))
	expect_that(percentOverlaps(gr2, gr1), equals(c(0, 4/8)))
})

if(Sys.getenv("IS_PBS") != "") {

makeGRangesFromDataFrameWithFirstThreeColumns = function(df) {
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
}


BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/"
library(circlize)
gr = generateRandomBed(nr = 1000)
gr = makeGRangesFromDataFrameWithFirstThreeColumns(gr)

files = dir(qq("@{BASE_DIR}/data/narrow_peaks"), pattern = "gz$")
peak_list = lapply(sample(files, 4), function(f) {
	df = read.table(qq("@{BASE_DIR}/data/narrow_peaks/@{f}"))
	qqcat("reading @{f}\n")
	makeGRangesFromDataFrameWithFirstThreeColumns(df)
})

###
annotate_to_genomic_features(gr, peak_list[[1]])
annotate_to_genomic_features(gr, peak_list[[1]], name = "peak")
annotate_to_genomic_features(gr, peak_list[[1]], name = "peak", type = "number")

annotate_to_genomic_features(gr, peak_list)
annotate_to_genomic_features(gr, peak_list, name = c("P1", "P2", "P3", "P4"))
annotate_to_genomic_features(gr, peak_list, prefix = "")

### build a transcriptDb object
library(GenomicFeatures)
txdb = loadDb(qq("@{BASE_DIR}/data/gen10.long.sqlite"))

gr = generateRandomBed(nr = 1000)
gr = makeGRangesFromDataFrameWithFirstThreeColumns(gr)
annotate_to_gene_models(gr, txdb, gene_model = "tx")
annotate_to_gene_models(gr, txdb, gene_model = "gene")
annotate_to_gene_models(gr, txdb, gene_model = "gene", annotation_type = "number")
annotate_to_gene_models(gr, txdb, gene_model = "gene", annotation_prefix = "")

}