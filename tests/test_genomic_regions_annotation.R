library(GenomicFeatures)
library(GetoptLong)
library(circlize)

source("/home/guz/project/development/epik/R/genomic_region_annotation.R")


gr1 = GRanges(seqname = "chr1", ranges = IRanges(start = c(4, 10), end = c(6, 16)))
gr2 = GRanges(seqname = "chr1", ranges = IRanges(start = c(7, 13), end = c(8, 20)))

percentOverlaps(gr1, gr2) #, equals(c(0, 4/7)))
percentOverlaps(gr2, gr1) #, equals(c(0, 4/8)))


makeGRangesFromDataFrameWithFirstThreeColumns = function(df) {
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
}


PROJECT_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/epik_roadmap"
gr = generateRandomBed(nr = 1000)
gr = makeGRangesFromDataFrameWithFirstThreeColumns(gr)

files = dir(qq("@{PROJECT_DIR}/data/narrow_peaks"), pattern = "gz$")
peak_list = lapply(sample(files, 4), function(f) {
	df = read.table(qq("@{PROJECT_DIR}/data/narrow_peaks/@{f}"))
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
txdb = loadDb(qq("@{PROJECT_DIR}/data/gen10.long.sqlite"))

gr = generateRandomBed(nr = 1000)
gr = makeGRangesFromDataFrameWithFirstThreeColumns(gr)
annotate_to_gene_models(gr, txdb, gene_model = "tx")
annotate_to_gene_models(gr, txdb, gene_model = "gene")
annotate_to_gene_models(gr, txdb, gene_model = "gene", annotation_type = "number")
annotate_to_gene_models(gr, txdb, gene_model = "gene", annotation_prefix = "")
