source("/home/guz/project/development/epik/tests/test_cr_head.R")
source("/home/guz/project/development/epik/R/common_utils.R")
source("/home/guz/project/development/epik/R/correlated_regions.R")
source("/home/guz/project/development/epik/R/genomic_region_merge.R")
source("/home/guz/project/development/epik/R/compare_meth.R")

library(GenomicFeatures)
library(circlize)
library(Rcpp)
sourceCpp("/home/guz/project/development/epik/src/extract_sites.cpp")

library(matrixStats)

# also apply t-test on subgroup1 and subgroup2
g = genes(TXDB)
chr = structure(as.vector(seqnames(g)), names = g$gene_id)
cr = correlated_regions(sample_id, EXPR[chr[rownames(EXPR)] == "chr21", ][1:2, ], 
	TXDB, chr = "chr21", window_size = 6, window_step = 3, 
	subgroup = subgroup, col = c("group1" = "red", "group2" = "blue"))

cr = correlated_regions(sample_id, EXPR, 
	TXDB, chr = "chr21", window_size = 6, window_step = 3, 
	subgroup = subgroup, col = c("group1" = "red", "group2" = "blue"))

cr[1:2]
c(cr[1], cr[2])
cr$ncpg
cr$corr_fdr = p.adjust(cr$corr_p, "BH")
mcols(cr)

cr = cr_add_fdr_column(cr)

cr2 = cr_reduce(cr, TXDB)
cr2 = cr_reduce(cr, TXDB, EXPR)

library(gtrellis)
compare_meth(g[seqnames(g) == "chr21"], cr)

correlated_regions(sample_id, EXPR["ENSG00000197978.8", , drop = FALSE], 
	TXDB, chr = "chr15", window_size = 6, window_step = 3, 
	subgroup = subgroup)
