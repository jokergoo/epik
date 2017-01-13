source("test_cr_head.R")

# also apply t-test on subgroup1 and subgroup2
g = genes(TXDB)
chr = structure(as.vector(seqnames(g)), names = g$gene_id)
cr = correlated_regions(sample_id, EXPR[chr[rownames(EXPR)] == "chr21", ][1:10, ], 
	TXDB, chr = "chr21", window_size = 6, window_step = 3, 
	subgroup = subgroup)

