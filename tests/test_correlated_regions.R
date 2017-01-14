source("test_cr_head.R")

# also apply t-test on subgroup1 and subgroup2
g = genes(TXDB)
chr = structure(as.vector(seqnames(g)), names = g$gene_id)
cr = correlated_regions(sample_id, EXPR[chr[rownames(EXPR)] == "chr21", ][1:2, ], 
	TXDB, chr = "chr21", window_size = 6, window_step = 3, 
	subgroup = subgroup)
cr = correlated_regions(sample_id, EXPR, 
	TXDB, chr = "chr21", window_size = 6, window_step = 3, 
	subgroup = subgroup)

cr[1:2]
c(cr[1], cr[2])
cr$ncpg
cr$corr_fdr = p.adjust(cr$corr_p, "BH")
mcols(cr)

cr_add_fdr_column(cr)

cr2 = reduce_cr(cr, TXDB)
cr2 = reduce_cr(cr, TXDB, EXPR)

compare_meth(g[seqnames(g) == "chr21"], cr)
