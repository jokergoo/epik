source("test_cr_head.R")

cr = correlated_regions(sample_id, EXPR, 
	TXDB, chr = "chr21", window_size = 6, window_step = 3, 
	subgroup = subgroup)

cr_hilbert_curve(cr)
cr_hilbert_curve(cr, legend = NULL)
cr_hilbert_curve(cr, merge_chr = FALSE)
cr_hilbert_curve(cr, add_chr_name = FALSE)

cr = add_fdr_column(cr)

cr = cr_enrichedheatmap(cr, TXDB, expr)

sig_cr_enrichedheatmap(cr, TXDB, fdr_cutoff = 0.1, meth_diff_cutoff = 0.1)
sig_cr_compare_cutoff(cr, TXDB)

cytoband_list = gtrellis_cr_genes(cr, TXDB, EXPR)
gtrellis_sig_cytoband(cr, TXDB, cytoband_list)
cr_genes_david(cr, david_user = "z.gu@dkfz.de")
