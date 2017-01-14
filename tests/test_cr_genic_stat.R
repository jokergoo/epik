source("test_cr_chr21.R")

sig_cr = cr[cr$corr_p < 0.05]

cr2 = reduce_cr(sig_cr, TXDB, EXPR)

cr_genic_stat(cr2, TXDB)
