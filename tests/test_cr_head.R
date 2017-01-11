source("test_head.R")

cat("load txdb...\n")
TXDB = loadDb(qq("@{BASE_DIR}/data/gen10_long_protein_coding_gene_adjusted.sqlite"))
genes = genes(TXDB)
map = structure(names(genes), names = gsub("\\.\\d+$", "", names(genes)))

####################################################
## expression data

cat("load expression...\n")
count = as.matrix(read.table(qq("@{BASE_DIR}/data/expression/57epigenomes.N.pc.gz"), row.names = 1, header = TRUE))
rpkm = as.matrix(read.table(qq("@{BASE_DIR}/data/expression/57epigenomes.RPKM.pc.gz"), row.names = 1, header = TRUE))
rownames(count) = map[rownames(count)]
rownames(rpkm) = map[rownames(rpkm)]
count = count[, sample_id]
rpkm = rpkm[, sample_id]

######################################################
## genes should have raw count > 0 in at least half samples
l = apply(count, 1, function(x) sum(x > 0) > length(x)/2)
expr = rpkm[l, , drop = FALSE]
EXPR = log2(expr + 1)   # log2(rpkm + 1)
EXPR = EXPR[intersect(rownames(EXPR), names(genes)), , drop = FALSE]
