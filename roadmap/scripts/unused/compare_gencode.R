
library(GenomicFeatures)
library(GetoptLong)

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
OUTPUT_DIR = BASE_DIR

cat("load txdb...\n")
TXDB = loadDb(qq("@{BASE_DIR}/data/gen10_long_protein_coding_gene_adjusted.sqlite"))
genes = genes(TXDB)
names(genes) = gsub("\\.\\d+$", "", names(genes))

TXDB_gencode19 = loadDb(qq("@{BASE_DIR}/data/gencode.v19.annotation_protein_coding.gtf.sqlite"))
genes_gencode19 = genes(TXDB_gencode19)
names(genes_gencode19) = gsub("\\.\\d+$", "", names(genes_gencode19))

cn = intersect(names(genes), names(genes_gencode19))
w1 = width(genes[cn]); names(w1) = cn
w2 = width(genes_gencode19[cn]); names(w2) = cn
# 14902 w1 == w2

s1 = start(promoters(genes[cn], upstream = 1, downstream = 0))
s2 = start(promoters(genes_gencode19[cn], upstream = 1, downstream = 0))
# 16758 s1 == s2

# 14895 s1 == s2 & w1 == w2

r = abs(w2 - w1)/pmin(w1, w2)
l = s1 == s2 & w1 == w2

cn = cn[l]



tx_list = as.list(transcriptsBy(TXDB, by = "gene"))
tx_list = lapply(tx_list, as.data.frame)
names(tx_list) = gsub("\\.\\d+$", "", names(tx_list))
tx_list = tx_list[cn]

tx_list2 = as.list(transcriptsBy(TXDB_gencode19, by = "gene"))
tx_list2 = lapply(tx_list2, as.data.frame)
names(tx_list2) = gsub("\\.\\d+$", "", names(tx_list2))
tx_list2 = tx_list2[cn]

g1 = lapply(tx_list, function(gr) {
	l = grepl("ENSG", gr$tx_name)
	gr = gr[!l, ]
	m = as.matrix(gr[, 2:3])
	data.frame(chr = gr[1, 1], start = min(m), end = max(m), strand = gr[1, 5], stringsAsFactors = FALSE)
})
g1 = do.call("rbind", g1)
g1 = GRanges(seqnames = g1[[1]], ranges = IRanges(g1[[2]], g1[[3]]), strand = g1[[4]], names = rownames(g1))
names(g1) = g1$names

g2 = lapply(tx_list2, function(gr) {
	l = grepl("ENSG", gr$tx_name)
	gr = gr[!l, ]
	m = as.matrix(gr[, 2:3])
	data.frame(chr = gr[1, 1], start = min(m), end = max(m), strand = gr[1, 5], stringsAsFactors = FALSE)
})
g2 = do.call("rbind", g2)
g2 = GRanges(seqnames = g2[[1]], ranges = IRanges(g2[[2]], g2[[3]]), strand = g2[[4]], names = rownames(g2))
names(g2) = g2$names

s1 = start(promoters(g1, upstream = 1, downstream = 0))
s2 = start(promoters(genes[names(g1)], upstream = 1, downstream = 0))

pdf(qq("@{OUTPUT_DIR}/plots/compare_gencode.pdf"))
plot(w1, w2, pch = 16, cex = 0.2, col = "#00000080", xlab = "gen10", ylab = "gen19", main = qq("gene width, @{length(cn)} in common"))
t1 = table(seqnames(genes[cn]))
t2 = table(seqnames(genes))
t1 = t1[paste0("chr", 1:22)]
t2 = t2[paste0("chr", 1:22)]
foo = barplot(t1, horiz = TRUE, names = rep("", 22), main = "#gene"); text(0, foo, paste0("chr", 1:22), adj = c(0, 0.5))
foo = barplot(t1/t2, horiz = TRUE, names = rep("", 22), main = "#gene/#all"); text(0, foo, paste0("chr", 1:22), adj = c(0, 0.5))

stat = c("0" = sum(s1 == s2)/length(s1),
	     "<=100" = sum(abs(s1 - s2) <= 100 & s1 != s2)/length(s1),
	     "<=500" = sum(abs(s1 - s2) <= 500 & abs(s1 - s2) > 100)/length(s1),
	     "<= 1000" = sum(abs(s1 - s2) <= 1000 & abs(s1 - s2) > 500)/length(s1),
	     "<= 5000" = sum(abs(s1 - s2) <= 5000 & abs(s1 - s2) > 1000)/length(s1),
	     ">5000" = sum(abs(s1 - s2) > 5000)/length(s1))
barplot(stat, main = "#genes with dist(gene_tss, first_tx_tss)")
dev.off()

