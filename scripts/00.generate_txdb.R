library(methods)
suppressPackageStartupMessages(library(GetoptLong))

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
OUTPUT_DIR = BASE_DIR

library(GenomicFeatures)

library(GTF)
gen10 = new("GTF")
gen10$read(qq("@{BASE_DIR}/data/gen10.long.gtf"))
l = sapply(gen10$gtf, function(g) g$type == "protein_coding")
gen10$gtf = gen10$gtf[l]
gen10$gtf = lapply(gen10$gtf, function(g) {
	tx_list = g$transcript
	l = sapply(tx_list, function(tx) tx$type == "protein_coding")
	if(sum(l) == 0) return(NULL)
	ttx_list = tx_list[l]
	s = sapply(tx_list, function(tx) tx$start)
	e = sapply(tx_list, function(tx) tx$end)
	g$transcript = tx_list
	g$start = min(s)
	g$end = max(e)
	g
})
gen10$gtf = gen10$gtf[!sapply(gen10$gtf, is.null)]

map = structure(names = gsub("\\.\\d+$", "", names(gen10$gtf)), names(gen10$gtf))
names(gen10$gtf) = gsub("\\.\\d+$", "", names(gen10$gtf))
chr = sapply(gen10$gtf, function(x) x$chr)

gen19 = new("GTF")
gen19$read(qq("@{BASE_DIR}/data/gencode.v19.annotation.gtf"))
l = sapply(gen19$gtf, function(g) g$type == "protein_coding")
gen19$gtf = gen19$gtf[l]
gen19$gtf = lapply(gen19$gtf, function(g) {
	tx_list = g$transcript
	l = sapply(tx_list, function(tx) tx$type == "protein_coding")
	if(sum(l) == 0) return(NULL)
	tx_list = tx_list[l]
	s = sapply(tx_list, function(tx) tx$start)
	e = sapply(tx_list, function(tx) tx$end)
	g$transcript = tx_list
	g$start = min(s)
	g$end = max(e)
	g
})
gen19$gtf = gen19$gtf[!sapply(gen19$gtf, is.null)]

names(gen19$gtf) = gsub("\\.\\d+$", "", names(gen19$gtf))

cn = intersect(names(gen10$gtf), names(gen19$gtf))
qqcat("filter by name: @{length(gen10$gtf)} -> @{length(cn)} <- @{length(gen19$gtf)}\n")
gen10$gtf = gen10$gtf[cn]
gen19$gtf = gen19$gtf[cn]

s10 = sapply(gen10$gtf, function(g) g$start)
s19 = sapply(gen19$gtf, function(g) g$start)
w10 = sapply(gen10$gtf, function(g) g$end - g$start + 1)
w19 = sapply(gen19$gtf, function(g) g$end - g$start + 1)

l = s10 == s19 & w10 == w19
qqcat("filter by tss and width: @{sum(l)}/@{length(l)}\n")
gen10$gtf = gen10$gtf[l]
names(gen10$gtf) = map[names(gen10$gtf)]

gen10$write(qq("@{OUTPUT_DIR}/data/gen10_long_protein_coding_gene_adjusted.gtf"))

txdb = makeTranscriptDbFromGFF(qq("@{OUTPUT_DIR}/data/gen10_long_protein_coding_gene_adjusted.gtf"),
    format = "gtf", gffGeneIdAttributeName = "gene_id")
saveDb(txdb, file = qq("@{OUTPUT_DIR}/data/gen10_long_protein_coding_gene_adjusted.sqlite"))


pdf(qq("@{OUTPUT_DIR}/plots/compare_gencode.pdf"), width = 8, height = 5)
par(mfrow = c(1, 2))
t1 = table(chr[cn[l]])
t2 = table(chr)
t1 = t1[paste0("chr", 1:22)]
t2 = t2[paste0("chr", 1:22)]
foo = barplot(t1, horiz = TRUE, names = rep("", 22), main = "#gene")
text(0, foo, paste0("chr", 1:22), adj = c(0, 0.5), cex = 0.6)
foo = barplot(t1/t2, horiz = TRUE, names = rep("", 22), main = "#gene/#all")
text(0, foo, paste0("chr", 1:22), adj = c(0, 0.5), cex = 0.6)
dev.off()

