
suppressPackageStartupMessages(library(GetoptLong))
cutoff = 0.05
GetoptLong("cutoff=f", "cutoff for filter cr")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))


neg_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3_fdr_less_than_0.05_methIQR_larger_than_0.2.rds")))
pos_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3_fdr_less_than_0.05_methIQR_larger_than_0.2.rds")))

n1 = sum(neg_cr$n)
n2 = sum(pos_cr$n)
w1 = sum(width(neg_cr))
w2 = sum(width(pos_cr))

pdf(qq("@{OUTPUT_DIR}/plots/cr_number.pdf"), width = 4, height = 6)
par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))
barplot(c("neg_cr" = n1, "pos_cr" = n2), beside = TRUE, col = c("green", "red"), 
    ylab = "#CpG", axes = FALSE)
axis(side = 2, at = c(0, 2, 4, 6, 8)*100000, labels = c("0K", "200K", "400K", "600K", "800K"))
barplot(c("neg_cr" = w1, "pos_cr" = w2), beside = TRUE, col = c("green", "red"), 
    ylab = "sum(width(cr))", axes = FALSE)
axis(side = 2, at = c(0, 10, 20, 30)*1000000, labels = c("0MB", "10MB", "20MB", "30MB"))
dev.off()


cr = c(neg_cr, pos_cr)
gm = genes(TXDB)
gm = gm[seqnames(gm) %in% CHROMOSOME]
g = gm[gm$gene_id %in% unique(cr$gene_id)]
strand(g) = "*"
start(g) = start(g) - 500000
end(g) = end(g) + 50000
start(g) = ifelse(start(g) > 1, start(g), 1)

r1 = epic::genomic_regions_correlation(list(neg_cr = neg_cr, pos_cr = pos_cr), list(cgi = CGI, shore = CGI_SHORE), 
	nperm = 1000, background = g, mc.cores = 2)
r2 = epic::genomic_regions_correlation(list(neg_cr = neg_cr, pos_cr = pos_cr), list(cgi = CGI, shore = CGI_SHORE), 
	nperm = 1000, mc.cores = 2)

lt = list(r1 = r1, r2 = r2)
saveRDS(lt, file = qq(paste0(OUTPUT_DIR, "/rds/cr_enrich_to_cgi.rds")))


# pdf(qq("@{OUTPUT_DIR}/hilbert_sig.pdf"), width = 8, height = 8)
# cr_hilbert(cr = cr_filtered, txdb = TXDB, merge_chr = TRUE)
# cr_hilbert(cr = cr_filtered, txdb = TXDB, merge_chr = FALSE)
# dev.off()


neg_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3.rds")))
pos_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3.rds")))
cr = c(neg_cr, pos_cr)

# pdf(qq("@{OUTPUT_DIR}/hilbert_all.pdf"), width = 14, height = 12)
# cr_hilbert(cr, chromosome = CHROMOSOME, merge_chr = TRUE)
# cr_hilbert(cr, chromosome = CHROMOSOME, merge_chr = FALSE)
# dev.off()


## how neg_cr enriched at tss

sample_id = attr(neg_cr, "sample_id")
expr = EXPR[unique(cr$gene_id), sample_id, drop = FALSE]
km = kmeans(expr, centers = 2)$cluster
gi1 = names(km[km == 1])
gi2 = names(km[km == 2])
if(mean(expr[gi1, ]) > mean(expr[gi2, ])) {
	gi3 = gi1
	gi1 = gi2
	gi2 = gi3
}
pdf(qq("@{OUTPUT_DIR}/plots/cr_tss_fdr_0.05_iqr_0.2.pdf"), width = 10, height = 4)
cr_enriched_at_tss(c(neg_cr[neg_cr$corr_fdr < 0.05 & neg_cr$meth_IQR > 0.2], 
	pos_cr[pos_cr$corr_fdr < 0.05 & pos_cr$meth_IQR > 0.2]), TXDB)
cr_enriched_at_tss(c(neg_cr[neg_cr$corr_fdr < 0.05 & neg_cr$meth_IQR > 0.2 & neg_cr$gene_id %in% gi1], 
	pos_cr[pos_cr$corr_fdr < 0.05 & pos_cr$meth_IQR > 0.2 & pos_cr$gene_id %in% gi1]), TXDB, main = qq("mean_expr = @{mean(expr[gi1, ])}"))
cr_enriched_at_tss(c(neg_cr[neg_cr$corr_fdr < 0.05 & neg_cr$meth_IQR > 0.2 & neg_cr$gene_id %in% gi2], 
	pos_cr[pos_cr$corr_fdr < 0.05 & pos_cr$meth_IQR > 0.2 & pos_cr$gene_id %in% gi2]), TXDB, main = qq("mean_expr = @{mean(expr[gi2, ])}"))
dev.off()

pdf(qq("@{OUTPUT_DIR}/plots/cr_genebody_fdr_0.05_iqr_0.2.pdf"), width = 10, height = 4)
cr_enriched_at_gene_body(c(neg_cr[neg_cr$corr_fdr < 0.05 & neg_cr$meth_IQR > 0.2], 
	pos_cr[pos_cr$corr_fdr < 0.05 & pos_cr$meth_IQR > 0.2]), TXDB, extend = 50000)
cr_enriched_at_gene_body(c(neg_cr[neg_cr$corr_fdr < 0.05 & neg_cr$meth_IQR > 0.2 & neg_cr$gene_id %in% gi1], 
	pos_cr[pos_cr$corr_fdr < 0.05 & pos_cr$meth_IQR > 0.2 & pos_cr$gene_id %in% gi1]), TXDB, main = qq("mean_expr = @{mean(expr[gi1, ])}"), extend = 50000)
cr_enriched_at_gene_body(c(neg_cr[neg_cr$corr_fdr < 0.05 & neg_cr$meth_IQR > 0.2 & neg_cr$gene_id %in% gi2], 
	pos_cr[pos_cr$corr_fdr < 0.05 & pos_cr$meth_IQR > 0.2 & pos_cr$gene_id %in% gi2]), TXDB, main = qq("mean_expr = @{mean(expr[gi2, ])}"), extend = 50000)
dev.off()



pdf(qq("@{OUTPUT_DIR}/plots/hilbert_all.pdf"), width = 16, height = 12)
cr_hilbert(cr, chromosome = CHROMOSOME, merge_chr = TRUE)
cr_hilbert(cr, chromosome = CHROMOSOME, merge_chr = FALSE)
dev.off()


q(save = "no")


# cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/correlated_regions_downstream.R")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=30:00:00,mem=30G -N correlated_regions_downstream' '@{cmd}'")
# system(cmd)

############### following can be deleted ##################


gm = genes(TXDB)
gm = gm[seqnames(gm) %in% CHROMOSOME]
GENOMIC_FEATURE_LIST = list()
GENOMIC_FEATURE_LIST$gene = gm
strand(GENOMIC_FEATURE_LIST$gene) = "*"
GENOMIC_FEATURE_LIST$exon = exons(TXDB)
strand(GENOMIC_FEATURE_LIST$exon) = "*"
GENOMIC_FEATURE_LIST$intron = setdiff(GENOMIC_FEATURE_LIST$gene, GENOMIC_FEATURE_LIST$exon)
GENOMIC_FEATURE_LIST$tss = promoters(gm, upstream = 1500, downstream = 500)
strand(GENOMIC_FEATURE_LIST$tss) = "*"
intergenic = gaps(sort(GENOMIC_FEATURE_LIST$gene))
intergenic = intergenic[ strand(intergenic) == "*"]
intergenic = intergenic[seqnames(intergenic) %in% CHROMOSOME]
GENOMIC_FEATURE_LIST$intergenic = intergenic



cr = GRanges()
for(chr in CHROMOSOME) {
	qqcat("loading @{chr}...\n")
	cr_tmp = readRDS_or_readRData(qq(paste0(OUTPUT_DIR, "/rds/@{chr}_cr.rds")))
	l = !is.na(cr_tmp$corr)
	cr = c(cr, cr_tmp[l])
}

cr_list = list()
cr_list["neg_0.7"] = cr[cr$corr < -0.7]
cr_list["neg_0.5"] = cr[cr$corr < -0.5]
cr_list["neg_0.3"] = cr[cr$corr < -0.3]
cr_list["neg_0.1"] = cr[cr$corr < -0.1]
cr_list["pos_0.1"] = cr[cr$corr > 0.1]
cr_list["pos_0.3"] = cr[cr$corr > 0.3]
cr_list["pos_0.5"] = cr[cr$corr > 0.5]
cr_list["pos_0.7"] = cr[cr$corr > 0.7]

library(epic)

g = gm[gm$gene_id %in% unique(cr$gene_id)]
strand(g) = "*"
start(g) = start(g) - 500000
end(g) = end(g) + 50000
start(g) = ifelse(start(g) > 1, start(g), 1)

res = genomic_regions_correlation(cr_list, GENOMIC_FEATURE_LIST, background = g, nperm = 20, )


########################################
cr = GRanges()
for(chr in CHROMOSOME) {
	qqcat("loading @{chr}...\n")
	cr_tmp = readRDS_or_readRData(qq(paste0(OUTPUT_DIR, "/rds/@{chr}_cr_w6s3.rds")))
	l = !is.na(cr_tmp$corr)
	cr = c(cr, cr_tmp[l])
}

# separate neg and pos
neg_cr = cr[cr$corr < 0]
pos_cr = cr[cr$corr > 0]

# add corr_fdr column
neg_cr$corr_fdr = p.adjust(neg_cr$corr_p, "BH")
pos_cr$corr_fdr = p.adjust(pos_cr$corr_p, "BH")

saveRDS(neg_cr, file = qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3.rds")))
saveRDS(pos_cr, file = qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3.rds")))

cr = c(neg_cr, pos_cr)
w neg_cr enriched at tss
pdf(qq("@{OUTPUT_DIR}/plots/cr_tss.pdf"), width = 10, height = 4)
cr_enriched_at_tss(c(neg_cr[neg_cr$corr_fdr < 0.01], pos_cr[pos_cr$corr_fdr < 0.01]), TXDB)
dev.off()

pdf(qq("@{OUTPUT_DIR}/plots/cr_qc.pdf"), width = 10, height = 4)
cr_qc(neg_cr, pos_cr, TXDB, chromosome = CHROMOSOME, region = "all")
cr_qc(neg_cr, pos_cr, TXDB, chromosome = CHROMOSOME, region = "intergenic")
cr_qc(neg_cr, pos_cr, TXDB, chromosome = CHROMOSOME, region = "tss")
cr_qc(neg_cr, pos_cr, TXDB, chromosome = CHROMOSOME, region = "gene")
dev.off()

#### test how many clusters on rows
submat = mat[sample(nrow(mat), 2000), sample(ncol(mat), 100)]

library(mclust)
d_clust <- Mclust(submat, G=1:10)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
plot(d_clust)

library(cluster)
x = clusGap(submat, kmeans, 10, B = 500)
plot(x)


d1 = density(neg_cr$corr)
d2 = density(pos_cr$corr)
plot(1, type = "n", xlim = range(c(d1$x, d2$x)), ylim = range(c(d1$y, d2$y)))
lines(d1$x, d1$y, col = "green")
lines(d2$x, d2$y, col = "red")

for(cutoff in 7:1* (-0.1)) {
	neg_cr2 = reduce_cr(neg_cr[neg_cr$corr < cutoff], TXDB)
	saveRDS(neg_cr2, file = qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3_corr_less_than_@{cutoff}.rds")))
}
for(cutoff in 7:1* (0.1)) {
	pos_cr2 = reduce_cr(pos_cr[pos_cr$corr > cutoff], TXDB)
	saveRDS(pos_cr2, file = qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3_corr_larger_than_@{cutoff}.rds")))
}

for(cutoff in c(0.1, 0.05, 0.01)) {
	neg_cr2 = reduce_cr(neg_cr[neg_cr$corr_fdr < cutoff], TXDB)
	saveRDS(neg_cr2, file = qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3_fdr_less_than_@{cutoff}.rds")))
}

for(cutoff in c(0.1, 0.05, 0.01)) {
	pos_cr2 = reduce_cr(pos_cr[pos_cr$corr_fdr < cutoff], TXDB)
	saveRDS(pos_cr2, file = qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3_fdr_less_than_@{cutoff}.rds")))
}

foo = function(region) {

n_neg = matrix(0, nrow = 7, ncol = 5)
IQR_cutoff = c(0, 0.1, 0.2, 0.3, 0.4)
colnames(n_neg) = paste0("IQR > ", IQR_cutoff)
cor_cutoff = 1:7/10
rownames(n_neg) = paste(">", cor_cutoff)
n_pos = n_neg
w_neg = n_neg
w_pos = n_neg
for(cor_cutoff in 7:1*0.1) {
	for(iqr_cutoff in c(0, 0.1, 0.2, 0.3, 0.4)) {
		i = paste0("> ", cor_cutoff)
		j = paste0("IQR > ", iqr_cutoff)

		qqcat("cor @{i}, @{j}\n")

		neg_cr2 = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/rds/all_neg_cr_w6s3_corr_less_than_-@{cor_cutoff}_methIQR_larger_than_@{iqr_cutoff}.rds"))
		# neg_cr2_list = split(neg_cr2, neg_cr2$gene_id)
		gi = unique(neg_cr2$gene_id)
		gl = width(gm[gi])
		names(gl) = gi

		gene_tss_dist = split(neg_cr2$gene_tss_dist, neg_cr2$gene_id)
		start = split(start(neg_cr2), neg_cr2$gene_id)
		end = split(end(neg_cr2), neg_cr2$gene_id)

		for(k in seq_along(gene_tss_dist)) {
			# neg_cr3 = neg_cr2_list[[k]]
			# if(region == "tss") {
			# 	neg_cr3 = neg_cr3[neg_cr3$gene_tss_dist > -1000 & neg_cr3$gene_tss_dist < 2000]
			# } else if(region == "gene") {
			# 	neg_cr3 = neg_cr3[neg_cr3$gene_tss_dist > 2000 & neg_cr3$gene_tss_dist < gl[gi[k]]]
			# } else if(region == "intergenic") {
			# 	neg_cr3 = neg_cr3[neg_cr3$gene_tss_dist < -1000 | neg_cr3$gene_tss_dist > gl[gi[k]]]
			# }
			x = gene_tss_dist[[k]]
			if(region == "tss") {
				l = x > -1000 & x < 2000
			} else if(region == "gene") {
				l = x > 2000 & x < gl[gi[k]]
			} else if(region == "intergenic") {
				l = x < -1000 | x > gl[gi[k]]
			}
			n_neg[i, j] = n_neg[i, j] + sum(l)
			w_neg[i, j] = w_neg[i, j] + sum(end[[k]][l] - start[[k]][l] + 1)
			qqcat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
			qqcat("  @{k}/@{length(gene_tss_dist)} neg_cr gene...")
		}
		cat("\n")

		pos_cr2 = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/rds/all_pos_cr_w6s3_corr_larger_than_@{cor_cutoff}_methIQR_larger_than_@{iqr_cutoff}.rds"))
		# pos_cr2_list = split(pos_cr2, pos_cr2$gene_id)
		gi = unique(pos_cr2$gene_id)
		gl = width(gm[gi])
		names(gl) = gi

		gene_tss_dist = split(pos_cr2$gene_tss_dist, pos_cr2$gene_id)
		start = split(start(pos_cr2), pos_cr2$gene_id)
		end = split(end(pos_cr2), pos_cr2$gene_id)

		for(k in seq_along(gene_tss_dist)) {
			# pos_cr3 = pos_cr2_list[[k]]
			# if(region == "tss") {
			# 	pos_cr3 = pos_cr3[pos_cr3$gene_tss_dist > -1000 & pos_cr3$gene_tss_dist < 2000]
			# } else if(region == "gene") {
			# 	pos_cr3 = pos_cr3[pos_cr3$gene_tss_dist > 2000 & pos_cr3$gene_tss_dist < gl[gi[k]]]
			# } else if(region == "intergenic") {
			# 	pos_cr3 = pos_cr3[pos_cr3$gene_tss_dist < -1000 | pos_cr3$gene_tss_dist > gl[gi[k]]]
			# }
			x = gene_tss_dist[[k]]
			if(region == "tss") {
				l = x > -1000 & x < 2000
			} else if(region == "gene") {
				l = x > 2000 & x < gl[gi[k]]
			} else if(region == "intergenic") {
				l = x < -1000 | x > gl[gi[k]]
			}
			n_pos[i, j] = n_pos[i, j] + sum(l)
			w_pos[i, j] = w_pos[i, j] + sum(end[[k]][l] - start[[k]][l])
			qqcat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b")
			qqcat("  @{k}/@{length(gene_tss_dist)} pos_cr gene...")
		}
		cat("\n")
	}
}
r1 = n_neg/n_pos
r2 = w_neg/w_pos
par(mfrow = c(1, 2))
matplot(r1, xlab = "abs_correlation", type = "b", pch = 1, ylab = "#neg/#pos", col = 1:5, axes = FALSE, main = qq("region = @{region}"))
legend("topleft", lty = 1, col = 1:5, legend = colnames(r1))
axis(side = 1, at = 1:7, labels = rownames(r1), cex = 0.8)
axis(side = 2)
box()
matplot(r2, xlab = "abs_correlation", type = "b", pch = 1, ylab = "width(neg)/width(pos)", col = 1:5, axes = FALSE, main = qq("region = @{region}"))
axis(side = 1, at = 1:7, labels = rownames(r1), cex = 0.8)
axis(side = 2)
box()	

return(list(n_neg, n_pos, w_neg, w_pos))
}

pdf(qq("@{OUTPUT_DIR}/plots/cr_qc.pdf"), width = 10, height = 5)
lt_tss = foo("tss")
lt_gene = foo("gene")
lt_intergenic = foo("intergenic")
dev.off()

res = list(lt_tss, lt_gene, lt_intergenic)
saveRDS(res, file = qq("@{OUTPUT_DIR}/rds/cr_compare.rds"))


pdf(qq("@{OUTPUT_DIR}/plots/cr_qc.pdf"), width = 10, height = 4)
cr_qc(neg_cr, pos_cr, TXDB, chromosome = CHROMOSOME, region = "all")
cr_qc(neg_cr, pos_cr, TXDB, chromosome = CHROMOSOME, region = "intergenic")
cr_qc(neg_cr, pos_cr, TXDB, chromosome = CHROMOSOME, region = "tss")
cr_qc(neg_cr, pos_cr, TXDB, chromosome = CHROMOSOME, region = "gene")
dev.off()

