library(methods)
library(GetoptLong)

cutoff = 0.05
meandiff = 0.2
GetoptLong("cutoff=f", "0.05",
	       "meandiff=s", "0")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

neg_cr = readRDS(file = qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3.rds")))
pos_cr = readRDS(file = qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3.rds")))

# make a general heatmap for the correlation landscape

cr = c(neg_cr, pos_cr)

gm = genes(TXDB)
gm = gm[seqnames(gm) %in% CHROMOSOME]
g = gm[gm$gene_id %in% unique(cr$gene_id)]
gl = width(g)
names(gl) = names(g)

cr2 = cr[cr$corr_fdr < cutoff & cr$meth_anova_fdr < cutoff & cr$meth_diameter > meandiff]
l_neg = cr2$corr < 0
l_pos = cr2$corr > 0

l_promoter = cr2$gene_tss_dist > -1000 & cr2$gene_tss_dist < 2000
l_gene = cr2$gene_tss_dist > 2000 & cr2$gene_tss_dist < gl[cr2$gene_id]
l_intergenic = cr2$gene_tss_dist < -1000 | cr2$gene_tss_dist > gl[cr2$gene_id]

meth_diff_list = list()
meth_diff_list[["promoter_neg"]] = cr2$meth_diameter[l_promoter & l_neg]
meth_diff_list[["promoter_pos"]] = cr2$meth_diameter[l_promoter & l_pos]
meth_diff_list[["gene_neg"]] = cr2$meth_diameter[l_gene & l_neg]
meth_diff_list[["gene_pos"]] = cr2$meth_diameter[l_gene & l_pos]
meth_diff_list[["intergenic_neg"]] = cr2$meth_diameter[l_intergenic & l_neg]
meth_diff_list[["intergenic_pos"]] = cr2$meth_diameter[l_intergenic & l_pos]

neg_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds")))
pos_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds")))

cr2 = c(neg_cr, pos_cr)
l_neg = c(rep(TRUE, length(neg_cr)), rep(FALSE, length(pos_cr)))
l_pos = c(rep(FALSE, length(neg_cr)), rep(TRUE, length(pos_cr)))

l_promoter = cr2$gene_tss_dist > -1000 & cr2$gene_tss_dist < 2000
l_gene = cr2$gene_tss_dist > 2000 & cr2$gene_tss_dist < gl[cr2$gene_id]
l_intergenic = cr2$gene_tss_dist < -1000 | cr2$gene_tss_dist > gl[cr2$gene_id]

length_list = list()
length_list[["promoter_neg"]] = width(cr2[l_promoter & l_neg])
length_list[["promoter_pos"]] = width(cr2[l_promoter & l_pos])
length_list[["gene_neg"]] = width(cr2[l_gene & l_neg])
length_list[["gene_pos"]] = width(cr2[l_gene & l_pos])
length_list[["intergenic_neg"]] = width(cr2[l_intergenic & l_neg])
length_list[["intergenic_pos"]] = width(cr2[l_intergenic & l_pos])

affected_genes = list()
affected_genes[["promoter_neg"]] = length(unique(cr2[l_promoter & l_neg]$gene_id))
affected_genes[["promoter_pos"]] = length(unique(cr2[l_promoter & l_pos]$gene_id))
affected_genes[["gene_neg"]] = length(unique(cr2[l_gene & l_neg]$gene_id))
affected_genes[["gene_pos"]] = length(unique(cr2[l_gene & l_pos]$gene_id))
affected_genes[["intergenic_neg"]] = length(unique(cr2[l_intergenic & l_neg]$gene_id))
affected_genes[["intergenic_pos"]] = length(unique(cr2[l_intergenic & l_pos]$gene_id))


cpg_density_list = list()
cpg_density_list[["promoter_neg"]] = cr2[l_promoter & l_neg]$n/width(cr2[l_promoter & l_neg])*1000
cpg_density_list[["promoter_pos"]] = cr2[l_promoter & l_pos]$n/width(cr2[l_promoter & l_pos])*1000
cpg_density_list[["gene_neg"]] = cr2[l_gene & l_neg]$n/width(cr2[l_gene & l_neg])*1000
cpg_density_list[["gene_pos"]] = cr2[l_gene & l_pos]$n/width(cr2[l_gene & l_pos])*1000
cpg_density_list[["intergenic_neg"]] = cr2[l_intergenic & l_neg]$n/width(cr2[l_intergenic & l_neg])*1000
cpg_density_list[["intergenic_pos"]] = cr2[l_intergenic & l_pos]$n/width(cr2[l_intergenic & l_pos])*1000


pdf(qq("@{OUTPUT_DIR}/plots/sig_cr_general_stat_fdr_@{cutoff}_methdiff_@{meandiff}.pdf"), width = 16, height = 4)
par(mfrow = c(1, 5), mar = c(7, 4, 4, 2))
boxplot(meth_diff_list, outline = FALSE, col = c("darkgreen", "red"), ylab = "meth mean diff", names = rep("", 6), main = "meth mean diff")
par(las = 3); axis(side = 1, at = 1:6, labels = names(meth_diff_list)); par(las = 0)
boxplot(length_list, outline = FALSE, col = c("darkgreen", "red"), ylab = "bp", names = rep("", 6), main = "length")
par(las = 3); axis(side = 1, at = 1:6, labels = names(length_list)); par(las = 0)
foo = barplot(unlist(affected_genes), col = c("darkgreen", "red"), ylab = "#gene", names = rep("", 6), main = "affected_genes")
par(las = 3); axis(side = 1, at = foo, labels = names(affected_genes)); par(las = 0)
foo = barplot(sapply(length_list, sum)/1000, col = c("darkgreen", "red"), ylab = "kb", names = rep("", 6), main = "sum(length)")
par(las = 3); axis(side = 1, at = foo, labels = names(length_list)); par(las = 0)
boxplot(cpg_density_list, outline = FALSE, col = c("darkgreen", "red"), ylab = "cpg/kb", names = rep("", 6), main = "#cpg per kb")
par(las = 3); axis(side = 1, at = 1:6, labels = names(cpg_density_list)); par(las = 0)
dev.off()


# for(cutoff in c(0.1, 0.05, 0.01)) {
#     for(meandiff in c(0, 0.1, 0.2, 0.3)) {
#         cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/04.correlated_regions_genic_stat.R --cutoff @{cutoff} --meandiff @{meandiff}")
#         cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=2:00:00,mem=5G -N correlated_regions_genic_stat_fdr_@{cutoff}_meandiff_@{meandiff}' '@{cmd}'")
#         system(cmd)
#     }
# }

