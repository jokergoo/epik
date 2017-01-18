library(methods)
library(GetoptLong)

cutoff = 0.05
meandiff = 0.1
rerun = FALSE
GetoptLong("cutoff=f", "0.05",
	       "meandiff=s", "0",
	       "rerun!", "")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))


sig_neg_cr = readRDS(qq("@{OUTPUT_DIR}/rds/all_neg_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds"))
sig_pos_cr = readRDS(qq("@{OUTPUT_DIR}/rds/all_pos_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds"))

sig_cr = c(sig_neg_cr, sig_pos_cr)
gm = genes(TXDB)
gm = gm[seqnames(gm) %in% CHROMOSOME]
g = gm[gm$gene_id %in% unique(sig_cr$gene_id)]
strand(g) = "*"
start(g) = start(g) - 500000
end(g) = end(g) + 50000
start(g) = ifelse(start(g) > 1, start(g), 1)

## correlation between cr and CGI/shore
rdata_file = qq("@{OUTPUT_DIR}/rds/cr_enrich_to_cgi_fdr_@{cutoff}_methdiff_@{meandiff}.rds")
if(file.exists(rdata_file) && !rerun) {
	lt = readRDS(rdata_file)
	r1 = lt$r1
	r2 = lt$r2
} else {
	r1 = epic::genomic_regions_correlation(list(sig_neg_cr = sig_neg_cr, sig_pos_cr = sig_pos_cr), list(cgi = CGI, shore = CGI_SHORE), 
		nperm = 1000, background = g, mc.cores = 2)

	r2 = epic::genomic_regions_correlation(list(sig_neg_cr = sig_neg_cr, sig_pos_cr = sig_pos_cr), list(cgi = CGI, shore = CGI_SHORE), 
		nperm = 1000, mc.cores = 2)

	lt = list(r1 = r1, r2 = r2)
	saveRDS(lt, file = qq("@{OUTPUT_DIR}/rds/cr_enrich_to_cgi_fdr_@{cutoff}_methdiff_@{meandiff}.rds"))
}

pdf(qq("@{OUTPUT_DIR}/plots/cr_enrich_to_cgi_fdr_@{cutoff}_methdiff_@{meandiff}.pdf"))
barplot(r1$fold_change, beside=TRUE, ylim = c(0, max(r1$fold_change)*1.1), ylab = "foldchage", main = "extended gene model as background", col = c("pink", "orange"))
legend("topright", pch = 15, legend = c("CGI", "shore"), col = c("pink", "orange"))
box()

barplot(r2$fold_change, beside=TRUE, ylim = c(0, max(r2$fold_change)*1.1), ylab = "foldchage", main = "extended gene model", col = c("pink", "orange"))
legend("topright", pch = 15, legend = c("CGI", "shore"), col = c("pink", "orange"))
box()
dev.off()

## how neg_cr enriched at tss/gene body
pdf(qq("@{OUTPUT_DIR}/plots/cr_tss_fdr_@{cutoff}_methdiff_@{meandiff}.pdf"), width = 10, height = 4)
cr_enriched_at_tss(sig_neg_cr, sig_pos_cr, TXDB)
dev.off()


pdf(qq("@{OUTPUT_DIR}/plots/cr_genebody_fdr_@{cutoff}_methdiff_@{meandiff}.pdf"), width = 10, height = 4)
cr_enriched_at_gene_body(sig_neg_cr, sig_pos_cr, TXDB)
dev.off()

## how many transcripts that have CR for each gene
tx = transcripts(TXDB)
tx_tss = promoters(tx, upstream = 1, downstream = 0)
tx_tss = start(tx_tss)
names(tx_tss) = tx$tx_name
names(tx) = tx$tx_name
strand = as.vector(strand(tx))
names(strand) = tx$tx_name

pdf(qq("@{OUTPUT_DIR}/plots/cr_gene_number_of_tx_@{cutoff}_methdiff_@{meandiff}.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

x = mcols(sig_neg_cr)[sig_neg_cr$tx_tss_dist > -5000 & sig_neg_cr$tx_tss_dist < 5000, c("gene_id", "nearest_tx_tss")]
x = unique(x)
x$nearest_tx_tss_pos = tx_tss[x$nearest_tx_tss]
x$strand = as.vector(strand(gm[x$gene_id]))
x$gene_tss = start(promoters(gm[x$gene_id], upstream = 1, downstream = 0))
x$diff = abs(x$gene_tss - x$nearest_tx_tss_pos)
x2 = unique(x[c(1, 3)])
foo = structure(x2$gene_id, names = x2$nearest_tx_tss_pos)
count_neg = table(table(foo))
barplot(count_neg, xlab = "#tx", ylab = "count", col = "darkgreen", main = "number of genes with k tx\nthat tx tss (-5000, 5000) overlap with neg_cr")

x = mcols(sig_pos_cr)[sig_pos_cr$tx_tss_dist > -5000 & sig_pos_cr$tx_tss_dist < 5000, c("gene_id", "nearest_tx_tss")]
x = unique(x)
x$nearest_tx_tss_pos = tx_tss[x$nearest_tx_tss]
x$strand = as.vector(strand(gm[x$gene_id]))
x$gene_tss = start(promoters(gm[x$gene_id], upstream = 1, downstream = 0))
x$diff = abs(x$gene_tss - x$nearest_tx_tss_pos)
x2 = unique(x[c(1, 3)])
foo = structure(x2$gene_id, names = x2$nearest_tx_tss_pos)
count_pos = table(table(foo))
barplot(count_pos, xlab = "#tx", ylab = "count", col = "red", main = "number of genes with k tx\nthat tx tss (-5000, 5000) overlap with pos_cr")

dev.off()




gene = genes(TXDB)
tss = promoters(gene, upstream = 1, downstream = 0)
tss = tss[names(tss) %in% sig_neg_cr$gene_id]
mat_neg = normalizeToMatrix(sig_neg_cr, tss, mapping_column = "gene_id",
		 extend = 5000, mean_mode = "absolute", w = 50, trim = 0)
mat_pos = normalizeToMatrix(sig_pos_cr, tss, mapping_column = "gene_id",
		 extend = 5000, mean_mode = "absolute", w = 50, trim = 0)

pdf(qq("@{OUTPUT_DIR}/plots/cr_only_heatmap_fdr_@{cutoff}_methdiff_@{meandiff}.pdf"), width = 8, height = 16)

l = rowSums(mat_neg) > 0
ht_list = EnrichedHeatmap(mat_neg[l, ], name = "neg", column_title = "neg_cr", col = c("0" = "white", "1" = "darkgreen"),
	top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "darkgreen"))), 
	top_annotation_height = unit(2, "cm"), use_raster = TRUE, raster_quality = 2) + 
EnrichedHeatmap(mat_pos[l, ], name = "pos", column_title = "pos_cr", col = c("0" = "white", "1" = "red"),
	top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "red"))), 
	top_annotation_height = unit(2, "cm"), use_raster = TRUE, raster_quality = 2)	

draw(ht_list, main_heatmap = "neg", cluster_rows = TRUE, show_row_dend = FALSE)

l = rowSums(mat_pos) > 0
ht_list = EnrichedHeatmap(mat_neg[l, ], name = "neg", column_title = "neg_cr", col = c("0" = "white", "1" = "darkgreen"),
	top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "darkgreen"))), 
	top_annotation_height = unit(2, "cm"), use_raster = TRUE, raster_quality = 2) + 
EnrichedHeatmap(mat_pos[l, ], name = "pos", column_title = "pos_cr", col = c("0" = "white", "1" = "red"),
	top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "red"))), 
	top_annotation_height = unit(2, "cm"), use_raster = TRUE, raster_quality = 2)	

draw(ht_list, main_heatmap = "pos", cluster_rows = TRUE, show_row_dend = FALSE)
dev.off()

# how many of them overlap to CGI

# for(cutoff in c(0.1, 0.05, 0.01)) {
#     for(meandiff in c(0, 0.1, 0.2, 0.3)) {      
# 		cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/03.correlated_regions_misc.R --cutoff @{cutoff} --meandiff @{meandiff}")
# 		cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=30:00:00,mem=10G -N correlated_regions_misc_fdr_@{cutoff}_methdiff_@{meandiff}' '@{cmd}'")
# 		system(cmd)
# 	}
# }


