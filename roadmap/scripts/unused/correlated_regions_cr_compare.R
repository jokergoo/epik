
suppressPackageStartupMessages(library(GetoptLong))
group = 1
GetoptLong("group=i", "1|2|3|4")



library(RColorBrewer)
BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

km = readRDS(qq(paste0(OUTPUT_DIR, "/rds/mat_all_cr_enriched_to_gene_extend_50000_target_0.6_row_order_km3_and_km4.rds")))[[2]]


gm = genes(TXDB)
gm = gm[seqnames(gm) %in% CHROMOSOME]



foo = function(region) {

n_neg = matrix(0, nrow = 6, ncol = 5)
IQR_cutoff = c(0, 0.1, 0.2, 0.3, 0.4)
colnames(n_neg) = paste0("IQR > ", IQR_cutoff)
cor_cutoff = 1:6/10
rownames(n_neg) = paste("|cor| >", cor_cutoff)
n_pos = n_neg
w_neg = n_neg
w_pos = n_neg
c_neg = n_neg
c_pos = n_neg
mean_w_neg = n_neg
q25_w_neg = n_neg
q75_w_neg = n_neg
mean_w_pos = n_pos
q25_w_pos = n_pos
q75_w_pos = n_pos
interdist_neg = n_neg
interdist_pos = n_pos
den_cpg_neg = n_neg
den_cpg_pos = n_neg
for(cor_cutoff in 6:1*0.1) {
	for(iqr_cutoff in c(0, 0.1, 0.2, 0.3, 0.4)) {
		i = paste0("|cor| > ", cor_cutoff)
		j = paste0("IQR > ", iqr_cutoff)

		qqcat("@{i}, @{j}\n")

		neg_cr2 = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/rds/all_neg_cr_w6s3_corr_less_than_-@{cor_cutoff}_methIQR_larger_than_@{iqr_cutoff}_all.rds"))
		neg_cr2 = neg_cr2[neg_cr2$gene_id %in% genes]
		neg_cr2 = neg_cr2[order(neg_cr2$gene_id, start(neg_cr2))]
		# neg_cr2_list = split(neg_cr2, neg_cr2$gene_id)
		gi = unique(neg_cr2$gene_id)
		gl = width(gm[gi])
		names(gl) = gi

		if(region == "tss") {
			l = neg_cr2$gene_tss_dist > -1000 & neg_cr2$gene_tss_dist < 2000
		} else if(region == "gene") {
			l = neg_cr2$gene_tss_dist > 2000 & neg_cr2$gene_tss_dist < gl[neg_cr2$gene_id]
		} else if(region == "intergenic") {
			l = neg_cr2$gene_tss_dist < -1000 | neg_cr2$gene_tss_dist > gl[neg_cr2$gene_id]
		}

		if(region == "tss") {
			w = 3000
		} else if(region == "gene") {
			w = abs(gl - 2000)
		} else if(region == "intergenic") {
			w = 50000*2 - 1000
		}

		if(region %in% c("tss", "gene")) {
			s = start(neg_cr2[l])
			e = end(neg_cr2[l])
			interdist_neg[i, j] = median(unlist(tapply(seq_along(s), neg_cr2[l]$gene_id, function(ind) {
					if(length(ind) == 1) {
						return(NA)
					} else {
						rainfallTransform(data.frame(s[ind], e[ind]))$dist
					}
				})), na.rm = TRUE)
		} else {
			l1 = neg_cr2$gene_tss_dist < -1000
			l2 = neg_cr2$gene_tss_dist > gl[neg_cr2$gene_id]
			s1 = start(neg_cr2[l1])
			e1 = end(neg_cr2[l1])
			s2 = start(neg_cr2[l2])
			e2 = end(neg_cr2[l2])
			interdist_neg[i, j] = median(c(unlist(tapply(seq_along(s1), neg_cr2[l1]$gene_id, function(ind) {
					if(length(ind) == 1) {
						return(NA)
					} else {
						rainfallTransform(data.frame(s1[ind], e1[ind]))$dist
					}
				})), unlist(tapply(seq_along(s2), neg_cr2[l2]$gene_id, function(ind) {
					if(length(ind) == 1) {
						return(NA)
					} else {
						rainfallTransform(data.frame(s2[ind], e2[ind]))$dist
					}
				}))), na.rm = TRUE)
		}

		neg_cr3 = neg_cr2[l]
		n_neg[i, j] = length(neg_cr3)
		w_neg[i, j] = sum(width(neg_cr3))
		den_cpg_neg[i, j] = sum(neg_cr3$n)/sum(width(neg_cr3))
		mean_w_neg[i, j] = median(width(neg_cr3))
		if(region == "gene") {
			c_neg[i, j] = median(tapply(width(neg_cr3), neg_cr3$gene_id, sum)/w[unique(neg_cr3$gene_id)], na.rm = TRUE)
		} else {
			c_neg[i, j] = median(tapply(width(neg_cr3), neg_cr3$gene_id, sum)/w, na.rm = TRUE)
		}

		pos_cr2 = readRDS(qq("/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/rds/all_pos_cr_w6s3_corr_larger_than_@{cor_cutoff}_methIQR_larger_than_@{iqr_cutoff}_all.rds"))
		pos_cr2 = pos_cr2[pos_cr2$gene_id %in% genes]
		pos_cr2 = pos_cr2[order(pos_cr2$gene_id, start(pos_cr2))]
		gi = unique(pos_cr2$gene_id)
		gl = width(gm[gi])
		names(gl) = gi

		if(region == "tss") {
			l = pos_cr2$gene_tss_dist > -1000 & pos_cr2$gene_tss_dist < 2000
		} else if(region == "gene") {
			l = pos_cr2$gene_tss_dist > 2000 & pos_cr2$gene_tss_dist < gl[pos_cr2$gene_id]
		} else if(region == "intergenic") {
			l = pos_cr2$gene_tss_dist < -1000 | pos_cr2$gene_tss_dist > gl[pos_cr2$gene_id]
		}

		if(region == "tss") {
			w = 3000
		} else if(region == "gene") {
			w = abs(gl - 2000)
		} else if(region == "intergenic") {
			w = 50000*2 - 1000
		}

		if(region %in% c("tss", "gene")) {
			s = start(pos_cr2[l])
			e = end(pos_cr2[l])
			interdist_pos[i, j] = median(unlist(tapply(seq_along(s), pos_cr2[l]$gene_id, function(ind) {
					if(length(ind) == 1) {
						return(NA)
					} else {
						rainfallTransform(data.frame(s[ind], e[ind]))$dist
					}
				})), na.rm = TRUE)
		} else {
			l1 = pos_cr2$gene_tss_dist < -1000
			l2 = pos_cr2$gene_tss_dist > gl[pos_cr2$gene_id]
			s1 = start(pos_cr2[l1])
			e1 = end(pos_cr2[l1])
			s2 = start(pos_cr2[l2])
			e2 = end(pos_cr2[l2])
			interdist_pos[i, j] = median(c(unlist(tapply(seq_along(s1), pos_cr2[l1]$gene_id, function(ind) {
					if(length(ind) == 1) {
						return(NA)
					} else {
						rainfallTransform(data.frame(s1[ind], e1[ind]))$dist
					}
				})), unlist(tapply(seq_along(s2), pos_cr2[l2]$gene_id, function(ind) {
					if(length(ind) == 1) {
						return(NA)
					} else {
						rainfallTransform(data.frame(s2[ind], e2[ind]))$dist
					}
				}))), na.rm = TRUE)
		}

		pos_cr3 = pos_cr2[l]
		n_pos[i, j] = length(pos_cr3)
		w_pos[i, j] = sum(width(pos_cr3))
		den_cpg_pos[i, j] = sum(pos_cr3$n)/sum(width(pos_cr3))
		mean_w_pos[i, j] = median(width(pos_cr3))
		if(region == "gene") {
			c_pos[i, j] = median(tapply(width(pos_cr3), pos_cr3$gene_id, sum)/w[unique(pos_cr3$gene_id)], na.rm = TRUE)
		} else {
			c_pos[i, j] = median(tapply(width(pos_cr3), pos_cr3$gene_id, sum)/w, na.rm = TRUE)
		}
	}
}

return(list(n_neg = n_neg, n_pos = n_pos, w_neg = w_neg, w_pos = w_pos, c_neg = c_neg, 
	c_pos = c_pos, mean_w_neg = mean_w_neg, mean_w_pos = mean_w_pos, interdist_neg = interdist_neg,
	interdist_pos = interdist_pos, den_cpg_neg = den_cpg_neg, den_cpg_pos = den_cpg_pos))
}

foo_plot = function(n_neg, n_pos, w_neg, w_pos, c_neg, c_pos, mean_w_neg, mean_w_pos,
	interdist_neg, interdist_pos, den_cpg_neg, den_cpg_pos, region) {

r1 = n_neg/n_pos
r2 = w_neg/w_pos
r3 = c_neg/c_pos
r4 = mean_w_neg/mean_w_pos
r5 = interdist_neg/interdist_pos
r6 = den_cpg_neg/den_cpg_pos
col4 = brewer.pal(12, "Paired")[c(6, 8, 4, 2, 10)]
col8 = brewer.pal(12, "Paired")[c(6, 8, 4, 2, 10, 5, 7, 3, 1, 9)]
matplot(r2, xlab = "abs_correlation", type = "b", pch = 1, ylab = "width(neg)/width(pos)", lty = 1,
	axes = FALSE, main = qq("region = @{region}"), col = col4)
axis(side = 1, at = 1:6, labels = rownames(r1), cex = 0.8)
axis(side = 2)
box()
matplot(r4, xlab = "abs_correlation", type = "b", pch = 1, ylab = "mean_w(neg)/mean_w(pos)", lty = 1,
	axes = FALSE, main = qq("region = @{region}"), col = col4)
axis(side = 1, at = 1:6, labels = rownames(r1), cex = 0.8)
axis(side = 2)
box()	
matplot(cbind(mean_w_neg, mean_w_pos), xlab = "abs_correlation", ylim = range(c(mean_w_neg, mean_w_pos)), log = "y", type = "b", pch = 1, ylab = "mean_w", lty = 1,
	axes = FALSE, main = qq("region = @{region}"), col = col8)
axis(side = 1, at = 1:6, labels = rownames(r1), cex = 0.8)
axis(side = 2)
box()	
matplot(r6, xlab = "abs_correlation", type = "b", pch = 1, ylab = "den_cpg(neg)/den_cpg(pos)", lty = 1,
	axes = FALSE, main = qq("region = @{region}"), col = col4)
axis(side = 1, at = 1:6, labels = rownames(r1), cex = 0.8)
axis(side = 2)
box()	
matplot(cbind(den_cpg_neg, den_cpg_pos), xlab = "abs_correlation", ylim = range(c(den_cpg_neg, den_cpg_pos)), log = "y", type = "b", pch = 1, ylab = "den_cpg", lty = 1,
	axes = FALSE, main = qq("region = @{region}"), col = col8)
axis(side = 1, at = 1:6, labels = rownames(r1), cex = 0.8)
axis(side = 2)
box()
matplot(r3, xlab = "abs_correlation", type = "b", pch = 1, ylab = "cov(neg)/cov(pos)", lty = 1,
	axes = FALSE, main = qq("region = @{region}"), col = col4)
axis(side = 1, at = 1:6, labels = rownames(r1), cex = 0.8)
axis(side = 2)
box()	
matplot(cbind(c_neg, c_pos), xlab = "abs_correlation", ylim = range(c(c_neg, c_pos)), log = "y", type = "b", pch = 1, ylab = "cov", lty = 1,
	axes = FALSE, main = qq("region = @{region}"), col = col8)
axis(side = 1, at = 1:6, labels = rownames(r1), cex = 0.8)
axis(side = 2)
box()	
}

for(group in unique(km)) {
genes = names(km[km == group])

lt_tss = foo("tss")
lt_gene = foo("gene")
lt_intergenic = foo("intergenic")

res = list(tss = lt_tss, gene = lt_gene, intergenic = lt_intergenic)
saveRDS(res, file = qq("@{OUTPUT_DIR}/rds/cr_compare_@{group}.rds"))
#res = readRDS(qq("@{OUTPUT_DIR}/rds/cr_compare_@{group}.rds"))

pdf(qq("@{OUTPUT_DIR}/plots/cr_compare_@{group}.pdf"), width = 20, height = 12)
par(mfrow = c(3, 7))
foo_plot(res$tss$n_neg, res$tss$n_pos, res$tss$w_neg, res$tss$w_pos, res$tss$c_neg, res$tss$c_pos, res$tss$mean_w_neg, res$tss$mean_w_pos,
	res$tss$interdist_neg, res$tss$interdist_pos, res$tss$den_cpg_neg, res$tss$den_cpg_pos, "tss")
foo_plot(res$gene$n_neg, res$gene$n_pos, res$gene$w_neg, res$gene$w_pos, res$gene$c_neg, res$gene$c_pos, res$gene$mean_w_neg, res$gene$mean_w_pos,
	res$gene$interdist_neg, res$gene$interdist_pos, res$gene$den_cpg_neg, res$gene$den_cpg_pos,"gene")
foo_plot(res$intergenic$n_neg, res$intergenic$n_pos, res$intergenic$w_neg, res$intergenic$w_pos, res$intergenic$c_neg, res$intergenic$c_pos, res$intergenic$mean_w_neg, res$intergenic$mean_w_pos,
	res$intergenic$interdist_neg, res$intergenic$interdist_pos, res$intergenic$den_cpg_neg, res$intergenic$den_cpg_pos,"intergenic")
dev.off()

}


# for(group in 1:4) {
#     cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/correlated_regions_cr_compare.R --group @{group}")
    
#     cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=20:00:00,mem=10G -N cr_compare_group_@{group}' '@{cmd}'")
#     system(cmd)
# }


