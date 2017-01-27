
cr_genic_stat = function(cr_reduced, txdb) {

	gm = genes(txdb)
	gl = width(gm)
	names(gl) = names(gm)

	cr2 = cr_reduced
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
	cpg_density_list[["promoter_neg"]] = cr2[l_promoter & l_neg]$ncpg/width(cr2[l_promoter & l_neg])*1000
	cpg_density_list[["promoter_pos"]] = cr2[l_promoter & l_pos]$ncpg/width(cr2[l_promoter & l_pos])*1000
	cpg_density_list[["gene_neg"]] = cr2[l_gene & l_neg]$ncpg/width(cr2[l_gene & l_neg])*1000
	cpg_density_list[["gene_pos"]] = cr2[l_gene & l_pos]$ncpg/width(cr2[l_gene & l_pos])*1000
	cpg_density_list[["intergenic_neg"]] = cr2[l_intergenic & l_neg]$ncpg/width(cr2[l_intergenic & l_neg])*1000
	cpg_density_list[["intergenic_pos"]] = cr2[l_intergenic & l_pos]$ncpg/width(cr2[l_intergenic & l_pos])*1000

	par(mfrow = c(1, 5), mar = c(7, 4, 4, 2))
	boxplot(meth_diff_list, outline = FALSE, col = c("darkgreen", "red"), ylab = "meth mean diff", names = rep("", 6), main = "meth mean diff")
	par(las = 3); axis(side = 1, at = 1:6, labels = names(meth_diff_list)); par(las = 0)
	boxplot(length_list, outline = FALSE, col = c("darkgreen", "red"), ylab = "bp", names = rep("", 6), main = "length")
	par(las = 3); axis(side = 1, at = 1:6, labels = names(length_list)); par(las = 0)
	foo = barplot(unlist(affected_genes), col = c("darkgreen", "red"), ylab = "#gene", names = rep("", 6), main = "affected_genes")
	box()
	par(las = 3); axis(side = 1, at = foo, labels = names(affected_genes)); par(las = 0)
	foo = barplot(sapply(length_list, sum)/1000, col = c("darkgreen", "red"), ylab = "kb", names = rep("", 6), main = "sum(length)")
	box()
	par(las = 3); axis(side = 1, at = foo, labels = names(length_list)); par(las = 0)
	boxplot(cpg_density_list, outline = FALSE, col = c("darkgreen", "red"), ylab = "cpg/kb", names = rep("", 6), main = "#cpg per kb")
	par(las = 3); axis(side = 1, at = 1:6, labels = names(cpg_density_list)); par(las = 0)

}
