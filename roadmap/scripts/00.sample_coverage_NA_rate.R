library(methods)
suppressPackageStartupMessages(library(GetoptLong))

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
OUTPUT_DIR = BASE_DIR

# this xls file is actually a plain text file
meta = read.table(qq("@{BASE_DIR}/data/EG.mnemonics.name.xls"), sep = "\t", stringsAsFactors = FALSE)
cn = meta[[1]]

n_all = matrix(0, nrow = 22, ncol = length(cn))
rownames(n_all) = paste0("chr", 1:22)
colnames(n_all) = cn
n_na = n_all

for(chr in paste0("chr", 1:22)) {
	qqcat("reading CpG coverage for @{chr}...\n")
	cov = read.table(qq("@{BASE_DIR}/data/ReadCoverage_Removed_E027_E064/@{chr}.rc"))

	cov = as.matrix(cov[, -1])

	n_all[chr, ] = nrow(cov)
	n_na[chr, ] = apply(cov, 2, function(x) sum(x < 0))
}

n_all = rbind(n_all, "all" = colSums(n_all))
n_na = rbind(n_na, "all" = colSums(n_na))
ratio = n_na/n_all

library(ComplexHeatmap)

pdf(qq("@{OUTPUT_DIR}/plots/all_sample_coverage_na_rate.pdf"))
ht = Heatmap(ratio[-23, ], name = "NA_proportion", cluster_rows = FALSE, cluster_columns = FALSE, col = c("white", "red"),
	top_annotation = HeatmapAnnotation(NA_prop_all = anno_barplot(ratio[23,], axis = FALSE)),
	top_annotation_height = unit(3, "cm"))
draw(ht, column_title = "Proportion of CpG sites with no methylation data (NA)")
decorate_annotation("NA_prop_all", {
	grid.lines(unit(c(0, 1), "npc"), unit(c(0.1, 0.1), "native"), gp = gpar(lty = 2, col = "grey"))
	grid.yaxis(main = FALSE, gp = gpar(fontsize = 8))
	grid.text("NA_prop_all", unit(1, "npc") + unit(10, "mm"), unit(0.5, "npc"), just = "bottom", rot = -90)
})
dev.off()
