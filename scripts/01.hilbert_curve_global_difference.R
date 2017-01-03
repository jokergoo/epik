library(methods)
library(animation)
library(Xmisc)
library(GetoptLong)

rerun = FALSE
GetoptLong("rerun!", "rerun")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

sample_id = rownames(SAMPLE)

sample_id_subgroup1 = intersect(sample_id, rownames(SAMPLE[SAMPLE$subgroup == "subgroup1", ]))
sample_id_subgroup2 = intersect(sample_id, rownames(SAMPLE[SAMPLE$subgroup == "subgroup2", ]))

generate_diff_color_fun = function(x) {
	q = quantile(x, c(0.05, 0.95), na.rm = TRUE)
	max_q = max(abs(q))
	colorRamp2(c(-max_q, 0, max_q), c("#3794bf", "#FFFFFF", "#df8640"))
}

chromInfo = getChromInfoFromUCSC(GENOME)
chromInfo = chromInfo[chromInfo$chrom %in% CHROMOSOME, ]
chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
chromGr_1kb_window = makeWindows(chromGr, w = 1000, short.keep = TRUE)
mcols(chromGr_1kb_window) = NULL

n_states = 5

rdata_file = qq("@{OUTPUT_DIR}/rds/global_difference.RData")
if(file.exists(rdata_file) && !rerun) {
	load(rdata_file)
} else {

	## difference of methylation
	# mean methylaiton by 1kb window, compared to histone modifications which are also segmented by 1kb window
	gr_meth = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, 
		gf_list = list(chromGr_1kb_window = chromGr_1kb_window), filter_fun = function(s) TRUE)[[1]]
	mat = mcols(gr_meth)
	mat = as.matrix(mat[, -ncol(mat)])
	mcols(gr_meth) = NULL
	gr_meth$mean1 = rowMeans(mat[, sample_id_subgroup1])
	gr_meth$mean2 = rowMeans(mat[, sample_id_subgroup2])
	gr_meth$diff = gr_meth$mean1 - gr_meth$mean2

	# there are some regions in chromGr_1kb_window has been removed
	# because there is no CpG in it. However, we need to recover it
	# which is important for HMM segmentation
	mtch = as.matrix(findOverlaps(gr_meth, chromGr_1kb_window))
	gr_meth2 = chromGr_1kb_window
	gr_meth2$mean1 = NA
	gr_meth2$mean2 = NA
	gr_meth2$diff = 0
	gr_meth2$mean1[mtch[, 2]] = gr_meth$mean1[mtch[,1]]
	gr_meth2$mean2[mtch[, 2]] = gr_meth$mean2[mtch[,1]]
	gr_meth2$diff[mtch[, 2]] = gr_meth$diff[mtch[,1]]

	gr_meth = gr_meth2

	qqcat("HMM segmentation with @{n_states} states for meth\n")
	x = gr_meth$diff
	x[x == 0] = rnorm(sum(x == 0), 0, 0.01)
	gr_meth = segment_by_hmm(gr_meth, x, nStates = n_states)

	# histone modifications
	# diff_list = vector("list", length(MARKS))
	gr_list = vector("list", length(MARKS))
	for(i in seq_along(MARKS)) {
		qqcat("making hilbert curve for @{MARKS[i]}\n")

		hm_sample_id = intersect(chipseq_hooks$sample_id(MARKS[i]), sample_id_subgroup1)
		hm_list = lapply(hm_sample_id, function(sid) chipseq_hooks$peak(MARKS[i], sid))
		names(hm_list) = hm_sample_id
		n1 = length(hm_list)

		den_mean = 0

		x1 = 0
		for(j in seq_along(hm_list)) {
			qqcat("  @{hm_sample_id[j]}, group 1\n")
	        x1 = x1 + average_in_window(chromGr_1kb_window, hm_list[[j]], hm_list[[j]]$density, chromosome = CHROMOSOME)
	        den_mean = den_mean + sum(as.numeric(width(hm_list[[j]]))*hm_list[[j]]$density)/sum(as.numeric(width(hm_list[[j]])))
		}
		
		hm_sample_id = intersect(chipseq_hooks$sample_id(MARKS[i]), sample_id_subgroup2)
		hm_list = lapply(hm_sample_id, function(sid) chipseq_hooks$peak(MARKS[i], sid))
		names(hm_list) = hm_sample_id
		n2 = length(hm_list)
		x2 = 0
		for(j in seq_along(hm_list)) {
			qqcat("  @{hm_sample_id[j]}, group 2\n")
			x2 = x2 + average_in_window(chromGr_1kb_window, hm_list[[j]], hm_list[[j]]$density, chromosome = CHROMOSOME)
			den_mean = den_mean + sum(as.numeric(width(hm_list[[j]]))*hm_list[[j]]$density)/sum(as.numeric(width(hm_list[[j]])))
		}

		gr = chromGr_1kb_window
		gr$abs_diff = x1/n1 - x2/n2
		den_mean = den_mean/(n1 + n2)
		gr$diff = gr$abs_diff/den_mean

		
		gr_list[[i]] = gr
	}
	# names(diff_list) = MARKS
	names(gr_list) = MARKS
	
	for(i in seq_along(gr_list)) {
		qqcat("HMM segmentation with @{n_states} states for @{names(gr_list)[i]}\n")
		
		x = gr_list[[i]]$abs_diff
		l1 = x > 0
		l2 = x < 0
		l3 = x == 0
		x[l1] = log(x[l1] + 1)
		x[l2] = -log(-x[l2] + 1)
		x[l3] = rnorm(sum(l3), 0, 0.01)
		
		gr = segment_by_hmm(gr_list[[i]], x, nStates = n_states)
		gr_list[[i]] = gr
	}
	gr_list$meth = gr_meth
	save(gr_list, file = rdata_file)
}


pdf(qq("@{OUTPUT_DIR}/plots/hilbert_all_diff.pdf"), width = 7.5, height = 7)

qqcat("making hilbert curve for methylation\n")

gr_meth = gr_list$meth
# l1 and l2 basically are identical
l1 = !is.na(gr_meth$mean1)
l2 = !is.na(gr_meth$mean2)

col_fun = generate_diff_color_fun(gr_meth$diff[l1 & l2])
cm = ColorMapping(col_fun = col_fun)
lgd = color_mapping_legend(cm, title = "diff", plot = FALSE)

hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "pixel", level = 10, title = qq("mean methylation difference, mean = @{sprintf('%.2f', mean(c(gr_meth$mean1, gr_meth$mean2), na.rm = TRUE))}"), legend = lgd)
hc_layer(hc, gr_meth[l1 & l2], col = col_fun(gr_meth$diff[l1 & l2]))
hc_map(hc, add = TRUE, fill = NA, border = "#808080")
seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "normal", level = 6, newpage = FALSE)
hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
cm = ColorMapping(col_fun = col_fun)
lgd = color_mapping_legend(cm, title = "mean1", plot = FALSE)
# mean methylation in group1
hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "pixel", level = 10, title = qq("mean methylation in group1"), legend = lgd)
hc_layer(hc, gr_meth[l1], col = col_fun(gr_meth$mean1[l1]))
hc_map(hc, add = TRUE, fill = NA, border = "#808080")
seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "normal", level = 6, newpage = FALSE)
hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
cm = ColorMapping(col_fun = col_fun)
lgd = color_mapping_legend(cm, title = "mean2", plot = FALSE)
# mean methylation in group 2
hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "pixel", level = 10, title = qq("mean methylation in group2"), legend = lgd)
hc_layer(hc, gr_meth[l2], col = col_fun(gr_meth$mean2[l2]))
hc_map(hc, add = TRUE, fill = NA, border = "#808080")
seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "normal", level = 6, newpage = FALSE)
hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))

# distribution of the methylation difference
# plot(density(meth_gr$meth_diff), xlab = "mean difference", main = "distribution of mean difference, methylation")
densityHeatmap(split(gr_meth$diff[l1 & l2], as.vector(seqnames(gr_meth[l1 & l2]))), cluster_columns = TRUE, 
	show_column_dend = TRUE, range = quantile(gr_meth$diff[l1 & l2], c(0.1, 0.9)),
	title = "density distribution of mean difference, methylation",
	ylab = "methylation difference")

gr = gr_meth
col_fun = generate_diff_color_fun(gr$diff)
for(j in 1:n_states) {
	qqcat("making hilbert curve for meth, state @{j}\n")
	l = gr$states == j
	cm = ColorMapping(col_fun = col_fun)
	lgd = color_mapping_legend(cm, title = "diff", plot = FALSE)

	hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "pixel", level = 10, title = qq("meth, state @{j}, mean = @{sprintf('%.1e', mean(gr$diff[l]))}"), legend = lgd)
	hc_layer(hc, gr[l], col = col_fun(gr$diff[l]))
	hc_map(hc, add = TRUE, fill = NA, border = "#808080")
	seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
	hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "normal", level = 6, newpage = FALSE)
	hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))
}
# den_x = matrix(nrow = 512, ncol = length(MARKS))
# den_y = den_x
# for(i in seq_along(diff_list)) {
# 	den  = density(diff_list[[i]])
# 	den_x[, i] = den$x
# 	den_y[, i] = den$y
# }
# matplot(den_x, den_y, type = "l", xlab = "mean density difference", ylab = "density", xlim = c(-0.5, 0.5),
# 	main = "density distribution of histone marks and DNase", col = seq_along(MARKS), lty = 1)
# legend("topleft", lty = 1, col = seq_along(MARKS), legend = MARKS)
# abline(v = 0, lty = 2, col = "#CCCCCC")

library(gtrellis)
for(i in seq_along(MARKS)) {

	gr = gr_list[[MARKS[i]]]
	col_fun = generate_diff_color_fun(gr$diff)
	cm = ColorMapping(col_fun = col_fun)
	lgd = color_mapping_legend(cm, title = qq("relative diff"), plot = FALSE)
	foo = which.max(gr$abs_diff)
	den_mean = gr$abs_diff[foo]/gr$diff[foo]
	hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "pixel", level = 10, title = qq("@{MARKS[i]} mean density difference, mean = @{sprintf('%.1e', mean(den_mean))}"), legend = lgd)
	hc_layer(hc, gr, col = col_fun(gr$diff), mean_mode = "w0")
	hc_map(hc, add = TRUE, fill = NA, border = "#808080")
	seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
	hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "normal", level = 6, newpage = FALSE)
	hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))

	# if den_diff ~ 0, it basically means there is no peaks in both groups
	l = abs(gr$diff) > 1e-6
	# diff_list[[i]] = gr$den_diff[l]

	densityHeatmap(split(gr$diff[l], as.vector(seqnames(gr))[l]), cluster_columns = TRUE, 
		show_column_dend = TRUE, range = quantile(gr$diff[l], c(0.1, 0.9)),
		title = qq("density distribution of mean density difference, @{MARKS[i]}"),
		ylab = "mean density difference")

	col_fun = generate_diff_color_fun(gr$diff)
	ylim = quantile(abs(gr$diff[l]), 0.99)
	for(j in 1:n_states) {
		qqcat("making hilbert curve for @{MARKS[i]}, state @{j}\n")
		l = gr$states == j
		cm = ColorMapping(col_fun = col_fun)
		lgd = color_mapping_legend(cm, title = "diff", plot = FALSE)

		hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "pixel", level = 10, title = qq("@{MARKS[i]}, state @{j}, mean = @{sprintf('%.1e', mean(gr$diff[l]))}"), legend = lgd)
		hc_layer(hc, gr[l], col = col_fun(gr$diff[l]))
		hc_map(hc, add = TRUE, fill = NA, border = "#808080")
		seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
		hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "normal", level = 6, newpage = FALSE)
		hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))

		print(mean(gr$diff[l]))
		gtrellis_layout(category = CHROMOSOME, nrow = 5, compact = TRUE, n_track = 1,
			track_ylim = c(-ylim, ylim), track_ylab = "rel_diff", title = qq("@{MARKS[i]}, state @{j}, mean = @{sprintf('%.1e', mean(gr$diff[l]))}"),
			add_ideogram_track = TRUE, add_name_track = TRUE)
		add_track(gr[l], panel_fun = function(gr2) {
			print(length(gr2))
			x = (start(gr2) + end(gr2))/2
			grid.points(x, gr2$diff, gp = gpar(fontsize = unit(0.5, "mm"), col = "#00000020"))
		})
	}
}

## expression
gene = genes(TXDB)
expr_diff = rowMeans(EXPR[, intersect(colnames(EXPR), sample_id_subgroup1)]) - 
            rowMeans(EXPR[, intersect(colnames(EXPR), sample_id_subgroup2)])
gene = gene[rownames(EXPR)]
expr_diff = expr_diff[rownames(EXPR)]

col_fun = generate_diff_color_fun(expr_diff)
cm = ColorMapping(col_fun = col_fun)
lgd = color_mapping_legend(cm, title = qq("expr_diff"), plot = FALSE)
hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "pixel", level = 10, title = qq("mean expression diff, mean = @{sprintf('%.1f', mean(EXPR[names(gene), ]))}"), legend = lgd)
hc_layer(hc, gene, col = col_fun(expr_diff), mean_mode = "w0")
hc_map(hc, add = TRUE, fill = NA, border = "#808080")
seekViewport(qq("hilbert_curve_@{HilbertCurve:::.ENV$I_PLOT}"))
hc = GenomicHilbertCurve(chr = CHROMOSOME, mode = "normal", level = 6, newpage = FALSE)
hc_map(hc, add = TRUE, fill = NA, border = NA, labels_gp = gpar(fontsize = 20))

dev.off()





q(save = "no")

n_gr = length(gr_list)
gr_name = names(gr_list)
gr_name_combine = outer(gr_name, gr_name, paste, sep = "_vs_")
# corr_mat = matrix(nrow = length(gr_name_combine), ncol = 4, dimnames = list(gr_name_combine, c("5/5", "5/1", "1/5", "1/1")))
corr_mat = matrix(nrow = length(gr_name_combine), ncol = 10, dimnames = list(gr_name_combine, 0:9*0.1))
jaccard_intersect = corr_mat
jaccard_union = corr_mat
jaccard_gr1 = corr_mat
jaccard_gr2 = corr_mat
for(k in 1:10) {
	for(i in 1:(n_gr-1)) {
		for(j in (i+1):n_gr) {
			gr1 = gr_list[[i]]
			gr2 = gr_list[[j]]
			# gr1 = gr1[gr1$states == ifelse(k %in% 1:2, 5, 1)]
			# gr2 = gr2[gr2$states == ifelse(k %in% c(1, 3), 5, 1)]
			nm = qq("@{gr_name[i]}_vs_@{gr_name[j]}")
			x = abs(gr1$diff)
			q = quantile(x[x > 1e-6], (k-1)*0.1)
			gr1 = gr1[gr1$diff < -q]
			x = abs(gr2$diff)
			q = quantile(x[x > 1e-6], (k-1)*0.1)
			gr2 = gr2[gr2$diff > q]

			r = sum(as.numeric(width(gr1)))/sum(as.numeric(width(gr2)))
			if(r > 1) r = 1/r
			corr_mat[nm, k] = epic::genomic_corr_jaccard(gr1, gr2)
			corr_mat[nm, k] = corr_mat[nm, k]/r
			
			jaccard_intersect[nm, k] = sum(as.numeric(width(intersect(gr1, gr2))))
			jaccard_union[nm, k] = sum(as.numeric(width(union(gr1, gr2))))
			jaccard_gr1[nm, k] = sum(as.numeric(width(gr1)))
			jaccard_gr2[nm, k] = sum(as.numeric(width(gr2)))

			qqcat("@{nm}, @{k}\n")
		}
	}
}

l = is.na(corr_mat[, 1])
corr_mat = corr_mat[!l, ]
jaccard_intersect = jaccard_intersect[!l, ]
jaccard_union = jaccard_union[!l, ]
jaccard_gr1 = jaccard_gr1[!l, ]
jaccard_gr2 = jaccard_gr2[!l, ]


# mat_55 = matrix(nrow = length(MARKS), ncol = length(MARKS), dimnames = list(MARKS, MARKS))
# mat_51 = mat_55
# mat_15=  mat_55
# mat_11 = mat_55
# for(i in seq_along(MARKS)) {
# 	for(j in seq_along(MARKS)) {
# 		if(i == j) {
# 			mat_55[i, j] = 0
# 			mat_51[i, j] = 0
# 			mat_15[i, j] = 0
# 			mat_11[i, j] = 0
# 			next
# 		}
# 		ind = qq("@{MARKS[i]}_vs_@{MARKS[j]}")
# 		if(ind %in% rownames(corr_mat)) {
# 			mat_55[i, j] = corr_mat[qq("@{MARKS[i]}_vs_@{MARKS[j]}"), 1]
# 			mat_51[i, j] = corr_mat[qq("@{MARKS[i]}_vs_@{MARKS[j]}"), 2]
# 			mat_15[i, j] = corr_mat[qq("@{MARKS[i]}_vs_@{MARKS[j]}"), 3]
# 			mat_11[i, j] = corr_mat[qq("@{MARKS[i]}_vs_@{MARKS[j]}"), 4]
# 		} else {
# 			mat_55[i, j] = corr_mat[qq("@{MARKS[j]}_vs_@{MARKS[i]}"), 1]
# 			mat_51[i, j] = corr_mat[qq("@{MARKS[j]}_vs_@{MARKS[i]}"), 2]
# 			mat_15[i, j] = corr_mat[qq("@{MARKS[j]}_vs_@{MARKS[i]}"), 3]
# 			mat_11[i, j] = corr_mat[qq("@{MARKS[j]}_vs_@{MARKS[i]}"), 4]
# 		}
# 	}
# }

# add_venn = function(x1, x2, union, max, x, y, w, h) {
# 	grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "white"))
# 	grid.rect(x - 0.5*w, y, width = x1/union*w, height = h, gp = gpar(fill = "#FF000040"), just = c("left", "center"))
# 	grid.rect(x + 0.5*w, y, width = x2/union*w, height = h, gp = gpar(fill = "#0000FF40"), just = c("right", "center"))
# 	r = union/max
# 	grid.lines(unit.c(x - 0.5*w, x - 0.5*w + r*w), unit.c(y - 0.5*h, y - 0.5*h), gp = gpar(col = "red", lwd = 4))
# }

# cell_fun = function(i, j, k, x, y, w, h) {
# 	if(i == j) return(NULL)
# 	qqcat("@{MARKS[i]}_vs_@{MARKS[j]}\n")
# 	ind = qq("@{MARKS[i]}_vs_@{MARKS[j]}")
# 	if(ind %in% rownames(corr_mat)) {
# 		x1 = jaccard_gr1[qq("@{MARKS[i]}_vs_@{MARKS[j]}"), k]
# 		x2 = jaccard_gr2[qq("@{MARKS[i]}_vs_@{MARKS[j]}"), k]
# 		union = jaccard_union[qq("@{MARKS[i]}_vs_@{MARKS[j]}"), k]
# 	} else {
# 		x1 = jaccard_gr1[qq("@{MARKS[j]}_vs_@{MARKS[i]}"), k]
# 		x2 = jaccard_gr2[qq("@{MARKS[j]}_vs_@{MARKS[i]}"), k]
# 		union = jaccard_union[qq("@{MARKS[j]}_vs_@{MARKS[i]}"), k]
# 	}
# 	add_venn(x1, x2, union, max(jaccard_union), x, y, w, h)
# }

# col_fun = colorRamp2(c(0, 1), c("white", "red"))
# Heatmap(rbind(mat_55, mat_51), col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, 
# 	# split = factor(paste("state =", rep(c(5, 1), each = length(MARKS))), levels = c("state = 5", "state = 1")), 
# 	show_row_names = FALSE, column_title = "state = 5",
# 	cell_fun = function(i, j, x, y, w, h, fill) {
# 		if(i  %% 7 + 1> j) cell_fun(i %% 7 + 1, j, ifelse(i > 7, 2, 1), x, y, w, h)
# 	}) +
# Heatmap(rbind(mat_15, mat_11), col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
# 	column_title = "state = 1",
# 	cell_fun = function(i, j, x, y, w, h, fill) {
# 		if(i %% 7 + 1 > j) cell_fun(i %% 7 + 1, j, ifelse(i > 7, 4, 3), x, y, w, h)
# 	})

# ht_list = NULL
# for(j in 1:4) {
# 	ht_list = ht_list + Heatmap(corr_mat[, j, drop = FALSE], col = col_fun, 
# 		show_row_names = j == 1, row_names_side = "left", row_names_gp = gpar(fontsize = 8),
# 		cluster_rows = FALSE) +
# 	rowAnnotation(venn = local({j = j; function(index) {
# 		print(j)
# 		n = length(index)
# 		pushViewport(viewport())
# 		for(i in index) {
# 			nm = rownames(corr_mat)[i]
# 			x1 = jaccard_gr1[nm, j]
# 			x2 = jaccard_gr2[nm, j]
# 			union = jaccard_union[nm, j]
# 			add_venn(x1, x2, union, max(jaccard_union), unit(0.5, "npc"), unit((n - i + 1 - 0.5)/n, "npc"), unit(1, "npc"), unit(1/n, "npc"))
# 		}
# 		upViewport()
# 	}}), width = unit(2, "cm"))
# }

max_u = max(jaccard_union)
foo_plot = function(k = 1, l1) {
	par(mar = c(5, 4, 4, 10))
	plot(NULL, xlim = c(0, 0.9), ylim = c(0, max(corr_mat, na.rm = TRUE)),
		xlab = "quantile of absolute density difference", ylab = "jaccard coefficient",
		main = "correlation between histone marks")
	# l1 = grepl("meth", rownames(corr_mat))
	for(i in which(l1)) {
		lines(0:9*0.1, corr_mat[i, ])
	}
	abline(v = k*0.1 - 0.1, col = "grey", lty = 1)
	par(xpd = NA)
	text(k*0.1-0.1, max(corr_mat, na.rm = TRUE), qq("q = @{k*0.1-0.1}"), adj = c(-0.1, 1), cex = 0.6, col = "grey")
	points(rep(k*0.1 - 0.1, sum(l1)), corr_mat[l1, k], pch = 16, cex = 0.8)
	if(all(c("H3K27me3_vs_H3K36me3", "H3K4me1_vs_H3K27ac") %in% rownames(corr_mat)[l1])) {
		lines(0:9*0.1, corr_mat["H3K27me3_vs_H3K36me3", ], col = "green", lwd = 2)
		lines(0:9*0.1, corr_mat["H3K4me1_vs_H3K27ac", ], col = "red", lwd = 2)
		points(k*0.1 - 0.1, corr_mat["H3K27me3_vs_H3K36me3", k], pch = 16, col = "green", cex = 0.8)
		points(k*0.1 - 0.1, corr_mat["H3K4me1_vs_H3K27ac", k], pch = 16, col = "red", cex = 0.8)
	}
	od = order(corr_mat[l1, 10])
	x = corr_mat[l1, 10][od]
	labels = rownames(corr_mat)[l1][od]
	text_height = strheight("foo", cex = 0.7) * (1 + 2)
	h1 = x - text_height * 0.5
	h2 = x + text_height * 0.5
	pos = smartAlign(h1, h2, c(0, max(corr_mat, na.rm = TRUE)))
	h = (pos[, 1] + pos[, 2])/2 - 0.1
	segments(0.91, x, 1, h)
	text(1.01, h, labels, adj = c(0, 0.5), cex = 0.7)
	w = max(strwidth(labels, cex = 0.7))
	w2 = jaccard_union[labels, k]/max_u*w
	rect(1.01, h - text_height*0.25, 1.01+w*(jaccard_gr1[labels,k]/jaccard_union[labels,k]), 
		h + text_height*0.25, col = "#FF000040")
	rect(1.01+w, h - text_height*0.25, 1.01+w*(1-jaccard_gr2[labels,k]/jaccard_union[labels,k]), 
		h + text_height*0.25, col = "#0000FF40")
	segments(1.01,  h - text_height*0.35,1.01+w2, h - text_height*0.35, lwd = 4, col = "red", lend = "butt")
	par(xpd = FALSE)
}
saveGIF(
	for(k in 1:10) {
	# par(mfrow = c(1, 2))
	foo_plot(k, !grepl("meth|DNase", rownames(corr_mat)))
	# foo_plot(k, grepl("meth|DNase", rownames(corr_mat)) & rownames(corr_mat) != "hyper_meth_vs_hypo_meth")
}, movie.name = qq("@{OUTPUT_DIR}/plots/epi_signal_global_jaccard_down_up.gif"), interval = 1, ani.width = 400, ani.height = 400)

pdf(qq("@{OUTPUT_DIR}/plots/chipseq_signal_global_jaccard_down_up.pdf"), width = 5, height = 5)
# par(mfrow = c(1, 2))
foo_plot(k = 10, !grepl("meth|DNase", rownames(corr_mat)))
# foo_plot(k = 10, grepl("meth|DNase", rownames(corr_mat)) & rownames(corr_mat) != "hyper_meth_vs_hypo_meth")
dev.off()





# cutoff_meth = quantile(abs(meth_gr$meth_diff), c(seq(0, 0.9, by = 0.1), c(0.925, 0.95, 0.975, 0.99, 0.995)))
cutoff_list = lapply(gr_list, function(gr) {
	gr = gr[abs(gr$diff) > 1e-6]
	quantile(abs(gr$diff), c(seq(0, 0.9, by = 0.1), c(0.925, 0.95, 0.975, 0.99, 0.995)))
})

foo_plot = function(i) {
	par(mfrow = c(2, 3))
	# l = abs(meth_gr$meth_diff) > cutoff_meth[i]
	# q = 1-sum(l)/length(l)
	# plot(density(meth_gr$meth_diff[l]), main = qq("methylation, |diff| > @{sprintf('%.2f', cutoff_meth[i])}, q=@{sprintf('%.4f', q)}\n@{sum(l)} (@{sprintf('%.1f', sum(l)/length(l)*100)}%) CpG sites"), xlim = c(-1, 1))
	# abline(v = 0, col = "grey", lty = 2)
	for(j in seq_along(gr_list)) {
		if(grepl("meth|DNase", gr_name[j])) next
		gr = gr_list[[j]]
		gr = gr[abs(gr$den_diff) > 1e-6]
		rg = range(gr$den_diff)
		l = abs(gr$den_diff) > cutoff_list[[j]][i]
		q = 1-sum(l)/length(l)
		par(mar = c(5.1, 4.1, 4.1, 2.1), new = FALSE)
		plot(density(gr$den_diff[l]), main = qq("@{gr_name[j]}, |rel_diff| > @{sprintf('%.2f', cutoff_list[[j]][i])}, q=@{sprintf('%.4f', q)}\n@{sum(l)} (@{sprintf('%.1f', sum(l)/length(l)*100)}%) 1kb windows"), xlim = rg)
		abline(v = 0, col = "grey", lty = 2)

		gr = gr[l]
		gr_neg = gr[gr$den_diff < 0]
		gr_pos = gr[gr$den_diff > 0]
		mtch1 = as.matrix(findOverlaps(gr_neg, gr_meth))
		x1 = gr_meth[unique(mtch1[, 2])]$meth_diff; if(sum(gr$den_diff < 0)/length(gr) < 0.1) x1 = NULL
		mtch2 = as.matrix(findOverlaps(gr_pos, gr_meth))
		x2 = gr_meth[unique(mtch2[, 2])]$meth_diff; if(sum(gr$den_diff > 0)/length(gr) < 0.1) x2 = NULL
		par(mar = c(18, 22, 4.1, 2.1), new = TRUE)
		boxplot(list(hypo = x1, hyper = x2), outline = FALSE, ylab = "methylatioh difference")
		abline(h = 0, lty = 2, col = "grey")
		if(length(x1)) qqcat("@{gr_name[j]}: hypo: [@{quantile(x1, 0.5)}, @{IQR(x1)}]\n")
		if(length(x2)) qqcat("@{gr_name[j]}: hyper: [@{quantile(x2, 0.5)}, @{IQR(x2)}]\n")
	}
}

saveGIF(
for(i in seq_along(cutoff_list[[1]])) {
	foo_plot(i)
}, movie.name = qq("@{OUTPUT_DIR}/plots/epi_signals_diff.gif"), interval = 1, ani.width = 900, ani.height = 600)

pdf(qq("@{OUTPUT_DIR}/plots/epi_signals_diff.pdf"), width = 12, height = 8)
foo_plot(10)
dev.off

gr_name = names(gr_list)
gr_meth = gr_list$meth
p_neg = vector("list", length(gr_list))
names(p_neg) = names(gr_list)
m_neg_three = p_neg
m_pos_three = p_neg
for(i in seq_along(gr_list)) {
	if(grepl("meth|DNase", gr_name[i])) next
	gr = gr_list[[i]]
	gr = gr[abs(gr$diff) > 1e-6]
	for(q_cutoff in 0:9*0.1) {
		qqcat("@{gr_name[i]}, q@{q_cutoff}\n")
		q = quantile(abs(gr$diff), q_cutoff)

		gr_neg = gr[gr$diff < -q]
		gr_pos = gr[gr$diff > q]
		mtch1 = as.matrix(findOverlaps(gr_neg, gr_meth))
		x1 = gr_meth[unique(mtch1[, 2])]$diff
		mtch2 = as.matrix(findOverlaps(gr_pos, gr_meth))
		x2 = gr_meth[unique(mtch2[, 2])]$diff
		p_neg[[i]][as.character(q_cutoff)] = sum(as.numeric(width(gr_neg)))/(sum(as.numeric(width(gr_neg))) + sum(as.numeric(width(gr_pos))))
		m_neg_three[[i]][[as.character(q_cutoff)]] = quantile(x1, c(0.1, 0.5, 0.9))
		m_pos_three[[i]][[as.character(q_cutoff)]] = quantile(x2, c(0.1, 0.5, 0.9))
	}
}

add_prop = function(x) {
	n = length(x)
	pushViewport(viewport(xscale = c(0, n), yscale = c(0, 1)))
	grid.rect(1:n-0.5, rep(0, n), width = 1, height = x, default.units = "native", gp = gpar(fill = "darkgreen"), just = "bottom")
	grid.rect(1:n-0.5, rep(1, n), width = 1, height = 1-x, default.units = "native", gp = gpar(fill = "red"), just = "top")
	upViewport()
}

add_boxplot = function(neg_lt, pos_lt) {
	n = length(neg_lt)
	neg_mat = do.call("cbind", neg_lt)
	pos_mat = do.call("cbind", pos_lt)
	rg = range(c(neg_mat, pos_mat))
	rg = max(abs(rg))
	pushViewport(viewport(xscale = c(0, n), yscale = c(-rg, rg)))
	grid.rect(1:n - 0.75, neg_mat[2, ], width = 0.5, height = neg_mat[3, ] - neg_mat[1, ], default.units = "native", gp = gpar(fill = "darkgreen"))
	grid.rect(1:n - 0.25, pos_mat[2, ], width = 0.5, height = pos_mat[3, ] - pos_mat[1, ], default.units = "native", gp = gpar(fill = "red"))
	upViewport()
}

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 6, nc = 1)))
for(i in 1:6) {
	pushViewport(viewport(layout.pos.row = i))
		pushViewport(viewport(layout = grid.layout(nr = 2, nc = 1)))
			pushViewport(viewport(layout.pos.row = 1))
			add_prop(p_neg[[i]])
			upViewport()
			pushViewport(viewport(layout.pos.row = 2))
			add_boxplot(m_neg_three[[i]], m_pos_three[[i]])
			upViewport()
		upViewport()
	upViewport()
}
upViewport()

# gr_meth = get_mean_methylation_in_genomic_features(sample_id, chromosome = CHROMOSOME, gf_list = list(chromGr_1kb_window = chromGr_1kb_window))[[1]]
# mat = mcols(gr_meth)
# mat = as.matrix(mat[, -ncol(mat)])
# meth_diff = rowMeans(mat[, sample_id_subgroup1]) - rowMeans(mat[, sample_id_subgroup2])
# gr_meth$meth_diff = meth_diff

# library(MASS)
# for(i in seq_along(gr_list)) {
# 	x = abs(gr_list[[i]]$den_diff)
# 	q = quantile(x, c(0.3, 0.95))
# 	l1 = x > q[1] & x > 1e-6 & x < q[2]
# 	gr = gr_list[[i]][l1]
# 	l2 = abs(gr_meth$meth_diff) > 0.1 & abs(gr_meth$meth_diff) < 0.3 & as.vector(seqnames(gr_meth) == "chr19")
# 	mtch = as.matrix(findOverlaps(gr, gr_meth[l2]))
# 	x = gr$den_diff[mtch[, 1]]
# 	y = gr_meth$meth_diff[l2][mtch[, 2]]
# 	image(kde2d(x, y, n = 200), col = rev(brewer.pal(11, "Spectral")), xlab = MARKS[i], ylab = "methylation diff", useRaster = TRUE)
# }

# for(i in seq_along(gr_list)[-length(gr_list)]) {
# 	for(j in (i+1):length(gr_list)) {
# 		x1 = abs(gr_list[[i]]$den_diff)
# 		x2 = abs(gr_list[[j]]$den_diff)
# 		q1 = quantile(x1, c(0.3, 0.95))
# 		q2 = quantile(x2, c(0.3, 0.95))
# 		l1 = x1 > q1[1] & x1 > 1e-6 & x1 < q1[2] & as.vector(seqnames(gr_list[[i]]) == "chr19")
# 		gr1 = gr_list[[i]][l1]
# 		l2 = x2 > q2[1] & x2 > 1e-6 & x2 < q2[2]
# 		gr2 = gr_list[[j]][l2]
# 		mtch = as.matrix(findOverlaps(gr1, gr2))
# 		x = gr1$den_diff[mtch[, 1]]
# 		y = gr2$den_diff[mtch[, 2]]
# 		image(kde2d(x, y, n = 200), col = rev(brewer.pal(11, "Spectral")), xlab = MARKS[i], ylab = MARKS[j], useRaster = TRUE)
# 	}
# }



# cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/01.hilbert_curve_global_difference.R --rerun")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=60:00:00,mem=10G -N hc_difference' '@{cmd}'")
# system(cmd)
