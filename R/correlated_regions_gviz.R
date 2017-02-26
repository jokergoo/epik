
if(!exists(".boxes")) {
	.boxes = Gviz.epik:::.boxes
}
if(!exists(".arrowBar")) {
	.arrowBar = Gviz.epik:::.arrowBar
}
if(!exists(".fontGp")) {
	.fontGp = Gviz.epik:::.fontGp
}

elementNROWS = function (x) {
    if (!is.list(x))
        x <- as.list(x)
    ans <- try(.Call2("sapply_NROW", x, PACKAGE = "S4Vectors"),
        silent = TRUE)
    if (!inherits(ans, "try-error")) {
        names(ans) <- names(x)
        return(ans)
    }
    return(vapply(x, NROW, integer(1)))
}

# == title
# Customized Gviz plot for a single gene
#
# == param
# -sig_cr correlated regions which show significant correlations
# -gi gene id
# -expr the expression matrix which was used in `correlated_regions`
# -txdb the transcriptome annotation which was used in `correlated_regions`
# -gf_list a list of `GenomicRanges::GRanges` objects which contains additional genomic annotations
# -hm_list a list of `GenomicRanges::GRanges` objects which contains histome modification peaks. The value is a two-layer
#    list. The first layer is histome modifications and the second layer is the peaks in each sample which has current histome
#    modification data. Name of the first layer is histome mark name and the name of the second layer is sample ID.
# -title title of the plot
#
# == details
# There are following Gviz tracks:
#
# - gene models. Multiple transcripts will also be plotted.
# - correlation between methylation and expression
# - heatmap for methylation
# - significant correlated regions
# - CpG density
# - annotation to other genomic features, if provided
# - histome modification data, if provided
#
# A modified version of Gviz (https://github.com/jokergoo/Gviz.epik ) is used to make the plot.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
cr_gviz = function(sig_cr, gi, expr, txdb, gf_list = NULL, hm_list = NULL, title = gi) {

	cr_param = metadata(sig_cr)$cr_param
	sample_id = cr_param$sample_id
	extend = cr_param$extend
	window_size = cr_param$window_size
	window_step = cr_param$window_step
	max_width = cr_param$max_width
	cor_method = cr_param$cor_method
	subgroup = cr_param$subgroup
	cov_filter = cr_param$cov_filter
	raw_meth = cr_param$raw_meth
	cov_cutoff = cr_param$cov_cutoff
	min_dp = cr_param$min_dp
	species = cr_param$species

	if(is.null(raw_meth)) raw_meth = FALSE
	if(is.null(cov_cutoff)) cov_cutoff = 0
	if(is.null(min_dp)) min_dp = 5
	if(!raw_meth) cov_cutoff = 0
	
	if(!gi %in% sig_cr$gene_id) {
		stop(paste0("cannot find ", gi, "in cr.\n"))
	}

	chr = as.vector(seqnames(sig_cr[sig_cr$gene_id == gi]))[1]
	sig_cr = sig_cr[sig_cr$gene_id == gi]

	methylation_hooks$set_chr(chr, verbose = FALSE)
	
	e = expr[gi, sample_id]

	gene = genes(txdb)
	gene_start = start(gene[gi])
	gene_end = end(gene[gi])

	gene_start = gene_start - extend
	gene_start = ifelse(gene_start > 0, gene_start, 1)
	gene_end = gene_end + extend

	site = start(methylation_hooks$gr)

	gm_site_index = extract_sites(gene_start, gene_end, site, TRUE, 0)
	gm_site = site[gm_site_index]
	gm_meth = methylation_hooks$meth[gm_site_index, sample_id, drop = FALSE]
	gm_cov = methylation_hooks$cov[gm_site_index, sample_id, drop = FALSE]

	if(!is.null(cov_filter)) {
		l = apply(gm_cov, 1, cov_filter)
		gm_site = gm_site[l]
		gm_meth = gm_meth[l, , drop = FALSE]
		gm_cov = gm_cov[l, , drop = FALSE]
	}

	message(qq("rescan on @{gi} to calculate correlation in @{window_size}-CpG window..."))
	
	gr = correlated_regions_by_window(gm_site, gm_meth, e, gm_cov, cov_cutoff = cov_cutoff, chr = chr,
			subgroup = subgroup, cor_method = cor_method, window_size = window_size, window_step = window_step, min_dp = min_dp,
			max_width = max_width)

	message("add transcripts to gviz tracks...")
	options(Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/")
	trackList = list()
	trackList = pushTrackList(trackList, GenomeAxisTrack())
	trackList = pushTrackList(trackList, IdeogramTrack(genome = species, chromosome = chr))
	grtrack = GeneRegionTrack(txdb, chromosome = chr, start = gene_start, end = gene_end, 
		name="Gene\nmodel", showId = TRUE, rotate.title = TRUE, col = NA, showTitle = FALSE,
		size = 0.5)
		
	.boxes_wrap = function(GdObject, offsets) {
		df = .boxes(GdObject, offsets)
		l = df$gene == gi
		# df$fill[l] = "pink"
		df$fill[!l] = paste0(df$fill[!l], "40")
		df
	}
	assignInNamespace(".boxes", .boxes_wrap, "Gviz.epik")
	

	tx_gene_mapping = structure(grtrack@range$gene, names = grtrack@range$transcript)

	.arrowBar_wrap = function(xx1, xx2, strand, coords, y=20, W=3, D=10, H, col, lwd, lty, alpha, barOnly=FALSE,
        diff=.pxResolution(coord="y"), min.height=3) {
		env = parent.frame()
		if("bar" %in% ls(envir = env)) {
			bar = get("bar", envir = env)
			arrow_col = ifelse(tx_gene_mapping[rownames(bar)] == gi, "darkgrey", "#00000020")
		}
		.arrowBar(xx1 = xx1, xx2 = xx2, strand = strand, coords = coords, y=y, W=W, D=D, H, col = arrow_col, lwd = lwd, lty = lty, alpha = alpha, barOnly=barOnly,
        	diff=diff, min.height=min.height)
	}
	assignInNamespace(".arrowBar", .arrowBar_wrap, "Gviz.epik")
	
	.fontGp_wrap = function(GdObject, subtype = NULL, ...) {
		gp = .fontGp(GdObject, subtype, ...)
		if(!is.null(subtype)) {
			if(subtype == "group") {
				env = parent.frame()
				if("bartext" %in% ls(envir = env)) {
					bartext = get("bartext", envir = env)
					tx_name = bartext$txt
					l = tx_gene_mapping[tx_name] == gi
					gp$col = ifelse(l, "#808080", "#00000020")
				}
			}
		}
		return(gp)
	}
	assignInNamespace(".fontGp", .fontGp_wrap, "Gviz.epik")
	

	trackList = pushTrackList(trackList, grtrack)

	## correlation track
	message("add correlation line to gviz tracks...")
	corr_mat = matrix(0, nrow = 2, ncol = length(gr$corr))
	corr_mat[1, gr$corr > 0] = gr$corr[gr$corr > 0]
	corr_mat[2, gr$corr < 0] = gr$corr[gr$corr < 0]
	trackList = pushTrackList(trackList, DataTrack(name = qq("Correlation\nCpG window = @{window_size}"),
								range = gr,
								genome = species,
								data = corr_mat,
								type = c("hist"),
								groups = c("pos", "neg"),
								fill.histogram = c("red", "green"),
								col.histogram = NA,
								ylim = c(-1, 1), legend = FALSE,
								panelFun = local({window_size = window_size; function() grid.text(qq("Correlation, CpG window = @{window_size}bp"), 0, unit(1, "npc") - unit(2, "mm"), gp = gpar(fontsize = 10), just = c("left", "top"))}),
								size = 1.5))

	message("add sig_cr to gviz tracks...")
	pos_cr = sig_cr[sig_cr$corr > 0]
	if(length(pos_cr))
		trackList = pushTrackList(trackList, constructAnnotationTrack(reduce(pos_cr), chr, name = "sig_pos_cr", fill = "red", col = NA, 
			rotate.title = TRUE, start = gene_start, end = gene_end, size = 0.5,
			panelFun = function() {grid.text("pos_cr", 0, unit(0.5, "npc"), gp = gpar(fontsize = 10), just = c("left", "center"))}))
	neg_cr = sig_cr[sig_cr$corr < 0]
	if(length(neg_cr))
		trackList = pushTrackList(trackList, constructAnnotationTrack(reduce(neg_cr), chr, name = "sig_neg_cr", fill = "darkgreen", col = NA, 
			rotate.title = TRUE, start = gene_start, end = gene_end, size = 0.5,
			panelFun = function() {grid.text("neg_cr", 0, unit(0.5, "npc"), gp = gpar(fontsize = 10), just = c("left", "center"))}))

	message("add methylation to gviz tracks...")
	meth_mat = as.matrix(mcols(gr)[, paste0("mean_meth_", sample_id)])
	colnames(meth_mat) = NULL
	if(is.null(subgroup)) {
		trackList = pushTrackList(trackList, DataTrack(name = "meth",
										start = start(gr),
										end = end(gr),
										chromosome = seqnames(gr),
										genome = species,
										data = t(meth_mat),
										type = "heatmap",
										showSampleNames = FALSE,
										gradient = c("blue", "white", "red"),
										size = 0.2*ncol(meth_mat),
										col = NA,
										panelFun = function() {grid.text("methylation", 0, unit(1, "npc") - unit(2, "mm"), gp = gpar(fontsize = 10), just = c("left", "top"))},))
	} else {
		for(t in unique(subgroup)) {
			mat = meth_mat[, subgroup == t]
			trackList = pushTrackList(trackList, DataTrack(name = t,
										start = start(gr),
										end = end(gr),
										chromosome = seqnames(gr),
										genome = species,
										data = t(mat),
										type = "heatmap",
										showSampleNames = FALSE,
										gradient = c("blue", "white", "red"),
										size = 0.2*ncol(mat),
										col = NA,
										panelFun = local({t = t; function() grid.text(qq("methylation, @{t}"), 0, unit(1, "npc") - unit(2, "mm"), gp = gpar(fontsize = 10), just = c("left", "top"))})
									))
		}
	}
	
	### CpG density per 100bp
	message("add cpg density to gviz tracks...")
	segment = seq(gm_site[1], gm_site[length(gm_site)], by = 100)
	start = segment[-length(segment)]
	end = segment[-1]-1
	num = sapply(seq_along(start), function(i) sum(gm_site >= start[i] & gm_site <= end[i]))
	trackList = pushTrackList(trackList, DataTrack(name = "#CpG\nper 100bp",
		                            start = start,
		                            end = end,
		                            chromosome = rep(chr, length(start)),
									genome = species,
									data = num,
									col = "black",
									type = "hist",
									rotate.title = TRUE,
									size = 1,
									col.histogram = "orange",
									fill = "orange",
									panelFun = function() {grid.text("CpG density, window = 100bp", 0, unit(1, "npc") - unit(2, "mm"), gp = gpar(fontsize = 10), just = c("left", "top"))},))
	
	message("add other genomic features to gviz tracks...")
	gf_name = names(gf_list)
	for(i in seq_along(gf_list)) {
		trackList = pushTrackList(trackList, constructAnnotationTrack(gf_list[[i]], chr, name = gf_name[i], rotate.title = TRUE, start = gene_start, end = gene_end, size = 0.5,
			panelFun = local({gf_name = gf_name[i]; function() grid.text(gf_name, 0, unit(0.5, "npc"), gp = gpar(fontsize = 10), just = c("left", "center"))})))
	}

	# show mean signal in each subgroup
	if(!is.null(hm_list)) {
		# hm_list is a list of list
		# mark->sid->gr
		all_colors = brewer.pal(8, "Set2")
		hm_name = names(hm_list)
		for(i in seq_along(hm_list)) {
			single_hm_list = hm_list[[i]]
		
			message(qq("add histome modification (@{hm_name[i]}) to gviz tracks..."))
			single_hm_list2 = lapply(single_hm_list, function(gr) {
				gr = gr[seqnames(gr) == chr]
				l = start(gr) > gene_start & end(gr) < gene_end
				gr[l]
			})

			hm_merged = NULL
			for(j in seq_along(single_hm_list2)) {
				if(length(single_hm_list2[[j]])) hm_merged = rbind(hm_merged, as.data.frame(single_hm_list2[[j]]))
			}
			hm_merged = GRanges(seqnames = hm_merged[[1]], ranges = IRanges(hm_merged[[2]], hm_merged[[3]]))
			if(length(hm_merged) > 0) {
				segments = as(coverage(hm_merged), "GRanges")[-1]
				# also add zero-coverage to the GRanges object
				gr_g = GRanges(seqnames = chr, ranges = IRanges(gene_start, gene_end))
				gr_diff = GRanges(seqnames = chr, ranges = setdiff(ranges(gr_g), ranges(segments)))
				if(length(gr_diff)) {
					gr_diff$score = 0
					segments = c(segments, gr_diff)
					segments = sort(segments)
				}

				# covert to matrix
				hm_mat = matrix(0, nrow = length(single_hm_list), ncol = length(segments))
				rownames(hm_mat) = names(single_hm_list)
				for(j in seq_along(single_hm_list2)) {
					mtch = as.matrix(findOverlaps(segments, single_hm_list2[[j]]))
					hm_mat[j, mtch[, 1]] = single_hm_list2[[j]][mtch[, 2]]$density
				}

				if(is.null(subgroup)) {
						mat = cbind(hm_mat, rep(0, nrow(hm_mat)))
						mat = hm_mat
						# mat[1, ncol(mat)] = max(hm_mat)
						mean_signal = colMeans(mat)
						trackList = pushTrackList(trackList, DataTrack(name = hm_name[i],
													start = start(segments),
													end = end(segments),
													chromosome = seqnames(segments),
													genome = species,
													data = mean_signal,
													type = "hist",
													size = 1,
													ylim = c(0, max(mean_signal)),
													col.histogram = all_colors[i],
													fill = all_colors[i],
													panelFun = local({hm_name = hm_name[i]; function() grid.text(hm_name, 0, unit(1, "npc") - unit(2, "mm"), gp = gpar(fontsize = 10), just = c("left", "top"))})))
				} else {
					mean_signal_list = list()
					for(t in unique(subgroup)) {
						mat = hm_mat[rownames(hm_mat) %in% sample_id[subgroup == t], , drop = FALSE]
						# mat = cbind(mat, rep(0, nrow(mat)))
						# mat[1, ncol(mat)] = max(hm_mat)
						mean_signal_list[[t]] = colMeans(mat)
					}
					ylim = c(0, max(unlist(mean_signal_list)))
					for(t in unique(subgroup)) {
						mat = hm_mat[rownames(hm_mat) %in% sample_id[subgroup == t], , drop = FALSE]
						# mat = cbind(mat, rep(0, nrow(mat)))
						# mat[1, ncol(mat)] = max(hm_mat)
						mean_signal = colMeans(mat)

						trackList = pushTrackList(trackList, DataTrack(name = qq("@{hm_name[i]}\n@{t}"),
													start = start(segments),
													end = end(segments),
													chromosome = seqnames(segments),
													genome = species,
													data = mean_signal,
													type = "hist",
													size = 1,
													ylim = ylim,
													col.histogram = all_colors[i],
						                            fill = all_colors[i],
						                            panelFun = local({hm_name = hm_name[i]; t = t; function() grid.text(qq("@{hm_name}, @{t}"), 0, unit(1, "npc") - unit(2, "mm"), gp = gpar(fontsize = 10), just = c("left", "top"))})))
					}
				}
			}
		}
	}

	message("draw gviz plot...")
	# convert fontsize to cex
	plotTracks(trackList, from = gene_start, to = gene_end, chromosome = chr, main = title, cex.main = 1, showTitle = FALSE)

	n_tx = length(unique(grtrack@range[grtrack@range$gene == gi]$transcript))
	size1 = length(hm_list)*length(unique(subgroup)) + 0.5*length(gf_list) + 1 + 0.2*length(sample_id) + 0.5 + 0.5 + 1.5
	# (1.5*n_tx + 3 + 4)*(2/3) + length(strsplit(title, "\n")[[1]]) + 1

	# the height of text with fontsize = 12 equals to 0.12 inches
	# 0.5 is the inches of one single histome mark track
	num = size1*0.5 + (7*2/3 + length(strsplit(title, "\n")[[1]]) + 1)*0.12
	coef = 1.5*2/3*0.12
	hh = coef*n_tx + num
	num = sprintf("%.2f", num)
	coef = sprintf("%.2f", coef)
	hh = sprintf("%.2f", hh)

	message(qq("The suggested height of the image is @{coef}*n_tx + @{num} inches, here n_tx = @{n_tx} and the height is @{hh} inches."))
	
	assignInNamespace(".boxes", .boxes, "Gviz.epik")
	assignInNamespace(".arrowBar", .arrowBar, "Gviz.epik")
	assignInNamespace(".fontGp", .fontGp, "Gviz.epik")

	return(invisible(NULL))
}

pushTrackList = function(trackList, track) {
	if(!is.null(track)) {
		trackList[[length(trackList) + 1]] = track
	}
	return(trackList)
}

constructAnnotationTrack = function(gr, chr, name = NULL, genome = "hg19", start = 0, end = Inf, ...) {
	gr2 = gr[seqnames(gr) == chr]
	gr2 = gr2[end(gr2) > start & start(gr2) < end]

	if(length(gr2)) {

		AnnotationTrack(name = name,
		                start = start(gr2),
		                end = end(gr2),
		                chromosome = seqnames(gr2),
		                genome = genome, 
		                stacking = "dense",
		                showTitle = TRUE, 
		                height = unit(5, "mm"),
		                ...)
	} else {
		NULL
	}
}

