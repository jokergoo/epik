
# == title
# Scatter plot between methylation and expression in a correlated region
#
# == param
# -cr correlated regions from `correlated_regions` or `filter_correlated_regions`
# -expr the expression matrix which is same as in `correlated_regions`
# -gi gene id
# -text_column which column in ``cr`` should be put as annotation text in the plot
# -xlab xlab in the plot
# -ylab ylab in the plot
#
# == details
# Scatterplot for all CRs corresponding to the gene will be made.
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
cr_scatterplot = function(cr, expr, gi = NULL, text_column,
	xlab = "Methylation", ylab = "Expression") {

	cr_param = metadata(cr)$cr_param
	
	annotation = attr(cr, "factor")
	annotation_color = attr(cr, "col")
	sample_id = attr(cr, "sample_id")

	subgroup = cr_param$subgroup
	subgroup_col = cr_param$col
	sample_id = cr_param$sample_id

	if(!is.null(gi)) {
		cr = cr[cr$gene_id %in% gi]
	}

	mmat = mcols(cr)
	mmat = as.matrix(mmat[, paste0("mean_meth_", sample_id)])
	colnames(mmat) = sample_id

	if(missing(text_column)) {
		text_column = colnames(mcols(cr))
		text_column = text_column[!grepl("mean_meth_", text_column)]
		text_column = intersect(text_column, c("corr_p", "meth_IQR", "meth_anova", "meth_diameter", "gene_tss_dist"))
	}

	for(i in seq_len(length(cr))) {
		cr_now = cr[i]
		gi = cr_now$gene_id
		chr = as.vector(seqnames(cr_now))
		
		v = mmat[i, ]
		e = expr[gi, colnames(mmat)]

		scatterplot_with_boxplot(v, e, subgroup, subgroup_col, 
			main = qq("@{chr}:@{start(cr_now)}-@{end(cr_now)}, @{gi},\ncor = @{sprintf('%.2f', cr_now$corr)}, n_CpG = @{cr_now$ncpg}"),
			xlab = xlab, ylab = ylab,
			text_list = unlist(mcols(cr_now)[, text_column, drop = FALSE]))
		if(i %% 50 == 0) {
			message(qq("@{i}/@{length(cr)} CRs finished"))
		}
	}
}

# == title
# Scatterplot with boxplots on both sides
#
# == param
# -x values on x-axis
# -y values on y-axis
# -annotation annotations which show groups of data points
# -annotation_color colors for annotation
# -main title for the plot
# -xlab labels on x-axis
# -ylab labels on y-axis
# -xlim range on x-axis
# -ylim range on y-axis
# -text_list additional text which is a named vector or list (if the text is mixed with character and numbers)
#
# == details
# On the left and bottom, there are boxplots and on the top right, there is the scatter plot.
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
scatterplot_with_boxplot = function(x, y, annotation = rep("unknown", length(x)), 
	annotation_color = structure(seq_along(levels(annotation)), names = levels(annotation)),
	main = NULL, xlab = NULL, ylab = NULL, xlim = range(x), ylim = range(y), text_list = NULL) {

	if(!is.factor(annotation)) {
		annotation = factor(annotation)
	}
	xrange = xlim
	yrange = ylim

	layout(rbind(1:2, 3:4), widths=c(1, 2), heights=c(2, 1))
    par(mar = c(0, 5, 5, 0))
    boxplot(y ~ annotation, ylim = yrange, axes = FALSE, ann = FALSE, col = annotation_color[levels(annotation)])
    box()
    axis(side = 2, cex.axis = 1.5)
    title(ylab = ylab, cex.lab = 1.5)

    par(mar = c(0, 0, 5, 5))
    plot(x, y, xlim = xrange, ylim = yrange, axes = FALSE, ann = FALSE, cex = 1.5, pch = 16, col = annotation_color[annotation])
    box()
    title(main, cex.main = 1)

   # lm.res = lm(y ~ x)
   # abline(lm.res, col = "#E41A1C")

    par(mar = c(5, 5, 0, 0), xpd = NA)
    plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, ann = FALSE)
    legend("bottomleft", legend = levels(annotation), pch = 16, col = annotation_color[levels(annotation)])
    if(length(text_list)) {
	    text_list_name = names(text_list)
	    text_list = as.character(text_list)
	    for(i in seq_along(text_list)) {
	    	if(!is.na(suppressWarnings(as.numeric(text_list[i])))) {
	    		if(nchar(text_list[i]) > 5) {
	    			text_list[i] = sprintf("%.2e", as.numeric(text_list[i]))
	    		}
	    	}
	    }
	    text(0, 1, qq("@{text_list_name} = @{text_list}\n"), adj = c(0, 1))	
	}

    par(mar = c(5, 0, 0, 5))
    boxplot(x ~ annotation, ylim = xrange, horizontal = TRUE, axes = FALSE, ann = FALSE, col = annotation_color[levels(annotation)])
    box()
    axis(side = 1, cex.axis = 1.5)
    #axis(side = 4, at = seq_along(levels(d$cate)), labels=levels(d$cate))
    title(xlab = xlab, cex.lab = 1.5)
}
