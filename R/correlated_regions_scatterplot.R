
# == title
# Scatter plot between methylation and expression in one correlated region
#
# == param
# -cr correlated regions
# -expr the expression matrix which was used in `correlated_regions`
# -gi gene id
# -text_column the column name in ``cr`` which will be plotted as text annotations.
# -xlab xlab in the plot
# -ylab ylab in the plot
#
# == details
# Scatterplot for all correlated regions corresponding to the gene will be made. If you want to make
# a subset of correlated regions, directly subset ``cr``.
#
# Internally it uses `scatterplot_with_boxplot`.
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
cr_scatterplot = function(cr, expr, gi = NULL, text_column, xlab = "Methylation", ylab = "Expression") {

	cr_param = metadata(cr)$cr_param
	
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
			text_list = unlist(as.data.frame(mcols(cr_now)[, text_column, drop = FALSE])))
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
# -subgroup groups of data points
# -subgroup_col colors for groups
# -main title for the plot
# -xlab labels on x-axis
# -ylab labels on y-axis
# -xlim range on x-axis
# -ylim range on y-axis
# -text_list additional text which is a named vector or list (if the text is mixed with character and numbers)
#
# == details
# Boxplots are on the left and bottom, the scatter plot is on the top right.
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# scatterplot_with_boxplot(rnorm(40), rnorm(40), subgroup = sample(letters[1:2], 40, replace = TRUE))
#
scatterplot_with_boxplot = function(x, y, subgroup = rep("unknown", length(x)), 
	subgroup_col = structure(seq_along(levels(subgroup)), names = levels(subgroup)),
	main = NULL, xlab = NULL, ylab = NULL, xlim = range(x), ylim = range(y), text_list = NULL) {

	if(!is.factor(subgroup)) {
		subgroup = factor(subgroup)
	}
	xrange = xlim
	yrange = ylim

	layout(rbind(1:2, 3:4), widths=c(1, 2), heights=c(2, 1))
    par(mar = c(0, 5, 5, 0))
    boxplot(y ~ subgroup, ylim = yrange, axes = FALSE, ann = FALSE, col = subgroup_col[levels(subgroup)])
    box()
    axis(side = 2, cex.axis = 1.5)
    title(ylab = ylab, cex.lab = 1.5)

    par(mar = c(0, 0, 5, 5))
    plot(x, y, xlim = xrange, ylim = yrange, axes = FALSE, ann = FALSE, cex = 1.5, pch = 16, col = subgroup_col[subgroup])
    box()
    title(main, cex.main = 1)

   # lm.res = lm(y ~ x)
   # abline(lm.res, col = "#E41A1C")

    par(mar = c(5, 5, 0, 0), xpd = NA)
    plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, ann = FALSE)
    legend("bottomleft", legend = levels(subgroup), pch = 16, col = subgroup_col[levels(subgroup)])
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
    boxplot(x ~ subgroup, ylim = xrange, horizontal = TRUE, axes = FALSE, ann = FALSE, col = subgroup_col[levels(subgroup)])
    box()
    axis(side = 1, cex.axis = 1.5)
    #axis(side = 4, at = seq_along(levels(d$cate)), labels=levels(d$cate))
    title(xlab = xlab, cex.lab = 1.5)
}
