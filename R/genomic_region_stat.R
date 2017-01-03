########################################################
# this file should contain functions which visualize general 
# distribution of genomic regions 


# == title
# Basic statistics on genomic regions
#
# == param
# -gr_list a list of `GenomicRanges::GRanges`.
# -annotation a vector which contains class of samples. If it has names which correspond to ``gr_list``, the order of this vector is automatically adjusted.
# -annotation_color colors corresponding to classes of annotations
# -title title of the plot
# -species species, necessary if ``type`` is set to ``proportion``.
# -type type of statistics
# -chromosome subset of chromosomes
# -by_chr take all chromosomes as a whole or calculate statistics for every chromosome
#
# == details
# The function makes barplot to visualize different statistics in all samples.
#
# For ``type`` settings:
#
# -proportion proportion of total length of regions compared to the whole genome. It is more unbiased.
# -number number of regions. Sometimes only looking at the number of regions gives biased estimation of
#     amount of regions if the width of regions are very viarable.
# -median_width median width of regions
#
# == value
# A data frame which contains statistics for each chromosome in each sample.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# set.seed(123)
# gr_list = lapply(1:10, function(i) {
# 	df = generateRandomBed(1000)[1:sample(100:1000, 1), ]
# 	GRanges(df[[1]], ranges = IRanges(df[[2]], df[[3]]))
# })
# names(gr_list) = paste0("sample_", 1:10)
# basic_genomic_regions_stat(gr_list)
# basic_genomic_regions_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
#     annotation_color = c("a" = "red", "b" = "blue"))
# basic_genomic_regions_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
#     annotation_color = c("a" = "red", "b" = "blue"), type = "number")
# basic_genomic_regions_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
#     annotation_color = c("a" = "red", "b" = "blue"), type = "median_width")
# basic_genomic_regions_stat(gr_list, annotation = rep(letters[1:2], each = 5), 
#     annotation_color = c("a" = "red", "b" = "blue"), by_chr = TRUE)
basic_genomic_regions_stat = function(gr_list, annotation = NULL, annotation_color = NULL, 
	title = paste0("Basic statistics for genomic regions"), species = "hg19", 
	type = c("proportion", "number", "median_width"),
	chromosome = paste0("chr", c(1:22, "X", "Y")), by_chr = FALSE) {

	type = match.arg(type)[1]
	chromInfo = read.chromInfo(species = species, chromosome.index = chromosome)
	
	if(inherits(gr_list, "GRanges")) {
		gr_list = list(gr_list)
		names(gr_list) = deparse(substitute(gr))
	}
	gr_list = lapply(gr_list, function(gr) {
		gr[seqnames(gr) %in% chromosome]
	})

	if(is.null(names(gr_list))) {
		names(gr_list) = seq_along(gr_list)
	}
	sid = names(gr_list)

    
	if(!is.null(annotation)) {
		if(is.null(names(annotation))) names(annotation) = sid
		annotation = annotation[sid]
	}

	n_gr = length(gr_list)
	if(by_chr) {
		n_chr = length(chromInfo$chromosome)
		if(type %in% "proportion") {
			# proportion in genome
			lt = lapply(chromInfo$chromosome, function(chr) {
					sapply(gr_list, function(gr) {
						gr = gr[ seqnames(gr) == chr ]
						if(length(gr) == 0) {
							return(0)
						} else {
							sum(width(gr))/chromInfo$chr.len[chromInfo$chromosome == chr]
						}
					})
				})
			x = unlist(lt)
			ylab = "proportion in genome"

		} else if(type %in% "number") {
			# number of regions
			# proportion in genome
			lt = lapply(chromInfo$chromosome, function(chr) {
					sapply(gr_list, function(gr) {
						gr = gr[ seqnames(gr) == chr ]
						length(gr)
					})	
				})
			x = unlist(lt)
			ylab = "numbers of regions"

		} else if(type %in% "median_width") {
			# median length
			# proportion in genome
			lt = lapply(chromInfo$chromosome, function(chr) {
					sapply(gr_list, function(gr) {
						gr = gr[ seqnames(gr) == chr ]
						if(length(gr) == 0) {
							return(0)
						} else {
							median(width(gr))
						}
					})	
				})
			x = unlist(lt)
			ylab = "median width of regions"

		}	
	} else {
		if(type %in% "proportion") {
			# proportion in genome
			genome.len = sum(chromInfo$chr.len)
			x = sapply(gr_list, function(gr) sum(width(gr))/genome.len)
			ylab = "proportion in genome"

		} else if(type %in% "number") {
			# number of regions
			x = sapply(gr_list, length)
			ylab = "numbers of regions"

		} else if(type %in% "median_width") {
			# median length
			x = sapply(gr_list, function(gr) median(width(gr)))
			ylab = "median width of regions"
			
		}
	}

	if(by_chr) {   # use lattice plot
		df = data.frame(x = x, 
			sid = rep(sid, times = length(chromInfo$chromosome)), 
			chr = rep(chromInfo$chromosome, each = length(sid)),
			stringsAsFactors = FALSE)
		df$chr = factor(df$chr, levels = chromInfo$chromosome)
		df$sid = factor(df$sid, levels = sid)	
	} else {   # use barplot
		sid = names(x)
		df = data.frame(x = x, sid = factor(sid, levels = sid))
	}

	if(is.null(annotation)) {
		gp = ggplot(df, aes(y = x, x = sid)) + geom_bar(stat = "identity")
	} else {
		df = cbind(df, annotation = annotation[df$sid])
		gp = ggplot(df, aes(y = x, x = sid)) + geom_bar(stat = "identity", aes(fill = annotation))
	}

	if(by_chr) gp = gp + facet_wrap( ~ chr)

	if(!is.null(annotation_color)) {	
		gp = gp + scale_fill_manual(values = annotation_color)
	}

	gp = gp + xlab("") + ylab(ylab) + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + ggtitle(title)
	print(gp)

	return(invisible(df))
}
