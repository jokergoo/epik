
# == title
# Add unit "bp" to the number
#
# == param
# -x a numeric vector. 
#
# == details
# It adds a new ``bp`` class to the vector.
#
# == value
# A same vector as input
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# bp(10)
# bp(10.1)
bp = function(x) {
	x = round(x)
	class(x) = c(class(x), "bp")
	x
}

# == title
# Add unit "kb" to the number
#
# == param
# -x a numeric vector.
#
# == details
# The input values are multiplied by 1000 and send to `bp`.
#
# == value
# A numeric vector measured in bp
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# kb(10)
# kb(10.01)
kb = function(x) {
	bp(x*1000)
}

# == title
# Add unit "mb" to the number
#
# == param
# -x a numeric vector.
#
# == details
# The input values are multiplied by 1000000 and send to `bp`.
#
# == value
# A numeric vector measured in bp
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# mb(10)
# mb(10.01)
mb = function(x) {
	bp(x*1000000)
}

print.bp = function(x, ...) {
	x = paste0(prettyNum(sprintf("%.0f", x), big.mark = ","), " bp")
	print(x, quote = FALSE)
}


# == title
# Merge genomic regions
#
# == param
# -gr a `GenomicRanges::GRanges` object
# -gap a numeric value means to extend each region by ``gap/2`` times of its width before merging. If ``gap`` represents
#      number of base pairs, use `bp`, `kb` or `mb` to wrap it. If ``gap`` represents absolute number of base pairs,
#      the functionality is same as `GenomicRanges::reduce` (``gap`` is sent to ``min.gapwidth``).
# -max_gap maximum distance to merge, measured in base pairs. Only work if ``gap`` is a ratio value.
# -.message internal use
# -.revmap internal use
# -... further arguments passed to `GenomicRanges::reduce` (exclude ``min.gapwidth`` and ``with.revmap``)
#
# == details
# `GenomicRanges::reduce` only merges regions with fixed gap width, but sometimes it is not reasonable to set gap
# to a same width for all regions. Assuming we have a list of differentially methylated regions (DMRs) and we want to reduce
# the number of DMRs by merging neighouring DMRs. DMRs distribute differently in different places in the genome, e.g. DMRs are dense
# and short in CpG-rich regions (e.g. CpG islands) while long in CpG-sparse regions (e.g. gene bodies and intergenic regions),
# thus the merging should be applied based to the width of every DMR itself. `reduce2` can merge regions by the width of every region itself.
# This type of merging is dynamic because after each iteration of merging, some regions are merged into a large region and 
# it will has longer extension. The whole merging will proceed iteratively unless there is no new merging.
#
# Note ``with.revmap`` is always set to ``TRUE`` when calling `GenomicRanges::reduce`, thus there is always a ``revmap``
# meta column in the returned `GenomicRanges::GRanges` object.
#
# == value
# a `GenomicRanges::GRanges` object with a ``revmap`` meta column.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 4, 8 ,16), 
#     end = c(2, 5, 10, 30)), value = 1:4)
# reduce2(gr, gap = bp(2))
# reduce2(gr, gap = 0.6)
# reduce2(gr, gap = 0.6, max_gap = 4)
reduce2 = function(gr, gap = 0.1, max_gap = Inf, .message = TRUE, .revmap = NULL, ...) {
	if(length(gr) == 0) {
		mcols(gr) = NULL
		gr$revmap = IntegerList()
		return(gr)
	}

	if(inherits(gap, "bp")) {
		message(qq("Regions are merged by an absolute distance: (<= @{gap} bp)"))
		return(reduce(gr, min.gapwidth = gap, with.revmap = TRUE, ...))
		
	} else {

		if(.message) message(qq("Regions are extended by @{gap}*width (maximum @{max_gap} bp)."))
		
		n = length(gr)
		gr_extend = gr
		gr_width = width(gr)
		
		start(gr_extend) = start(gr) - round(pmin(gr_width*gap/2, max_gap/2))
		end(gr_extend) = end(gr) + round(pmin(gr_width*gap/2, max_gap/2))

		gr_reduced = reduce(gr_extend, with.revmap = TRUE, ...)
		revmap = mcols(gr_reduced)[, "revmap"]
		s = start(gr)
		e = end(gr)
		os = sapply(revmap, function(ind) min(s[ind]))
		oe = sapply(revmap, function(ind) max(e[ind]))
		start(gr_reduced) = os
		end(gr_reduced) = oe

		if(is.null(.revmap)) {
			.revmap = lapply(1:n, function(x) x)
		}

		.revmap = lapply(revmap, function(ind) unlist(.revmap[ind]))

		n2 = length(gr_reduced)
		if(n2 == 1) {
			gr_reduced$revmap = as(.revmap, "IntegerList")
			return(gr_reduced)
		} else if(n == n2) {
			gr_reduced$revmap = as(.revmap, "IntegerList")
			return(gr_reduced)
		} else {
			message(qq("regions have been reduced from @{n} to @{n2}..."))
			gr = reduce2(gr_reduced, max_gap = max_gap, gap = gap, .message = FALSE, .revmap = .revmap, ...)
			return(gr)
		}
	}
}
