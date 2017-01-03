

# == title
# Read methylation dataset
#
# == param
# -... All the arguments do not make sense here, see 'details' section.
# -RESET .
# -READ.ONLY .
# -LOCAL .
#
# == detail
# Methylation dataset from whole genome bisulfite sequencing is always huge and it does not
# make sense to read them all into the memory. Normally, the methylation dataset is stored
# by chromosome and this hook function can be set to read methylation data in a per-chromosome
# manner.
#
# Generally, for methylation dataset, there are methylation rate, CpG coverage and genomic positions
# for CpG sites. Sometimes there is also smoothed methylation rate. All these datasets can be set
# by defining a proper ``methylation_hooks$get_by_chr``. The value for ``methylation_hooks$get_by_chr``
# is a function with only one argument which is the chromosome name. This function defined how to
# read methylation dataset for a single chromosome. The function must return a list which contains
# following mandatory elements:
#
# -gr a `GenomicRanges::GRanges` object which contains genomic positions for CpG sites.
# -meth a matrix which contains methylation rate. This will be the main methylation dataset the epik
#       package uses, so it should be smoothed methylation rate if the CpG coverage is not high.
#       Note, this matrix must have column names which is sample names and will be used to match
#       other datasets (e.g. RNASeq)
#
# It can also contain some optional elements and they are not needed for the core analysis:
#
# -cov a matrix which contains CpG coverage.
# -raw a matrix which contains unsmoothed methylation rate (or the original methylation rate calculatd
#      as the fraction of methylated CpG in a CpG site)
#
# Note each row in above datasets should correspond to the same CpG site. 
#
# After ``methylation_hooks$get_by_chr`` is set, the "current chromosome" for the methylation dataset
# can be set by ``methylation_hooks$set_chr(chr)`` where ``chr`` is the chromosome name you want to go.
# After validating the dataset, following variables can be used directly:
#
# - ``methylation_hooks$gr``
# - ``methylation_hooks$meth``
# - ``methylation_hooks$sample_id``
# - ``methylation_hooks$cov`` if available
# - ``methylation_hooks$raw`` if available
#
# == value
# Hook functions
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
methylation_hooks = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}

METH_OBJ = list(meth = NULL, raw = NULL, cov = NULL, gr = NULL, chr = NULL)
methylation_hooks = setGlobalOptions(
	get_by_chr = list(.value = NULL, .class = "function"),
	meth = list(.value = function() METH_OBJ$meth, .private = TRUE, .visible = FALSE),
	raw = list(.value = function() METH_OBJ$raw, .private = TRUE, .visible = FALSE),
	cov = list(.value = function() METH_OBJ$cov, .private = TRUE, .visible = FALSE),
	gr = list(.value = function() METH_OBJ$gr, .private = TRUE, .visible = FALSE),
	sample_id = list(.value = function() METH_OBJ$sample_id, .class = "character", .private = TRUE, .visible = FALSE),
	set_chr = list(.value = NULL, .class = "function", .private = TRUE, .visible = FALSE)
)

methylation_hooks$set_chr = function(chr, verbose = TRUE) {
	previous_chr = METH_OBJ$chr
	if(!is.null(previous_chr)) {
		if(previous_chr == chr) {
			if(verbose) qqcat("[@{chr}] @{chr} is already set.\n")
	        return(invisible(NULL))
	    }
	}

	obj = methylation_hooks$get_by_chr(chr)

	# validate obj
	if(is.null(obj$meth)) {
		stop("The list which is returned by `methylation_hooks$get_by_chr` should contain `meth`.")
	}
	if(is.null(obj$gr)) {
		stop("The list which is returned by `methylation_hooks$get_by_chr` should contain `gr`.")
	}
	if(is.null(obj$raw)) {
		cat("The list which is returned by `methylation_hooks$get_by_chr` has no `raw`.")
	}
	if(is.null(obj$cov)) {
		cat("The list which is returned by `methylation_hooks$get_by_chr` has no `cov`.")
	}

	sample_id = colnames(obj$meth)
	if(is.null(sample_id)) {
		stop("The methylation matrix which is represented as `meth` must have column names which represents as sample ids.")
	}

	if(inherits(obj$meth, "data.frame")) {
		obj$meth = as.matrix(obj$meth)
	}
	if(!inherits(obj$meth, "matrix")) {
		stop("`meth` should be a matrix.")
	}
	if(!inherits(obj$gr, "GRanges")) {
		stop("`gr` should be a GRanges object.")
	}
	if(!is.null(obj$raw)) {
		if(inherits(obj$raw, "data.frame")) {
			obj$raw = as.matrix(obj$raw)
		}
		if(!inherits(obj$raw, "matrix")) {
			stop("`raw` should be a matrix.")
		}
	}
	if(!is.null(obj$cov)) {
		if(inherits(obj$cov, "data.frame")) {
			obj$cov = as.matrix(obj$cov)
		}
		if(!inherits(obj$cov, "matrix")) {
			stop("`cov` should be a matrix.")
		}
	}

	if(length(obj$gr) != nrow(obj$meth)) {
		stop("Number of rows in `meth` should be the same as the length of `gr`.")
	}
	if(!is.null(obj$raw)) {
		if(length(obj$gr) != nrow(obj$raw)) {
			stop("Number of rows in `raw` should be the same as the length of `gr`.")
		}
		if(!identical(colnames(obj$meth), colnames(obj$raw))) {
			stop("Column names of `raw` should be identical to the column names of `meth`.")
		}
	}
	if(!is.null(obj$cov)) {
		if(length(obj$gr) != nrow(obj$cov)) {
			stop("Number of rows in `cov` should be the same as the length of `gr`.")
		}
		if(!identical(colnames(obj$meth), colnames(obj$cov))) {
			stop("Column names of `cov` should be identical to the column names of `meth`.")
		}
	}

	METH_OBJ = obj
	METH_OBJ$sample_id = sample_id
	METH_OBJ$chr = chr
	METH_OBJ <<- METH_OBJ

	if(verbose) {
		qqcat("Following methylation datasets have been set for @{chr}:\n")
		qqcat("- `methylation_hooks$gr`: a GRanges object which contains positions of CpG sites.\n")
		qqcat("- `methylation_hooks$meth`: methylation matrix\n")
		if(!is.null(obj$raw)) qqcat("- `methylation_hooks$raw`: raw methylation matrix (unsmoothed)\n")
		if(!is.null(obj$cov)) qqcat("- `methylation_hooks$cov`: CpG coverage matrix\n")
		qqcat("There are @{length(obj$gr)} CpG sites, @{length(sample_id)} samples.\n")
	}
	
	return(invisible(NULL))
}

class(methylation_hooks) = c("methylation_hooks", class(methylation_hooks))

# == title
# Print the general information for the methylation hooks
#
# == param
# -x methylation hook function
# -... other arguments
#
# == value
# Hook functions
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
print.methylation_hooks = function(x, ...) {
	if(is.null(methylation_hooks$get_by_chr)) {
		cat("Importing data by setting `methylation_hooks$get_by_chr`. The value is a function with\n")
		cat("only one argument which is the chromosome name. The returned value should be a list which\n")
		cat("contains:\n")
		qqcat("- `gr`: a GRanges object which contains positions of CpG sites.\n")
		qqcat("- `meth`: methylation matrix\n")
		qqcat("- `raw`: raw methylation matrix (unsmoothed), optional.\n")
		qqcat("- `cov`: CpG coverage matrix, optional.\n")
	} else {
		if(is.null(METH_OBJ$chr)) {
			qqcat("`methylation_hooks$get_by_chr` has been set. Use `methylation_hooks$set_chr() to set a chromosome.\n")
		} else {
			qqcat("`methylation_hooks$get_by_chr` has been set. Current chromosome is @{METH_OBJ$chr}\n")
		}
	}
}

# == title
# Read ChIP-Seq dataset
#
# == param
# -... All the arguments do not make sense here, see 'details' section.
# -RESET .
# -READ.ONLY .
# -LOCAL .
#
# == details
# Unlike methylation dataset which is always stored as matrix, ChIP-Seq dataset is stored
# as a list of peak regions that each one corresponds to peaks in a sample. In many cases, 
# there are ChIP-Seq datasets for multiple histone marks that each mark does not include all
# samples sequenced in e.g. whole genome bisulfite sequencing or RNA-Seq, thus, to import
# such type of flexible data format, users need to define following hook functions:
#
# -sample_id This self-defined function returns a list of sample IDs given the name of a histome mark.
# -peak This function should return a `GenomicRanges::GRanges` object which are peaks for a given
#       histome mark in a given sample. The `GenomicRanges::GRanges` object should better have a meta column 
#       which is the intensity of the histome modification signals.
# -chromHMM This hook is optional. If chromatin segmentation by chromHMM is avaialble, this hoook
#           can be defined as a function which accepts one single sample ID as argument and returns
#           a `GenomicRanges::GRanges` object. The `GenomicRanges::GRanges` object should have a meta column named
#           "states" which is the chromatin states inferred by chromHMM.
#
# == value
# Hook functions
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
chipseq_hooks = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}
chipseq_hooks = setGlobalOptions(
	sample_id = list(.value = function(mark) stop("you need to define `sample_id` hook"),
		             .class = "function",
		             .validate = function(f) length(as.list(f)) == 2,
		             .failed_msg = "The function should only have one argument which is the name of the histome mark."),
	peak = list(.value = function(mark, sid) stop("you need to define `peak` hook"),
		        .class = "function",
		        .validate = function(f) length(as.list(f)) == 3,
		        .failed_msg = "The function should only have two argument which are the name of the histome mark and a sample ID."),
	chromHMM = list(.value = function(sid) stop("you need to define `chromHMM` hook"),
		            .class = "function",
		            .validate = function(f) length(as.list(f)) == 2,
		            .failed_msg = "The function should only have one argument which is the sample ID.")
)

# == title
# Get a list of peak regions for a given histome mark
#
# == param
# -mark name of the histome mark
# -sample_id a vector of sample IDs
#
# == details
# It works after `chipseq_hooks` is set.
#
# == value
# A list of `GenomicRanges::GRanges` objects.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
get_peak_list = function(mark, sample_id = chipseq_hooks$sample_id(mark)) {
    peak_list = lapply(sample_id, function(sid) chipseq_hooks$peak(mark, sid))
    names(peak_list) = sample_id
    peak_list
}
