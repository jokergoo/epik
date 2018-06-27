

# == title
# Read methylation dataset
#
# == param
# -... please ignore, see 'details' section.
# -RESET remove all hooks
# -READ.ONLY please ignore
# -LOCAL please ignore
#
# == detail
# Methylation dataset from whole genome bisulfite sequencing is always huge and it does not
# make sense to read them all into the memory. Normally, the methylation dataset is stored
# by chromosome and this hook function can be set to read methylation data in a per-chromosome
# manner. In the package, there are many functions use it internally to read methylation datasets.
#
# Generally, for methylation dataset, there are methylation rate (ranging from 0 to 1), CpG coverage and genomic positions
# for CpG sites. Sometimes there is also smoothed methylation rate. All these datasets can be set
# by defining a proper ``methylation_hooks$get_by_chr``. The value for ``methylation_hooks$get_by_chr``
# is a function with only one argument which is the chromosome name. This function defines how to
# read methylation dataset for a single chromosome. The function must return a list which contains
# following mandatory elements:
#
# -gr a `GenomicRanges::GRanges` object which contains genomic positions for CpG sites. Positions should be sorted.
# -meth a matrix which contains methylation rate. This will be the main methylation dataset the epik
#       package uses, so it should be smoothed methylation rate if the CpG coverage is not high.
#       Note, this matrix must have column names which is sample names and will be used to match
#       other datasets (e.g. RNASeq)
# -cov a matrix which contains CpG coverage.
#
# It can also contain some optional elements and they are not needed for the core analysis:
#
# -raw a matrix which contains unsmoothed methylation rate (or the original methylation rate calculatd
#      as the fraction of methylated CpG in a CpG site)
#
# Note each row in above datasets should correspond to the same CpG site. 
#
# In following example code, assume the methylation data has been processed by bsseq package and saved as
# ``path/bsseq_$chr.rds``, then the definition of ``methylation_hooks$get_by_chr`` is:
#
#   methylation_hooks$get_by_chr = function(chr) {
#       obj = readRDS(paste0("path/bsseq_", chr, ".rds"))
#       lt = list(gr   = granges(obj),
#                 raw  = getMeth(obj, type = "raw"),
#                 cov  = getCoverage(obj, type = "Cov"),
#                 meth = getMeth(obj, type = "smooth")
#       return(lt)
#   }
#
# After ``methylation_hooks$get_by_chr`` is properly set, the "current chromosome" for the methylation dataset
# can be set by ``methylation_hooks$set_chr(chr)`` where ``chr`` is the chromosome name you want to go.
# After validating the dataset, following variables can be used directly:
#
# - ``methylation_hooks$gr``
# - ``methylation_hooks$meth``
# - ``methylation_hooks$sample_id``
# - ``methylation_hooks$cov``
# - ``methylation_hooks$raw`` if available
#
# ``methylation_hooks$set_chr(chr)`` tries to reload the data only when the current chromosome changes.
#
# == value
# Hook functions
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
methylation_hooks = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}

METH_OBJ = new.env()
assign("meth", NULL, envir = METH_OBJ)
assign("raw", NULL, envir = METH_OBJ)
assign("cov", NULL, envir = METH_OBJ)
assign("gr", NULL, envir = METH_OBJ)
assign("chr", NULL, envir = METH_OBJ)
assign("sample_id", NULL, envir = METH_OBJ)

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
			if(verbose) message(qq("[@{chr}] @{chr} is already set."))
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
	if(is.unsorted(start(obj$gr))) {
		stop("`gr` should be sorted.")
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
	
	assign("meth", obj$meth, envir = METH_OBJ)
	assign("raw", obj$raw, envir = METH_OBJ)
	assign("cov", obj$cov, envir = METH_OBJ)
	assign("gr", obj$gr, envir = METH_OBJ)
	assign("chr", chr, envir = METH_OBJ)
	assign("sample_id", sample_id, envir = METH_OBJ)

	if(verbose) {
		message(qq("Following methylation datasets have been set for @{chr}:"))
		message(qq("- `methylation_hooks$gr`: a GRanges object which contains positions of CpG sites."))
		message(qq("- `methylation_hooks$meth`: methylation matrix"))
		if(!is.null(obj$raw)) message(qq("- `methylation_hooks$raw`: raw methylation matrix (unsmoothed)"))
		if(!is.null(obj$cov)) message(qq("- `methylation_hooks$cov`: CpG coverage matrix"))
		message(qq("There are @{length(obj$gr)} CpG sites, @{length(sample_id)} samples."))
	}
	
	return(invisible(NULL))
}

class(methylation_hooks) = c("methylation_hooks", "GlobalOptionsFun")

# == title
# Print the methylation_hooks object
#
# == param
# -x a `methylation_hooks` objects
# -... additional arguments
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
print.methylation_hooks = function(x, ...) {
	if(is.null(methylation_hooks$get_by_chr)) {
		str = "Please set `methylation_hooks$get_by_chr` to import methylation dataset. The value is a function with only one argument which is the chromosome name. The returned value should be a list which contains:\n"
		qqcat(str, strwrap = TRUE)
		qqcat("- `gr`: a GRanges object which contains positions of CpG sites.\n")
		qqcat("- `meth`: methylation matrix (mainly used in the package)\n")
		qqcat("- `raw`: raw methylation matrix (unsmoothed), optional.\n")
		qqcat("- `cov`: CpG coverage matrix.\n")
	} else {
		if(is.null(METH_OBJ$chr)) {
			qqcat("`methylation_hooks$get_by_chr` has been set. Use `methylation_hooks$set_chr() to set a chromosome.\n")
		} else {
			qqcat("`methylation_hooks$get_by_chr` has been set. Current chromosome is @{METH_OBJ$chr}\n")
		}
	}
}

strwrap2 = function(x) {
	paste0(strwrap(x), collapse = "\n")
}


# == title
# Read ChIP-Seq dataset
#
# == param
# -... please ignore, see 'details' section.
# -RESET remove all hooks
# -READ.ONLY please ignore
# -LOCAL please ignore
#
# == details
# Unlike methylation dataset which is always stored as matrix, ChIP-Seq dataset is stored
# as a list of peak regions that each one corresponds to peaks in one sample. In many cases, 
# there are ChIP-Seq datasets for multiple histone marks that each mark does not include all
# samples sequenced in e.g. whole genome bisulfite sequencing or RNA-Seq, thus, to import
# such type of flexible data format, users need to define following hook functions:
#
# -sample_id This self-defined function returns a list of sample IDs given the name of a histone mark.
# -peak This function should return a `GenomicRanges::GRanges` object which are peaks for a given
#       histone mark in a given sample. The `GenomicRanges::GRanges` object should better have a meta column named "density"
#       which is the density of the histone modification signals. (**Note when you want to take the histone
#       modification signals as quatitative analysis, please make sure they are properly normalized between samples**)
# -chromHMM This hook is optional. If chromatin segmentation by chromHMM is avaialble, this hook
#           can be defined as a function which accepts sample ID as argument and returns
#           a `GenomicRanges::GRanges` object. The `GenomicRanges::GRanges` object should have a meta column named
#           "states" which is the chromatin states inferred by chromHMM.
#
# The ``chipseq_hooks$peak()`` must have two arguments ``mark`` and ``sid`` which are the name of the histone mark
# and the sample id. There can also be more arguments such as chromosomes.
#
# As an example, let's assume the peak files are stored in a format of ``path/$sample_id/$mark.bed``, then we can define
# hooks functions as:
#
#   # here `qq` is from GetoptLong package which allows simple variable interpolation
#   chipseq_hooks$sample_id = function(mark) {
#       peak_files = scan(pipe(qq("ls path/*/@{mark}.bed")), what = "character")
#       sample_id = gsub("^path/(.*?)/.*$", "\\1", peak_files)
#       return(sample_id)
#   }
#
#   # here ... is important that the epik package will pass more arguments to it
#   chipseq_hooks$peak = function(mark, sid, ...) {
#       peak_file = qq("path/@{sid}/@{mark}.bed")
#       df = read.table(peak_file, sep = "\t", stringsAsFactors = FALSE)
#       GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), density = df[[5]])
#   }
#
# Normally ``chipseq_hooks$peak()`` are not directly used, it is usually used by `get_peak_list` to read peaks from all samples as a list.
# You can also add more arguments when defining ``chipseq_hooks$peak()`` that these arguments can be passed from `get_peak_list` as well. 
# For example, you can add chromosome name as the third argument that you do not need to read the full dataset at a time:
#
#   # to make it simple, let's assume it only allows one single chromosome
#   chipseq_hooks$peak = function(mark, sid, chr) {
#       peak_file = qq("path/@{sid}/@{mark}.bed")
#       df = read.table(pipe(qq("awk '$1==\"@{chr}\"' @{peak_file}")), sep = "\t", stringsAsFactors = FALSE)
#       GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), density = df[[5]])
#   }
#
# then you can call `get_peak_list` as:
#
#   get_peak_list(mark, chr = "chr1")
#
# The ``chipseq_hooks$chromHMM()`` must have one argument ``sid`` which is the sample id, also there can be more arguments such as chromosomes.
# The usage for the additional argumetns are same as ``chipseq_hooks$peak()``.
#
# == value
# Hook functions
#
# == seealso
# `get_peak_list`, `get_chromHMM_list`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
chipseq_hooks = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}
chipseq_hooks = setGlobalOptions(
	sample_id = list(.value = function(mark) stop("you need to define `sample_id` hook"),
		             .class = "function",
		             .validate = function(f) length(as.list(f)) == 2,
		             .failed_msg = strwrap2("The function should only have one argument which is the name of the histone mark.")),
	peak = list(.value = function(mark, sid, ...) stop("you need to define `peak` hook"),
		        .class = "function",
		        .validate = function(f) length(as.list(f)) >= 4,
		        .failed_msg = strwrap2("The function should have more than two arguments which are the name of the histone mark, sample id and other stuff. If you only use the first two, simply add `...` as the third argument.")),
	chromHMM = list(.value = NULL,
		            .class = "function",
		            .validate = function(f) length(as.list(f)) >= 3,
		            .failed_msg = strwrap2("The function should have more than one arguments which are the sample id and other stuff. If you only use the first one, simply add `...` as the second argument."))
)


# == title
# Get a list of peak regions for a given histone mark
#
# == param
# -mark name of the histone mark
# -sample_id a vector of sample IDs. If not defined, it is the total samples that are available for this histone mark.
# -... more arguments pass to `chipseq_hooks`$peak().
#
# == details
# It works after `chipseq_hooks` is set.
#
# == value
# A list of `GenomicRanges::GRanges` objects.
#
# If you e.g. set "chr" as the third argument when defining `chipseq_hooks`$peak(), "chr" can also be passed here through ``...``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
get_peak_list = function(mark, sample_id = chipseq_hooks$sample_id(mark), ...) {
    peak_list = lapply(sample_id, function(sid) {
    	oe = try(gr <- chipseq_hooks$peak(mark, sid, ...))
    	if(inherits(oe, "try-error")) {
			return(NULL)
		} else {
			return(gr)
		}
    })
    names(peak_list) = sample_id
    peak_list[!sapply(peak_list, is.null)]
}

# == title
# Get a list of chromatin segmentation regions 
#
# == param
# -sample_id a vector of sample IDs.
# -merge if the sample IDs specified are from a same subgroup and user wants to merge them as consensus states
# -window window 
# -min minimal overlap
# -... more arguments pass to `chipseq_hooks`$chromHMM().
#
# == details
# It works after `chipseq_hooks` is set.
#
# == value
# A list of `GenomicRanges::GRanges` objects.
#
# If you e.g. set "chr" as the third argument when defining `chipseq_hooks`$peak(), "chr" can also be passed here through ``...``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
get_chromHMM_list = function(sample_id, merge = FALSE, window = NULL, min = 0.5, ...) {
	lt = lapply(sample_id, function(sid) {
		oe = try(gr <- chipseq_hooks$chromHMM(sid, ...))
		if(inherits(oe, "try-error")) {
			return(NULL)
		} else {
			return(gr)
		}
	})
	names(lt) = sample_id
	lt = lt[!sapply(lt, is.null)]

	if(merge) {
		gr_list_1 = lt
		if(is.null(window)) {
			window = numbers::mGCD(c(sapply(gr_list_1, function(gr) numbers::mGCD(unique(width(gr))))))
			message(qq("window is set to @{window}"))
			if(window == 1) {
				message("when converting bed files to GRanges objects, be careful with the 0-based and 1-based coordinates.")
			}
		}

		# if(is.null(all_states)) {
			all_states = unique(c(unlist(lapply(gr_list_1, function(gr) {unique(as.character(mcols(gr)[, 1]))}))))
			all_states = sort(all_states)
			message(qq("@{length(all_states)} states in total"))
		# }

		message("extracting states")
		m1 = as.data.frame(lapply(gr_list_1, function(gr) {
			k = round(width(gr)/window)
			s = as.character(mcols(gr)[, 1])
			as.integer(factor(rep(s, times = k), levels = all_states))
		}))

		m1 = as.matrix(m1)

		gr = makeWindows(gr_list_1[[1]], w = window)

		message("counting for each state")
		t1 = rowTabulates(m1)
		l = rowMaxs(t1) >= floor(min*length(lt))
		t1 = t1[l, , drop = FALSE]
		gr = gr[l]

		message("determine states")
		states1 = rowWhichMax(t1)

		mcols(gr) = NULL
		gr$states = all_states[states1]
		return(gr)
	} else {
		return(lt)
	}
}
