
methylation_hooks = setGlobalOptions(
	set = list(.value = NULL, .class = "function"),
	get_data = list(.value = function(chr) stop("you need to define `get_data`"),
	                             .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 2
	                             	}),
	meth  = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `meth` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	raw  = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `raw` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	site         = list(.value = function(obj)  stop("you need to define `site` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 3
	                             	}),
	coverage     = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `coverage` hook"),
		                       .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	GRanges      = list(.value =function(obj)  stop("you need to define `GRanges` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 2
	                             	}),
	obj = NULL,
	sample_id = list(.value = NULL, .class = "character")
)

methylation_hooks$set = function(chr) {

    if(!is.null(methylation_hooks$obj)) {
        if(attr(methylation_hooks$obj, "chr") == chr) {
            qqcat("[@{chr}] @{chr} is already set.\n")
            return(invisible(NULL))
        }
    }
    
    obj = methylation_hooks$get_data(chr)
    attr(obj, "chr") = chr

    methylation_hooks$obj = obj

    methylation_hooks$sample_id = colnames(methylation_hooks$meth())

    return(invisible(NULL))
}

chipseq_hooks = setGlobalOptions(
	sample_id = list(.value = function(mark) stop("you need to define `sample_id` hook"),
		             .class = "function",
		             .validate = function(f) length(as.list(f)) == 2),
	peak = list(.value = function(mark, sid) stop("you need to define `peak` hook"),
		        .class = "function",
		        .validate = function(f) length(as.list(f)) == 3),
	chromHMM = list(.value = function(sid) stop("you need to define `chromHMM` hook"),
		            .class = "function",
		            .validate = function(f) length(as.list(f)) == 2)
)

# peak list only for one mark
get_peak_list = function(mark, sample_id = chipseq_hooks$sample_id(mark)) {
    peak_list = lapply(sample_id, function(sid) chipseq_hooks$peak(mark, sid))
    names(peak_list) = sample_id
    peak_list
}
