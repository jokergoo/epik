#######################################################
#           mandatory variables
#######################################################

# SAMPLE should be a data frame, row names must be sample id and there must be a
#   'class' column which corresponds to classes of samples. There can also
#   be other annotation columns.
SAMPLE = data.frame(class = ..., ...)
rownames(SAMPLE) = ...

# COLOR: a list of color settings corresponding to annotation column in 
#   'SAMPLE'. The value of each list element must be either a named vector
#   or a color mapping function. 'COLOR$class' must be defined or random color
#   will be assigned. Names of other color settings should be same as
#   corresponding columns in 'SAMPLE'.
COLOR = c(class = c(...), ...)

# CHROMOSOME: a vector of chromosome names.
CHROMOSOME = paste0("chr", 1:22)

# GENOME: abbreviation of species.
GENOME = "hg19"

# OUTPUT_DIR: path of output directory. Several sub directories will be created.
OUTPUT_DIR = ...

# GENOMIC_FEATURE_LIST: a list of genomic features as GRanges objects. There
#   must be a element named 'cgi' and 'cgi_shore' will be added accordingly. 
#   If `TXDB` is provided, gene, exon, intron, tss, intergenic
#   will be added automatically.
GENOMIC_FEATURE_LIST = list(
    cgi = ...
)

# how to get the object which stores data by chromosome
methylation_hooks$get_data = function(chr) {
	obj = ...
	return(obj)
})

methylation_hooks$meth = function(obj = methylation_hooks$obj, row_index = NULL, col_index = NULL) {
	if(is.null(row_index) && is.null(col_index)) {

    } else if(is.null(row_index)) {

    } else if(is.null(col_index)) {

    } else {

    }
})

methylation_hooks$raw = function(obj = methylation_hooks$obj, row_index = NULL, col_index = NULL) {
	if(is.null(row_index) && is.null(col_index)) {

    } else if(is.null(row_index)) {

    } else if(is.null(col_index)) {

    } else {

    }
})

methylation_hooks$site = function(obj = methylation_hooks$obj, index = NULL) {

})

methylation_hooks$GRanges = function(obj = methylation_hooks$obj) {

})

methylation_hooks$coverage = function(obj = methylation_hooks$obj, row_index = NULL, col_index = NULL) {
	if(is.null(row_index) && is.null(col_index)) {

    } else if(is.null(row_index)) {

    } else if(is.null(col_index)) {

    } else {

    }
})

#######################################################
#           optional variables
#######################################################


# TXDB (optional): a `GenomicFeatures::TxDb` object.
TXDB = loadDb(...)

# GTF file which is used to build 'TXDB'. If it is null, `metadata(TXDB)[3, "value"]` will be used
GTF_FILE = ...

# gene_type in GTF files, use `extract_field_from_gencode(GTF_FILE, level = "gene", primary_key = "gene_id", field = "gene_type")`
# to get all available gene types
GENE_TYPE = c("protein_coding", ...)

# EXPR (optional): a matrix which contains expression values. Column names 
#   should be sample id and row names should be gene ids. Note gene ids in the 
#   expression matrix should be same type as genes in `GenomicFeatures::TxDb`.
EXPR = ...

# MARKS (optional): a vector of ChIP-Seq markers.
MAKRS = c(...)

# chipseq_hooks() is optional unless you want to do integrative analysis.
chipseq_hooks$sample_id = function(mark) {

})

chipseq_hooks$peak = function(mark, sid) {

})

chipseq_hooks$chromHMM = function(sid) {

})

CR_CUTOFF = 0.01
CGI_SHORE_EXTEND = 2000
