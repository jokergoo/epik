

# == title
# Load and validate configuration file
#
# == param
# -config_file path of configuration file
# -export_env environment where to export variables
# -validate whether do validation
# -use_std_dir whether to create folders
#
# == details
# To run functions in epic package smoothly, users can validate their data by `load_config` function.
# All necessary variables are initialized in a configuration file.
# The configuration file should provide following variables:
#
# ``SAMPLE``: a data frame, row names must be sample id and there must be a
#   'class' column which corresponds to classes of samples. There can also
#   be other annotation columns.
#
# ``COLOR``: a list of color settings corresponding to annotation column in 
#   'SAMPLE'. The value of each list element must be either a named vector
#   or a color mapping function. 'COLOR$class' must be defined or random color
#   will be assigned. Names of other color settings should be same as
#   corresponding columns in 'SAMPLE'.
#
# ``TXDB`` (optional): a `GenomicFeatures::TxDb` object.
#
# ``GTF_FILE`` (optional): GTF file which is used to built ``TXDB``. If it is not specified, the function tries to extract from ``TXDB``.
#
# ``EXPR`` (optional): a matrix which contains expression values. Column names 
#   should be sample id and row names should be gene ids. Note gene ids in the 
#   expression matrix should be same type as genes in `GenomicFeatures::TxDb`.
#
# ``CHROMOSOME``: a vector of chromosome names.
#
# ``GENOME``: abbreviation of GENOME.
#
# ``PROJECT_DIR``: path of output directory. Several sub directories will be created.
#
# ``GENOMIC_FEATURE_LIST``: a list of genomic features as GRanges objects. There
#   must be a element named 'cgi'.
#
# ``MARKS`` (optional): a vector of ChIP-Seq markers.
#
# ``methylation_hooks()`` must be defined.
#
# ``chipseq_hooks()`` is optional unless you want to do integrative analysis.
#
# ``CGI_SHORE_EXTEND``: extension of cgi, by default it is 2kb both upstream and downstream.
#
# == value
# No value is returned.
# 
# == seealso
# `epik`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
load_config = function(config_file, export_env = parent.frame(), validate = TRUE, use_std_dir = FALSE) {
	
	SAMPLE = NULL
	COLOR = NULL
	TXDB = NULL
	EXPR = NULL
	CHROMOSOME = NULL
	GENOME = NULL
	PROJECT_DIR = NULL
	GENOMIC_FEATURE_LIST = NULL
	MARKS = NULL
	GTF_FILE = NULL
	CGI_SHORE_EXTEND = 2000

	message("sourcing ", config_file)
	sys.source(config_file, envir = environment())
	
	if(is.null(GENOMIC_FEATURE_LIST)) {
		GENOMIC_FEATURE_LIST = list()
	}
	if(is.null(CGI_SHORE_EXTEND)) {
		CGI_SHORE_EXTEND = 2000
	}

	if(is.null(SAMPLE)) {
		stop("`SAMPLE` should be defined in the configuration file.")
	}
	if(is.null(COLOR)) {
		stop("`COLOR` should be defined in the configuration file.")
	}
	if(is.null(CHROMOSOME)) {
		stop("`CHROMOSOME` should be defined in the configuration file.")
	}
	if(is.null(GENOME)) {
		stop("`GENOME` should be defined in the configuration file.")
	}
	if(is.null(PROJECT_DIR) && use_std_dir) {
		stop("`PROJECT_DIR` should be defined in the configuration file.")
	}

	if(!validate) {
		assign("SAMPLE", SAMPLE, envir = export_env)
		assign("COLOR", COLOR, envir = export_env)
		assign("TXDB", TXDB, envir = export_env)
		assign("EXPR", EXPR, envir = export_env)
		assign("CHROMOSOME", CHROMOSOME, envir = export_env)
		assign("GENOME", GENOME, envir = export_env)
		assign("PROJECT_DIR", OUTPUT_DIR, envir = export_env)
		assign("GENOMIC_FEATURE_LIST", GENOMIC_FEATURE_LIST, envir = export_env)
		assign("MARKS", MARKS, envir = export_env)
		assign("GTF_FILE", GTF_FILE, envir = export_env)
		assign("CGI_SHORE_EXTEND", CGI_SHORE_EXTEND, envir = export_env)

		if(exists("CHROMATIN_STATES_COLOR")) {
			assign("CHROMATIN_STATES_COLOR", CHROMATIN_STATES_COLOR, envir = export_env)
		}
		if(exists("CHROMATIN_STATES_ORDER")) {
			assign("CHROMATIN_STATES_ORDER", CHROMATIN_STATES_ORDER, envir = export_env)
		}
	}

	# test SAMPLE
	if(is.null(SAMPLE$subgroup)) {
		SAMPLE$subgroup = rep("subgroup_1", nrow(SAMPLE))
		warning("There should be a 'subgroup' column in 'SAMPLE'.")
	}
	
	if(is.null(rownames(SAMPLE))) {
		qq_stop("'SAMPLE must have row names (which are sample IDs and will be used to match other datasets).")
	}

	if(is.null(COLOR$subgroup)) {
		subgroup_level = unique(as.vector(SAMPLE$subgroup))
		COLOR$subgroup = structure(rand_color(length(subgroup_level)), names = subgroup_level)
		warning("`COLOR$subgroup` should be defined.")
	}

	cn = intersect(names(COLOR), colnames(SAMPLE))
	cn_diff = setdiff(colnames(SAMPLE), cn)
	if(length(cn_diff)) {
		qq_message("Following columns in `SAMPLE` are dropped because there is no corresponding color defined in `COLOR`: @{paste(cn_diff, collapse = ', ')}")
	}
	COLOR = COLOR[cn]
	SAMPLE = SAMPLE[, cn, drop = FALSE]

	qq_message("There are @{nrow(SAMPLE)} samples defined in `SAMPLE`.")
	qq_message("@{length(unique(SAMPLE$subgroup))} subgroups: @{paste0(unique(SAMPLE$subgroup))}")
	qq_message("Following sample annotations will be used: @{paste0(cn, collapse = ',')}")

	sample_id = rownames(SAMPLE)

	# initialize the folder structure
	if(use_std_dir) {
		initialize_project_directory(PROJECT_DIR)
	}

	# check chromosome and GENOME
	chrom_info = getChromInfoFromUCSC(GENOME)
	chr_cn = intersect(CHROMOSOME, chrom_info[, 1])
	if(length(chr_cn) == 0) {
		stop("Cannot match 'GENOME' to 'CHROMOSOME'.")
	}
	chrom_info = chrom_info[chrom_info[, 1] %in% chr_cn]
	CHROMOSOME = chr_cn
	qq_message("@{length(chr_cn)} chromosomes after match to @{GENOME}")

	if(!is.null(TXDB)) {
		# sqlite_path = TXDB$conn@dbname

		# TEMP_DIR = paste0(PROJECT_DIR, "/temp")
		# dir.create(TEMP_DIR, showWarnings = FALSE)
		# .__txt_temp_file__. = tempfile(tmpdir = TEMP_DIR, fileext = ".sqlite")
		# file.copy(sqlite_path, .__txt_temp_file__.)

		# TXDB = loadDb(.__txt_temp_file__.)

		# .Last = function() {
		#     file.remove(.__txt_temp_file__.)
		# }
		# assign(".Last", .Last, envir = .GlobalEnv)

		if(is.null(GTF_FILE)) GTF_FILE = metadata(TXDB)[3, "value"]
		if(!file.exists(GTF_FILE)) {
			stop("cannot find ", GTF_FILE)
		}
		if(basename(metadata(TXDB)[3, "value"]) != basename(GTF_FILE)) {
			qq_warning("Base name is not the same for 'GTF_FILE' (@{basename(GTF_FILE)}) and the one used to build 'TXDB' (@{basename(metadata(TXDB)[3, 'value'])})")
		}
	}

	# genomic features
	if(is.character(GENOMIC_FEATURE_LIST)) {
		message("`GENOMIC_FEATURE_LIST` is provided as a character vector, assume they are pathes for bed files.")
		file_exist = file.exists(GENOMIC_FEATURE_LIST)
		if(!all(file_exist)) {
			message("Following bed files do not exist:")
			for(i in which(!file_exist)) {
				qq_message("  @{GENOMIC_FEATURE_LIST[i]}")
			}
			GENOMIC_FEATURE_LIST = GENOMIC_FEATURE_LIST[file_exist]
		}
		if(sum(file_exist) > 0) {
			gn = names(GENOMIC_FEATURE_LIST)
			if(is.null(gn)) {
				message("No 'name' for `GENOMIC_FEATURE_LIST`, use base names as names.")
				gn = basename(GENOMIC_FEATURE_LIST)
			}
			GENOMIC_FEATURE_LIST = lapply(GENOMIC_FEATURE_LIST, function(x) {
				message("reading ", x)
				df = read.table(x, sep = "\t", stringsAsFactors = FALSE)
				GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
			})
			names(GENOMIC_FEATURE_LIST) = gn
		} else {
			GENOMIC_FEATURE_LIST = list()
		}
	}
	if(!is.list(GENOMIC_FEATURE_LIST)) {
		stop("'GENOMIC_FEATURE_LIST' should be a list.")
	}
	if(is.null(GENOMIC_FEATURE_LIST$cgi)) {
		if(exists("CGI")) {
			if(inherits(CGI, c("GRanges", "data.frame"))) {
				qq_message("There is no `GENOMIC_FEATURE_LIST$cgi`, but found `CGI` which is a GRanges object or a data frame, assign it with `CGI`")
				GENOMIC_FEATURE_LIST$cgi = CGI
			}
		}
		stop("'GENOMIC_FEATURE_LIST' must have an element 'cgi'.")
	}
	if(is.null(GENOMIC_FEATURE_LIST$cgi_shore)) {
		if(exists("CGI_SHORE")) {
			if(inherits(CGI_SHORE, c("GRanges", "data.frame"))) {
				qq_message("There is no `GENOMIC_FEATURE_LIST$cgi_shore`, but found `CGI_SHORE` which is a GRanges object or a data frame, assign it with `CGI_SHORE`")
				GENOMIC_FEATURE_LIST$cgi_shore = CGI_SHORE
			}
		} else {
			message("Generate `GENOMIC_FEATURE_LIST$cgi_shore`")
			extended_cgi = GENOMIC_FEATURE_LIST$cgi
			start(extended_cgi) = start(extended_cgi) - CGI_SHORE_EXTEND
			end(extended_cgi) = end(extended_cgi) + CGI_SHORE_EXTEND
			shore = setdiff(extended_cgi, GENOMIC_FEATURE_LIST$cgi)
			GENOMIC_FEATURE_LIST$cgi_shore = shore
		}
	}

	gf_name = names(GENOMIC_FEATURE_LIST)
	for(i in seq_along(GENOMIC_FEATURE_LIST)) {
		if(is.character(GENOMIC_FEATURE_LIST[[i]])) {
			message("reading ", GENOMIC_FEATURE_LIST[[i]])
			GENOMIC_FEATURE_LIST[[i]] = read.table(GENOMIC_FEATURE_LIST[[i]], sep = "\t", stringsAsFactors = FALSE)
		}
		if(is.data.frame(GENOMIC_FEATURE_LIST[[i]])) {
			df = GENOMIC_FEATURE_LIST[[i]]
			GENOMIC_FEATURE_LIST[[i]] = makeGRangesFromDataFrame(df[1:3], 
                    seqnames.field = "V1",
                    start.field = "V2",
                    end.field = "V3")
		}
		GENOMIC_FEATURE_LIST[[i]] = GENOMIC_FEATURE_LIST[[i]][seqnames(GENOMIC_FEATURE_LIST[[i]]) %in% CHROMOSOME]
	}
	
	if(!is.null(TXDB)) {
		message("generate gene level annotations")
		gm = genes(TXDB)
		gm = set_proper_seqlengths(gm, GENOME)
		gm = gm[seqnames(gm) %in% CHROMOSOME]
		if(is.null(GENOMIC_FEATURE_LIST$gene)) {
			message("add `GENOMIC_FEATURE_LIST$gene`")
			GENOMIC_FEATURE_LIST$gene = gm
			strand(GENOMIC_FEATURE_LIST$gene) = "*"
		}
		if(is.null(GENOMIC_FEATURE_LIST$exon)) {
			message("add `GENOMIC_FEATURE_LIST$exon`")
			GENOMIC_FEATURE_LIST$exon = exons(TXDB)
			strand(GENOMIC_FEATURE_LIST$exon) = "*"
		}
		if(is.null(GENOMIC_FEATURE_LIST$intron)) {
			message("add `GENOMIC_FEATURE_LIST$intron`")
			GENOMIC_FEATURE_LIST$intron = setdiff(GENOMIC_FEATURE_LIST$gene, GENOMIC_FEATURE_LIST$exon)
		}
		if(is.null(GENOMIC_FEATURE_LIST$tss)) {
			message("add `GENOMIC_FEATURE_LIST$tss`")
			GENOMIC_FEATURE_LIST$tss = promoters(gm, upstream = 1500, downstream = 500)
			strand(GENOMIC_FEATURE_LIST$tss) = "*"
		}
		if(is.null(GENOMIC_FEATURE_LIST$intergenic)) {
			message("add `GENOMIC_FEATURE_LIST$intergenic`")
			intergenic = gaps(sort(GENOMIC_FEATURE_LIST$gene))
			intergenic = intergenic[ strand(intergenic) == "*"]
			intergenic = intergenic[seqnames(intergenic) %in% CHROMOSOME]
			GENOMIC_FEATURE_LIST$intergenic = intergenic
		}
	}
	qq_message("@{length(GENOMIC_FEATURE_LIST)} genomic features: @{paste0(gf_name, collapse = ', ')}")

	# test methylation_hooks()
	message("pick one chromosome to test methylation_hooks.")
	chr = chrom_info[which.min(chrom_info[, 2])[1], 1]
	methylation_hooks$set_chr(chr, verbose = FALSE)
	meth = methylation_hooks$meth[1:5, ]
	methylation_hooks$gr[1:2]
	meth_sample_id = methylation_hooks$sample_id

	if(length(setdiff(meth_sample_id, sample_id)) > 0) {
		stop("column names in methylation data should be all included in 'SAMPLE'.")
	}
	qq_message("There are @{length(meth_sample_id} samples have methylation data.")

	if(!is.null(EXPR) && !is.null(TXDB)) {
		EXPR = as.matrix(EXPR)
		# test txdb and expr
		if(length(setdiff(colnames(EXPR), sample_id)) > 0) {
			stop("column names in expression data should be all included in 'SAMPLE'.")
		}
		qq_message("There are @{ncol(EXPR)} samples have expression data.\n")
		if(length(intersect(colnames(meth), colnames(EXPR))) == 0) {
			stop("Columns in 'EXPR' have nothing in common to columns in methylation data.")
		}
		EXPR = EXPR[, intersect(colnames(EXPR), meth_sample_id), drop = FALSE]

		qq_message("There are @{ncol(EXPR)} samples both having methylation and expression data.")
		
		genes = genes(TXDB)
		if(length(intersect(names(genes), rownames(EXPR))) == 0) {
			stop("Cannot match genes in 'EXPR' to 'TXDB'.")
		}

		chr = as.vector(seqnames(genes))
		names(chr) = names(genes)
		EXPR = EXPR[chr[rownames(EXPR)] %in% CHROMOSOME, , drop = FALSE]
		message("take @{nrow(EXPR)} genes into analysis.")
	}
	

	# validate chipseq data input
	if(!is.null(MARKS)) {
		for(mk in MARKS) {
			sample_id = chipseq_hooks$sample_id(mk)
			if(length(sample_id) == 0) {
				stop(qq("No ChIP-Seq sample detected for mark @{mk}."))
			}
			if(length(intersect(rownames(SAMPLE), sample_id)) == 0) {
				stop(qq("No ChIP-Seq sample in 'SAMPLE' for mark @{mark}."))
			}
			if(length(setdiff(sample_id, rownames(SAMPLE)))) {
				qq_message("re-define `chipseq_hooks$sample_id()` for @{mk} to only cover samples in `SAMPLE`.")
				.__old_sample_id_hook__. = chipseq_hooks$sample_id
				chipseq_hooks$sample_id = function(mark) {
					sample_id = .__old_sample_id_hook__.(mark)
					intersect(sample_id, rownames(SAMPLE))
				}
				sample_id = chipseq_hooks$sample_id(mk)
			}

			message(mk, ": ", length(sample_id), "samples.")

			sid = sample(sample_id, 1)
			qqcat("random pick one sample: @{sid}\n")
			peak = chipseq_hooks$peak(mk, sid)
			if(!inherits(peak, "GRanges")) {
				stop("chipseq_hooks$peak() should return a GRanges object.")
			}
			# if(length(setdiff(as.character(seqnames(peak)), CHROMOSOME))) {
			# 	.__old_peak_hook__. = chipseq_hooks$peak
			# 	chipseq_hooks$peak = function(mark, sid) {
			# 		gr = .__old_peak_hook__.(mark, sid)
			# 		gr[seqnames(gr) %in% CHROMOSOME]
			# 	}
			# }
			if(!"density" %in% colnames(mcols(peak))) {
				stop(qq("Signal of peaks must be in the 'density' column. Please modify the 'peak' hook."))
			}
			
		}
	}

	assign("SAMPLE", SAMPLE, envir = export_env)
	assign("COLOR", COLOR, envir = export_env)
	assign("TXDB", TXDB, envir = export_env)
	assign("EXPR", EXPR, envir = export_env)
	assign("CHROMOSOME", CHROMOSOME, envir = export_env)
	assign("GENOME", GENOME, envir = export_env)
	assign("PROJECT_DIR", PROJECT_DIR, envir = export_env)
	assign("GENOMIC_FEATURE_LIST", GENOMIC_FEATURE_LIST, envir = export_env)
	assign("MARKS", MARKS, envir = export_env)
	assign("GTF_FILE", GTF_FILE, envir = export_env)
	assign("CGI_SHORE_EXTEND", CGI_SHORE_EXTEND, envir = export_env)

	if(exists("CHROMATIN_STATES_COLOR")) {
		assign("CHROMATIN_STATES_COLOR", CHROMATIN_STATES_COLOR, envir = export_env)
	}
	if(exists("CHROMATIN_STATES_ORDER")) {
		assign("CHROMATIN_STATES_ORDER", CHROMATIN_STATES_ORDER, envir = export_env)
	}

	gc(verbose = FALSE)

	qq_message("\nValidation passed and following global variables are imported: SAMPLE, COLOR, TXDB, EXPR, CHROMOSOME, GENOME, PROJECT_DIR, GENOMIC_FEATURE_LIST, MARKS, GTF_FILE\n")
}

initialize_project_directory = function(x) {
	dir.create(x, showWarnings = FALSE)
	dir.create(paste0(x, "/image"), showWarnings = FALSE)
	dir.create(paste0(x, "/gviz"), showWarnings = FALSE)
	dir.create(paste0(x, "/rds"), showWarnings = FALSE)
	dir.create(paste0(x, "/rds_cr"), showWarnings = FALSE)
	dir.create(paste0(x, "/temp"), showWarnings = FALSE)
	dir.create(paste0(x, "/enriched_heatmap"), showWarnings = FALSE)

	qq_message("creating following folders:")
	qq_message("  @{x}/image/ for general plots")
	qq_message("  @{x}/gviz/ for Gviz plots")
	qq_message("  @{x}/rds/ for rds files")
	qq_message("  @{x}/rds_cr/ for CR-related rds files")
	qq_message("  @{x}/temp/ for temporary files")
	qq_message("  @{x}/enriched_heatmap/")
}

