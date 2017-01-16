

# == title
# Load and validate configuration file
#
# == param
# -config_file path of configuration file
# -export_env environment where to export variables
# -validate whether do validation
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
# ``GENOME``: abbreviation of species.
#
# ``OUTPUT_DIR``: path of output directory. Several sub directories will be created.
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
# ``CR_CUTOFF``: cutoff for correlation significance of cr.
#
# == value
# No value is returned.
# 
# == seealso
# `epic`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
load_config = function(config_file, export_env = parent.frame(), validate = TRUE) {
	
	SAMPLE = NULL
	COLOR = NULL
	TXDB = NULL
	EXPR = NULL
	CHROMOSOME = NULL
	GENOME = NULL
	GENE_TYPE = NULL
	OUTPUT_DIR = NULL
	GENOMIC_FEATURE_LIST = NULL
	MARKS = NULL
	GTF_FILE = NULL
	CGI_SHORE_EXTEND = 2000
	CR_CUTOFF = 0.01

	cat("sourcing", config_file, "\n")
	sys.source(config_file, envir = environment())

	
	if(is.null(CGI_SHORE_EXTEND)) CGI_SHORE_EXTEND = 2000
	if(is.null(CR_CUTOFF)) CR_CUTOFF = 0.01

	if(!validate) {
		assign("SAMPLE", SAMPLE, envir = export_env)
		assign("COLOR", COLOR, envir = export_env)
		assign("TXDB", TXDB, envir = export_env)
		assign("EXPR", EXPR, envir = export_env)
		assign("CHROMOSOME", CHROMOSOME, envir = export_env)
		assign("GENOME", GENOME, envir = export_env)
		assign("OUTPUT_DIR", OUTPUT_DIR, envir = export_env)
		assign("GENOMIC_FEATURE_LIST", GENOMIC_FEATURE_LIST, envir = export_env)
		assign("MARKS", MARKS, envir = export_env)
		assign("GTF_FILE", GTF_FILE, envir = export_env)
		assign("CGI_SHORE_EXTEND", CGI_SHORE_EXTEND, envir = export_env)
		assign("CR_CUTOFF", CR_CUTOFF, envir = export_env)
	}

	# test SAMPLE
	if(is.null(SAMPLE$class)) {
		SAMPLE$class = rep("sample", nrow(SAMPLE))
		warning("There should be a 'class' column in 'SAMPLE'.")
	}
	
	if(is.null(rownames(SAMPLE))) {
		stop("'SAMPLE must have row names.")
	}

	if(is.null(COLOR$class)) {
		class = unique(COLOR$class)
		COLOR$class = structure(rand_color(length(class)), names = class)
		warning("'COLOR$class' should be defined.")
	}

	cn = intersect(names(COLOR), colnames(SAMPLE))
	COLOR = COLOR[cn]
	SAMPLE = SAMPLE[, cn, drop = FALSE]

	if(is.null(SAMPLE$class)) {
		SAMPLE$class = rep("sample", nrow(SAMPLE))
		warning("There should be a 'class' column in 'SAMPLE'.")
	}
	if(is.null(rownames(SAMPLE))) {
		stop("'SAMPLE must have row names.")
	}

	if(is.null(COLOR$class)) {
		class = unique(COLOR$class)
		COLOR$class = structure(rand_color(length(class)), names = class)
		warning("'COLOR$class' should be defined.")
	}

	cn = intersect(names(COLOR), colnames(SAMPLE))
	COLOR = COLOR[cn]
	SAMPLE = SAMPLE[, cn, drop = FALSE]

	cat(nrow(SAMPLE), "samples\n")
	cat(length(unique(SAMPLE$class)), "classes: ", paste0(unique(SAMPLE$class)), "\n")
	cat("Following sample annotations will be used:", paste0(cn, collapse = ","), "\n")
	cat("\n")

	sample_id = rownames(SAMPLE)

	# initialize the folder structure
	dir.create(OUTPUT_DIR, showWarnings = FALSE)
	dir.create(paste0(OUTPUT_DIR, "/gviz"), showWarnings = FALSE)
	dir.create(paste0(OUTPUT_DIR, "/rds"), showWarnings = FALSE)
	dir.create(paste0(OUTPUT_DIR, "/temp"), showWarnings = FALSE)
	dir.create(paste0(OUTPUT_DIR, "/enriched_heatmap"), showWarnings = FALSE)

	cat("creating following folders:\n")
	cat("  ", OUTPUT_DIR, " for general plots", "\n", sep = "")
	cat("  ", OUTPUT_DIR, "/gviz/ for Gviz plots", "\n", sep = "")
	cat("  ", OUTPUT_DIR, "/rds/ for rds files", "\n", sep = "")
	cat("  ", OUTPUT_DIR, "/temp/ for temporary files", "\n", sep = "")
	cat("  ", OUTPUT_DIR, "/enriched_heatmap/ for enriched heatmaps", "\n", sep = "")
	cat("\n")

	# check chromosome and species
	if(length(intersect(CHROMOSOME, getChromInfoFromUCSC(GENOME)[, 1])) == 0) {
		stop("Cannot match 'GENOME' to 'CHROMOSOME'.")
	}

	if(!is.null(TXDB)) {
		sqlite_path = TXDB$conn@dbname

		TEMP_DIR = paste0(OUTPUT_DIR, "/temp")
		dir.create(TEMP_DIR, showWarnings = FALSE)
		.__txt_temp_file__. = tempfile(tmpdir = TEMP_DIR, fileext = ".sqlite")
		file.copy(sqlite_path, .__txt_temp_file__.)

		TXDB = loadDb(.__txt_temp_file__.)

		.Last = function() {
		    file.remove(.__txt_temp_file__.)
		}
		assign(".Last", .Last, envir = .GlobalEnv)

		if(is.null(GTF_FILE)) GTF_FILE = metadata(TXDB)[3, "value"]
		if(!file.exists(GTF_FILE)) {
			stop("cannot find ", GTF_FILE, "\n")
		}
		if(basename(metadata(TXDB)[3, "value"]) != basename(GTF_FILE)) {
			warning(qq("Base name is not the same for 'GTF_FILE' (@{basename(GTF_FILE)}) and the one used to build 'TXDB' (@{basename(metadata(TXDB)[3, 'value'])})"))
		}
	}

	# genomic features
	if(is.character(GENOMIC_FEATURE_LIST)) {
		gn = names(GENOMIC_FEATURE_LIST)
		GENOMIC_FEATURE_LIST = lapply(GENOMIC_FEATURE_LIST, function(x) {
			cat("reading", x, "\n")
			read.table(x, sep = "\t", stringsAsFactors = FALSE)
		})
		names(GENOMIC_FEATURE_LIST) = gn
	}
	if(!is.list(GENOMIC_FEATURE_LIST)) {
		stop("'GENOMIC_FEATURE_LIST' should be a list.")
	}
	if(is.null(GENOMIC_FEATURE_LIST$cgi)) {
		stop("'GENOMIC_FEATURE_LIST' must have an element 'cgi'.")
	}
	gf_name = names(GENOMIC_FEATURE_LIST)
	for(i in seq_along(GENOMIC_FEATURE_LIST)) {
		if(is.character(GENOMIC_FEATURE_LIST[[i]])) {
			cat("reading", GENOMIC_FEATURE_LIST[[i]], "\n")
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
	if(is.null(GENOMIC_FEATURE_LIST$cgi_shore)) {
		extended_cgi = GENOMIC_FEATURE_LIST$cgi
		start(extended_cgi) = start(extended_cgi) - CGI_SHORE_EXTEND
		end(extended_cgi) = end(extended_cgi) + CGI_SHORE_EXTEND
		shore = setdiff(extended_cgi, GENOMIC_FEATURE_LIST$cgi)
		GENOMIC_FEATURE_LIST$cgi_shore = shore
	}
	if(!is.null(TXDB)) {
		qqcat("generate gene level annotations\n")
		gm = genes(TXDB)
		gm = gm[seqnames(gm) %in% CHROMOSOME]
		if(is.null(GENOMIC_FEATURE_LIST$gene)) {
			GENOMIC_FEATURE_LIST$gene = gm
			strand(GENOMIC_FEATURE_LIST$gene) = "*"
		}
		if(is.null(GENOMIC_FEATURE_LIST$exon)) {
			GENOMIC_FEATURE_LIST$exon = exons(TXDB)
			strand(GENOMIC_FEATURE_LIST$exon) = "*"
		}
		if(is.null(GENOMIC_FEATURE_LIST$intron)) {
			GENOMIC_FEATURE_LIST$intron = setdiff(GENOMIC_FEATURE_LIST$gene, GENOMIC_FEATURE_LIST$exon)
		}
		if(is.null(GENOMIC_FEATURE_LIST$tss)) {
			GENOMIC_FEATURE_LIST$tss = promoters(gm, upstream = 1500, downstream = 500)
			strand(GENOMIC_FEATURE_LIST$tss) = "*"
		}
		if(is.null(GENOMIC_FEATURE_LIST$intergenic)) {
			intergenic = gaps(sort(GENOMIC_FEATURE_LIST$gene))
			intergenic = intergenic[ strand(intergenic) == "*"]
			intergenic = intergenic[seqnames(intergenic) %in% CHROMOSOME]
			GENOMIC_FEATURE_LIST$intergenic = intergenic
		}
	}
	cat(length(GENOMIC_FEATURE_LIST), "genomic features:", paste0(names(GENOMIC_FEATURE_LIST), collapse = ", "), "\n")
	cat("\n")

	# test methylation_hooks()
	cat("Randomly pick one chromosome to test methylation_hooks.\n")
	chr = sample(CHROMOSOME, 1)
	methylation_hooks$set(chr)
	meth = methylation_hooks$meth(row_index = 1:5)
	cov = methylation_hooks$coverage(row_index = 1:5)
	methylation_hooks$site()[1:2]
	methylation_hooks$GRanges()[1:2]
	methylation_hooks$sample_id = colnames(meth)

	if(length(intersect(sample_id, colnames(meth))) == 0) {
		stop("Cannot match column names in methylation data to sample ids in 'SAMPLE'.")
	}
	if(length(intersect(sample_id, colnames(cov))) == 0) {
		stop("Cannot match column names in coverage data to sample ids in 'SAMPLE'.")
	}

	cat("There are", length(intersect(sample_id, colnames(meth))), "samples have methylation data.\n")
	cat("\n")

	if(!is.null(EXPR) && !is.null(TXDB)) {
		EXPR = as.matrix(EXPR)
		# test txdb and expr
		if(length(intersect(sample_id, colnames(EXPR))) == 0) {
			stop("Cannot match column names in 'EXPR' to sample ids in 'SAMPLE'.")
		}
		cat("There are", length(intersect(sample_id, colnames(EXPR))), "samples have expression data.\n")
		if(length(intersect(colnames(meth), colnames(EXPR))) == 0) {
			stop("Cannot match column names in 'EXPR' to column names in methylation data.")
		}
		EXPR = EXPR[, intersect(colnames(EXPR), sample_id), drop = FALSE]

		cat("There are", length(intersect(colnames(meth), colnames(EXPR))), "samples both have methylation and expression data.\n")
		
		genes = genes(TXDB)
		if(length(intersect(names(genes), rownames(EXPR))) == 0) {
			stop("Cannot match genes in 'EXPR' to 'TXDB'.")
		}

		if(!is.null(GENE_TYPE)) {
			gt = extract_field_from_gencode(GTF_FILE, level = "gene", primary_key = "gene_id", field = "gene_type")
			if(any(!GENE_TYPE %in% unique(gt))) {
				stop(qq("'@{paste(setdiff(GENE_TYPE, unique(gt)), collapse = ',')}' is/are not in gene types in GTF file."))
			}
			gt = gt[gt %in% GENE_TYPE]
			EXPR = EXPR[intersect(rownames(EXPR), names(gt)), , drop = FALSE]
		}

		chr = as.vector(seqnames(genes))
		names(chr) = names(genes)
		EXPR = EXPR[chr[rownames(EXPR)] %in% CHROMOSOME, , drop = FALSE]
		cat(nrow(EXPR), "genes into analysis.\n")
		cat("\n")
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
				.__old_sample_id_hook__. = chipseq_hooks$sample_id
				chipseq_hooks$sample_id = function(mark) {
					sample_id = .__old_sample_id_hook__.(mark)
					intersect(sample_id, rownames(SAMPLE))
				}
				sample_id = chipseq_hooks$sample_id(mk)
			}

			cat(mk, ":", length(sample_id), "samples.\n")

			sid = sample(sample_id, 1)
			qqcat("random pick one sample: @{sid}\n")
			peak = chipseq_hooks$peak(mk, sid)
			if(!inherits(peak, "GRanges")) {
				stop("chipseq_hooks$peak() should return a GRanges object.")
			}
			if(length(setdiff(as.character(seqnames(peak)), CHROMOSOME))) {
				.__old_peak_hook__. = chipseq_hooks$peak
				chipseq_hooks$peak = function(mark, sid) {
					gr = .__old_peak_hook__.(mark, sid)
					gr[seqnames(gr) %in% CHROMOSOME]
				}
			}
			if(!"density" %in% colnames(mcols(peak))) {
				stop(qq("Signal of peaks should be in the 'density' column. Please modify the 'peak' hook."))
			}
			
		}
	}

	assign("SAMPLE", SAMPLE, envir = export_env)
	assign("COLOR", COLOR, envir = export_env)
	assign("TXDB", TXDB, envir = export_env)
	assign("EXPR", EXPR, envir = export_env)
	assign("CHROMOSOME", CHROMOSOME, envir = export_env)
	assign("GENOME", GENOME, envir = export_env)
	assign("OUTPUT_DIR", OUTPUT_DIR, envir = export_env)
	assign("GENOMIC_FEATURE_LIST", GENOMIC_FEATURE_LIST, envir = export_env)
	assign("MARKS", MARKS, envir = export_env)
	assign("GTF_FILE", GTF_FILE, envir = export_env)
	assign("CGI_SHORE_EXTEND", CGI_SHORE_EXTEND, envir = export_env)
	assign("CR_CUTOFF", CR_CUTOFF, envir = export_env)

	gc(verbose = FALSE)

	cat("\nValidation passed and following global variables are imported: SAMPLE, COLOR, TXDB, EXPR, CHROMOSOME, GENOME, OUTPUT_DIR, GENOMIC_FEATURE_LIST, MARKS, GTF_FILE\n")
}
