

import_gencode_as_txdb = function(gtf, filter = NULL, include_mappings = FALSE) {

	qqcat("import @{gtf} as a GRanges object.\n")
	gr = import(gtf)
	l = eval(substitute(filter), envir = gr@elementMetadata@listData)
	if(!is.null(l)) {
		gr = gr[l]
	}
	cat("Making txdb\n")
	txdb = makeTxDbFromGRanges(gr)

	if(include_mappings) {
		cat("Generate mappings\n")
		l1 = !duplicated(gr$gene_id)
		l2 = !duplicated(gr$transcript_id)

		gene_id2name = NULL
		if(!is.null(gr$gene_name)) gene_id2name = structure(gr$gene_name[l1], names = gr$gene_id[l1])
		
		gene_id2type = NULL
		if(!is.null(gr$gene_type)) gene_id2type = structure(gr$gene_type[l1], names = gr$gene_id[l1])
		
		tx_id2name = NULL
		if(!is.null(gr$transcript_name)) tx_id2name = structure(gr$transcript_name[l2], names = gr$transcript_id[l2])
		
		tx_id2type = NULL
		if(!is.null(gr$transcript_type)) tx_id2type = structure(gr$transcript_type[l2], names = gr$transcript_id[l2])
		
		mapping = list(gene_id2name = gene_id2name,
			           gene_id2type = gene_id2type,
			           tx_id2name = tx_id2name,
			           tx_id2type = tx_id2type)

		return(list(txdb = txdb, mapping = mapping))
	} else {
		return(txdb)
	}
}


# == title 
# Extract field from gencode GTF file
#
# == param
# -file the input GTF file
# -level level of the annotation (e.g. gene, transcript, exon, ...)
# -primary_key primary field
# -field field to be retrieved
#
# == details
# Although GTF file can be imported by `GenomicFeatures::makeTranscriptDbFromGFF`, some information
# in the original GTF file will not be imported. This function aims to extract additionally information
# from GTF file.
#
# The function calls external perl script, so you need to have perl installed.
#
# == value
# A vector in which 'primary_key' corresponds to the name and 'field' corresponds to the value.
#
# == seealso
# `available_gencode_fields` lists all possible values for ``primary_key`` and ``field``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# download.file("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
#     destfile = "gencode.v19.annotation.gtf.gz")
# extract_field_from_gencode("gencode.v19.annotation.gtf.gz")
# extract_field_from_gencode("gencode.v19.annotation.gtf.gz", field = "gene_type")
# }
# NULL
extract_field_from_gencode = function(file, level = "gene", 
	primary_key = "gene_id", field = "gene_name") {
	df = read.table(pipe(qq("perl \"@{system.file(package = 'epic')}/perl_scripts/extract_field_from_gencode.pl\" @{file} @{level} @{primary_key} @{field}")), 
		stringsAsFactors = FALSE)
	return(structure(df[[2]], names = df[[1]]))
}

# == title
# Returns all supported fields in GTF data
#
# == param
# -file the input GTF file
# -level level of the annotation (e.g. gene, transcript, exon, ...)
#
# == details
# These fields are stored in the 9th column in the GTF file.
#
# == value
# A vector of available fields.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# \dontrun{
# download.file("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
#     destfile = "gencode.v19.annotation.gtf.gz")
# available_gencode_fields("gencode.v19.annotation.gtf.gz", level = "gene")
# available_gencode_fields("gencode.v19.annotation.gtf.gz", level = "transcript")
# }
# NULL
available_gencode_fields = function(file, level = "gene") {
	con = file(file, "r")
	while(1) {
		line_str = readLines(con, n = 1)
		if(grepl("^#", line_str)) next

		line = strsplit(line_str, "\t")[[1]]
		if(line[3] != level) {
			next
		} else {
			pair = strsplit(line[9], ";")[[1]]
			field = gsub("^\\s*(\\w+)\\s.*$", "\\1", pair)
			close(con)
			return(field)
		}
	}
	close(con)
	return(NULL)
}

