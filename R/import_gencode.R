

import_gencode_as_txdb = function(gtf, filter = NULL) {

	message(qq("import @{gtf} as a GRanges object."))
	gr = import(gtf)
	l = eval(substitute(filter), envir = gr@elementMetadata@listData)
	if(!is.null(l)) {
		gr = gr[l]
	}
	message("Making txdb")
	txdb = makeTxDbFromGRanges(gr)

	if(FALSE) {
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

match_by_gencode = function(gtf1, gtf2, filter = NULL) {

	message(qq("import @{gtf1} as a GRanges object."))
	gencode1 = import(gtf1)
	message(qq("import @{gtf2} as a GRanges object."))
	gencode2 = import(gtf2)

	g1 = unique(gencode1$gene_id)
	g2 = unique(gencode2$gene_id)
	map1 = structure(names = gsub("\\.\\d+$", "", g1), g1)
	map2 = structure(names = gsub("\\.\\d+$", "", g2), g2)
	
	l = eval(substitute(filter), envir = gencode1@elementMetadata@listData)
	if(!is.null(l)) {
		gencode1 = gencode1[l]
	}

	l = eval(substitute(filter), envir = gencode2@elementMetadata@listData)
	if(!is.null(l)) {
		gencode2 = gencode2[l]
	}
	
	
	message("construct a full gene bases on its transcripts")
	l = gencode1$type == "transcript"
	g1 = gencode1[l]
	g1 = GRanges(seqnames = gsub("\\.\\d+$", "", g1$gene_id), ranges = ranges(g1), strand = strand(g1))
	g1 = reduce(g1, min = 1e8)
	s1 = start(promoters(g1, upstream = 1, downstream = 0))
	names(s1) = as.vector(seqnames(g1))
	w1 = width(g1)
	names(w1) = as.vector(seqnames(g1))

	l = gencode2$type == "transcript"
	g2 = gencode2[l]
	g2 = GRanges(seqnames = gsub("\\.\\d+$", "", g2$gene_id), ranges = ranges(g2), strand = strand(g2))
	g2 = reduce(g2, min = 1e8)
	s2 = start(promoters(g2, upstream = 1, downstream = 0))
	names(s2) = as.vector(seqnames(g2))
	w2 = width(g2)
	names(w2) = as.vector(seqnames(g2))

	cn = intersect(names(s1), names(s2))

	s1 = s1[cn]
	s2 = s2[cn]
	w1 = w1[cn]
	w2 = w2[cn]

	l = s1 == s2 & w1 == w2

	nm = cn[l]
	message(qq("@{length(nm)}/@{length(map1)} genes are kept."))
	g1 = gencode1[gencode1$gene_id %in% map1[nm]]
	txdb = makeTxDbFromGRanges(g1)
	return(txdb)
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

