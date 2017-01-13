
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

