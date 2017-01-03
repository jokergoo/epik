# getChromInfoFromUCSC = function(species) {
# 	dir = tempdir()

# 	op = qq.options(READ.ONLY = FALSE)
# 	qq.options(code.pattern = "@\\{CODE\\}")

# 	filename = qq("@{species}_getChromInfoFromUCSC")
# 	if(file.exists(qq("@{dir}/@{filename}"))) {
# 		df = read.table(qq("@{dir}/@{filename}"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 	} else {
# 		suppressMessages(df <- GenomicFeatures::getChromInfoFromUCSC(species))
# 		write.table(df, file = qq("@{dir}/@{filename}"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# 	}
	
# 	qq.options(op)

# 	return(df)
# }

getChromInfoFromUCSC = memoise(GenomicFeatures::getChromInfoFromUCSC)
transcriptsBy = memoise(GenomicFeatures::transcriptsBy)
genes = memoise(GenomicFeatures::genes)
transcripts = memoise(GenomicFeatures::transcripts)
exons = memoise(GenomicFeatures::exons)
intronsByTranscript = memoise(GenomicFeatures::intronsByTranscript)
fiveUTRsByTranscript = memoise(GenomicFeatures::fiveUTRsByTranscript)
threeUTRsByTranscript = memoise(GenomicFeatures::threeUTRsByTranscript)
