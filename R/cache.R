transcriptsBy = memoise(GenomicFeatures::transcriptsBy)
genes = memoise(GenomicFeatures::genes)
transcripts = memoise(GenomicFeatures::transcripts)
exons = memoise(GenomicFeatures::exons)
intronsByTranscript = memoise(GenomicFeatures::intronsByTranscript)
fiveUTRsByTranscript = memoise(GenomicFeatures::fiveUTRsByTranscript)
threeUTRsByTranscript = memoise(GenomicFeatures::threeUTRsByTranscript)

getChromInfoFromUCSC = function(genome, goldenPath_url = "http://hgdownload.cse.ucsc.edu/goldenPath") {
	oe = try(df <- GenomicFeatures::getChromInfoFromUCSC(genome, goldenPath_url))
	if(inherits(oe, "try-error")) {
		data = readRDS(system.file("extdata", "chrom_info_list.rds", package = "epik"))
		return(data[[genome]])
	} else {
		return(df)
	}
}
if(!is.memoised(getChromInfoFromUCSC)) {
	getChromInfoFromUCSC = memoise(getChromInfoFromUCSC)
}

read.chromInfo = function(chromInfo = paste0(system.file(package = "circlize"),
    "/extdata/chromInfo.txt"), species = NULL, chromosome.index = NULL,
    sort.chr = TRUE) {
	oe = try(res <- circlize::read.chromInfo(chromInfo, species, chromosome.index, sort.chr))
	if(inherits(oe, "try-error")) {
		data = readRDS(system.file("extdata", "chrom_info_list.rds", package = "epik"))
		res = circlize::read.chromInfo(chromInfo = data[[species]])
	} 
	return(res)
}
if(!is.memoised(read.chromInfo)) {
	read.chromInfo = memoise(read.chromInfo)
}
