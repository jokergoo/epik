
library(GetoptLong)
library(bsseq)
BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/"

subgroup = c(rep("group1", 10), rep("group2", 17))
names(subgroup) = c("E003", "E004", "E005", "E006", "E007", "E011", "E012", "E013", "E016",
              "E024", "E050", "E065", "E066", "E071", "E079", "E094", "E095", "E096",
              "E097", "E098", "E100", "E104", "E105", "E106", "E109", "E112", "E113")
sample_id = names(subgroup)

methylation_hooks$get_by_chr = function(chr) {
    obj = readRDS(qq("@{BASE_DIR}/rds_methylation/@{chr}_roadmap_merged_bsseq.rds"))
    sample_id = c("E003", "E004", "E005", "E006", "E007", "E011", "E012", "E013", "E016",
              "E024", "E050", "E065", "E066", "E071", "E079", "E094", "E095", "E096",
              "E097", "E098", "E100", "E104", "E105", "E106", "E109", "E112", "E113")

    obj2 = list(gr = obj@rowData,
                raw = (obj@assays$data$M/obj@assays$data$Cov)[, sample_id],
                cov = obj@assays$data$Cov[, sample_id],
                meth = obj@trans(obj@assays$data$coef)[, sample_id]
                )
    return(obj2)
}

####### hooks for reading chipseq datasets #########

MARKS = c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3", "DNase.macs2")

chipseq_hooks$sample_id = function(mark) {
    sample_id = dir(qq("@{BASE_DIR}/data/narrow_peaks"), pattern = qq("E\\d+-@{mark}.narrowPeak.gz"))
    sample_id = gsub(qq("-@{mark}.narrowPeak.gz"), "", sample_id)
	return(sample_id)	
}

chipseq_hooks$peak = function(mark, sid, ...) {
    qqcat("reading peaks: @{sid}, @{mark}\n")
    df = read.table(qq("@{BASE_DIR}/data/narrow_peaks/@{sid}-@{mark}.narrowPeak.gz"), stringsAsFactors = FALSE)
    GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]] + 1, df[[3]]), density = df[[5]])
}

chipseq_hooks$chromHMM = function(sid) {
    qqcat("reading chromHMM file for @{sid}\n")
    f = qq("@{BASE_DIR}/data/chromatin_states/@{sid}_15_coreMarks_mnemonics.bed.gz")
    gr = read.table(f, sep = "\t", stringsAsFactors = FALSE)
    GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]] + 1, gr[[3]]), states = gr[[4]])
}

CGI = read.table(qq("@{BASE_DIR}/data/cpgIslandExt.bed"), stringsAsFactors = FALSE)
CGI = GRanges(seqnames = CGI[[1]], ranges = IRanges(CGI[[2]], CGI[[3]]))
CGI_SHORE = setdiff(flank(CGI, 2000, both = TRUE), CGI)
