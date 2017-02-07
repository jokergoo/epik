
library(GetoptLong)
library(bsseq)
PROJECT_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/epik_roadmap"

library(GlobalOptions)
source("/home/guz/project/development/epik/R/read_data_hooks.R")

####### hooks for reading methylation dataset ######
methylation_hooks
methylation_hooks$get_by_chr = function(chr) {
    obj = readRDS(qq("@{PROJECT_DIR}/rds_methylation/@{chr}_roadmap_merged_bsseq.rds"))
    obj2 = list(gr = obj@rowData,
                raw = obj@assays$data$M/obj@assays$data$Cov,
                cov = obj@assays$data$Cov,
                meth = obj@trans(obj@assays$data$coef)
                )
    return(obj2)
}
methylation_hooks
methylation_hooks$meth
methylation_hooks$set_chr("chr21")
methylation_hooks
head(methylation_hooks$meth)
head(methylation_hooks$gr)
methylation_hooks$sample_id

methylation_hooks$set_chr("chr22")


####### hooks for reading chipseq datasets #########

MARKS = c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3", "DNase.macs2")

chipseq_hooks$sample_id = function(mark) {
    sample_id = dir(qq("@{PROJECT_DIR}/data/narrow_peaks"), pattern = qq("E\\d+-@{mark}.narrowPeak.gz"))
    sample_id = gsub(qq("-@{mark}.narrowPeak.gz"), "", sample_id)
	return(sample_id)	
}

chipseq_hooks$peak = function(mark, sid, ...) {
    qqcat("reading peaks: @{sid}, @{mark}\n")
    df = read.table(qq("@{PROJECT_DIR}/data/narrow_peaks/@{sid}-@{mark}.narrowPeak.gz"), stringsAsFactors = FALSE)
    GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]] + 1, df[[3]]), density = df[[5]])
}

chipseq_hooks$chromHMM = function(sid, ...) {
    qqcat("reading chromHMM file for @{sid}\n")
    f = qq("@{PROJECT_DIR}/data/chromatin_states/@{sid}_15_coreMarks_mnemonics.bed.gz")
    gr = read.table(f, sep = "\t", stringsAsFactors = FALSE)
    GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]] + 1, gr[[3]]), states = gr[[4]])
}

for(m in MARKS) {
	qqcat("Sample id for @{m}:\n")
	sid = chipseq_hooks$sample_id(m)
	print(sid)
	print(head(chipseq_hooks$peak(m, sid[1])))
}

chipseq_hooks$chromHMM("E001")
