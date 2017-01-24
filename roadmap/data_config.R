
library(GetoptLong)
library(bsseq)
PROJECT_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/epik_roadmap/"

SUBGROUP = c(rep("group1", 10), rep("group2", 17))
names(SUBGROUP) = c("E003", "E004", "E005", "E006", "E007", "E011", "E012", "E013", "E016",
              "E024", "E050", "E065", "E066", "E071", "E079", "E094", "E095", "E096",
              "E097", "E098", "E100", "E104", "E105", "E106", "E109", "E112", "E113")
SAMPLE_ID = names(SUBGROUP)
SUBGROUP_COLOR = c("group1" = "red", "group2" = "blue")

methylation_hooks$get_by_chr = function(chr) {
    obj = readRDS(qq("@{PROJECT_DIR}/rds_methylation/@{chr}_roadmap_merged_bsseq.rds"))
    obj2 = list(gr = obj@rowData,
                raw = (obj@assays$data$M/obj@assays$data$Cov)[, SAMPLE_ID],
                cov = obj@assays$data$Cov[, SAMPLE_ID],
                meth = obj@trans(obj@assays$data$coef)[, SAMPLE_ID]
                )
    return(obj2)
}

####### hooks for reading chipseq datasets #########

MARKS = c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3", "DNase.macs2")

chipseq_hooks$sample_id = function(mark) {
    sample_id = dir(qq("@{PROJECT_DIR}/data/narrow_peaks"), pattern = qq("E\\d+-@{mark}.narrowPeak.gz"))
    sample_id = gsub(qq("-@{mark}.narrowPeak.gz"), "", sample_id)
    intersect(sample_id, names(SUBGROUP))	
}

chipseq_hooks$peak = function(mark, sid, ...) {
    df = read.table(qq("@{PROJECT_DIR}/data/narrow_peaks/@{sid}-@{mark}.narrowPeak.gz"), stringsAsFactors = FALSE)
    GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]] + 1, df[[3]]), density = df[[5]])
}

chipseq_hooks$chromHMM = function(sid, ...) {
    f = qq("@{PROJECT_DIR}/data/chromatin_states/@{sid}_15_coreMarks_mnemonics.bed.gz")
    gr = read.table(f, sep = "\t", stringsAsFactors = FALSE)
    GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]] + 1, gr[[3]]), states = gr[[4]])
}

TXDB = loadDb(qq("@{PROJECT_DIR}/txdb/gen10_long_protein_coding_gene_adjusted.sqlite"))

GENE = genes(TXDB)
TSS = promoters(GENE, upstream = 1, downstream = 0)
PROMOTER = promoters(GENE, upstream = 1500, downstream = 500)

CGI = read.table(qq("@{PROJECT_DIR}/data/cpgIslandExt.bed"), stringsAsFactors = FALSE)
CGI = GRanges(seqnames = CGI[[1]], ranges = IRanges(CGI[[2]], CGI[[3]]))
CGI_SHORE = setdiff(flank(CGI, 2000, both = TRUE), CGI)

CHROMOSOME = paste0("chr", 1:22)
GENOME = "hg19"

map = structure(names(GENE), names = gsub("\\.\\d+$", "", names(GENE)))

####################################################
## expression data
expression = list()
count = as.matrix(read.table(qq("@{PROJECT_DIR}/data/expression/57epigenomes.N.pc.gz"), row.names = 1, header = TRUE))
rpkm = as.matrix(read.table(qq("@{PROJECT_DIR}/data/expression/57epigenomes.RPKM.pc.gz"), row.names = 1, header = TRUE))
rownames(count) = map[rownames(count)]
rownames(rpkm) = map[rownames(rpkm)]

count = count[, SAMPLE_ID, drop = FALSE]
rpkm = rpkm[, SAMPLE_ID, drop = FALSE]

######################################################
## genes should have raw count > 0 in at least half samples
l = apply(count, 1, function(x) sum(x > 0) > length(x)/2)
expr = rpkm[l, , drop = FALSE]
EXPR = log2(expr + 1)   # log2(rpkm + 1)
EXPR = EXPR[intersect(rownames(EXPR), names(GENE)), , drop = FALSE]
EXPR = EXPR[as.vector(seqnames(GENE[rownames(EXPR)])) %in% CHROMOSOME, , drop = FALSE]

