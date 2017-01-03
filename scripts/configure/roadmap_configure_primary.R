
library(GlobalOptions)
library(GenomicFeatures)
library(GetoptLong)
library(matrixStats)
library(HilbertCurve)
library(EnrichedHeatmap)
library(circlize)
library(Gviz)
library(RColorBrewer)
library(dendextend)

source(qq("@{BASE_DIR}/scripts/lib/hooks.R"))
source(qq("@{BASE_DIR}/scripts/lib/methylation_cr.R"))
source(qq("@{BASE_DIR}/scripts/lib/methylation_cr_hc.R"))
source(qq("@{BASE_DIR}/scripts/lib/methylation_cr_enriched.R"))
source(qq("@{BASE_DIR}/scripts/lib/methylation_cr_gviz.R"))
source(qq("@{BASE_DIR}/scripts/lib/chromatin_states_transitions.R"))
source(qq("@{BASE_DIR}/scripts/lib/methylation_distribution.R"))

library(Rcpp)
sourceCpp(qq("@{BASE_DIR}/scripts/lib/extract_sites.cpp"))
sourceCpp(qq("@{BASE_DIR}/scripts/lib/dist_by_closeness.cpp"))


## bsseq_1.2.0
# library(bsseq, lib.loc = "/ibios/tbi_cluster/13.1/x86_64/R/R-3.1.2/lib64/R/library")

#######################################################
##  how to read and get methylation data

methylation_hooks$get_data = function(chr) {

    qqcat("[@{chr}] loading @{BASE_DIR}/rds_methylation/@{chr}_roadmap_merged_bsseq.rds\n")

    obj = readRDS(qq("@{BASE_DIR}/rds_methylation/@{chr}_roadmap_merged_bsseq.rds"))
    obj2 = list(gr = obj@rowData,
                raw = obj@assays$data$M/obj@assays$data$Cov,
                cov = obj@assays$data$Cov,
                meth = obj@trans(obj@assays$data$coef)
                # meth = getMeth(obj, type = "smooth"), 
                # cov = getCoverage(obj, type = "Cov"),
                # raw = getMeth(obj, type = "raw")
                )
    # l = !colnames(obj@assays$data$M) %in% c("E053", "E058")
    # obj2$raw = obj2$raw[, l, drop = FALSE]
    # obj2$cov = obj2$cov[, l, drop = FALSE]
    # obj2$meth = obj2$meth[, l, drop = FALSE]
    return(obj2)
}

methylation_hooks$meth = function(obj = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

    if(is.null(row_index) && is.null(col_index)) {
        obj$meth[, , drop = FALSE]
    } else if(is.null(row_index)) {
        obj$meth[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        obj$meth[row_index, , drop = FALSE]
    } else {
        obj$meth[row_index, col_index, drop = FALSE]
    }

}

methylation_hooks$raw = function(obj = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

   if(is.null(row_index) && is.null(col_index)) {
        obj$raw[, , drop = FALSE]
    } else if(is.null(row_index)) {
        obj$raw[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        obj$raw[row_index, , drop = FALSE]
    } else {
        obj$raw[row_index, col_index, drop = FALSE]
    }
}

methylation_hooks$site = function(obj = methylation_hooks$obj, index = NULL) {
    if(is.null(index))
        start(obj$gr)
    else start(obj$gr[index])
}

methylation_hooks$GRanges = function(obj = methylation_hooks$obj) {
    obj$gr
}

methylation_hooks$coverage = function(obj = methylation_hooks$obj,
    row_index = NULL, col_index = NULL) {

    if(is.null(row_index) && is.null(col_index)) {
        obj$cov[, , drop = FALSE]
    } else if(is.null(row_index)) {
        obj$cov[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        obj$cov[row_index, , drop = FALSE]
    } else {
        obj$cov[row_index, col_index, drop = FALSE]
    }
}

methylation_hooks$set("chr21")
sample_id = colnames(methylation_hooks$meth(row_index = 1:2))


###################################################
## annotations and colors

df = read.table(textConnection(
"eid    group    color    sample_type    sample_type_col
E017    IMR90    #E41A1C    CellLine    #80b2d5
E002    ESC    #924965    PrimaryCulture    #8cd2c7
E008    ESC    #924965    PrimaryCulture    #8cd2c7
E001    ESC    #924965    PrimaryCulture    #8cd2c7
E015    ESC    #924965    PrimaryCulture    #8cd2c7
E014    ESC    #924965    PrimaryCulture    #8cd2c7
E016    ESC    #924965    PrimaryCulture    #8cd2c7
E003    ESC    #924965    PrimaryCulture    #8cd2c7
E024    ESC    #924965    PrimaryCulture    #8cd2c7
E020    iPSC    #69608A    PrimaryCulture    #8cd2c7
E019    iPSC    #69608A    PrimaryCulture    #8cd2c7
E018    iPSC    #69608A    PrimaryCulture    #8cd2c7
E021    iPSC    #69608A    PrimaryCulture    #8cd2c7
E022    iPSC    #69608A    PrimaryCulture    #8cd2c7
E007    ES-deriv    #4178AE    ESCDerived    #bfb9db
E009    ES-deriv    #4178AE    ESCDerived    #bfb9db
E010    ES-deriv    #4178AE    ESCDerived    #bfb9db
E013    ES-deriv    #4178AE    ESCDerived    #bfb9db
E012    ES-deriv    #4178AE    ESCDerived    #bfb9db
E011    ES-deriv    #4178AE    ESCDerived    #bfb9db
E004    ES-deriv    #4178AE    ESCDerived    #bfb9db
E005    ES-deriv    #4178AE    ESCDerived    #bfb9db
E006    ES-deriv    #4178AE    ESCDerived    #bfb9db
E062    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E034    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E045    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E033    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E044    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E043    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E039    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E041    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E042    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E040    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E037    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E048    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E038    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E047    Blood_&_T-cell    #55A354    PrimaryCell    #faf7b4
E029    HSC_&_B-cell    #678C69    PrimaryCell    #faf7b4
E031    HSC_&_B-cell    #678C69    PrimaryCell    #faf7b4
E035    HSC_&_B-cell    #678C69    PrimaryCell    #faf7b4
E051    HSC_&_B-cell    #678C69    PrimaryCell    #faf7b4
E050    HSC_&_B-cell    #678C69    PrimaryCell    #faf7b4
E036    HSC_&_B-cell    #678C69    PrimaryCell    #faf7b4
E032    HSC_&_B-cell    #678C69    PrimaryCell    #faf7b4
E046    HSC_&_B-cell    #678C69    PrimaryCell    #faf7b4
E030    HSC_&_B-cell    #678C69    PrimaryCell    #faf7b4
E026    Mesench    #B65C73    PrimaryCulture    #8cd2c7
E049    Mesench    #B65C73    PrimaryCulture    #8cd2c7
E025    Mesench    #B65C73    PrimaryCulture    #8cd2c7
E023    Mesench    #B65C73    PrimaryCulture    #8cd2c7
E052    Myosat    #E67326    PrimaryCulture    #8cd2c7
E055    Epithelial    #FF9D0C    PrimaryCulture    #8cd2c7
E056    Epithelial    #FF9D0C    PrimaryCulture    #8cd2c7
E059    Epithelial    #FF9D0C    PrimaryCulture    #8cd2c7
E061    Epithelial    #FF9D0C    PrimaryCulture    #8cd2c7
E057    Epithelial    #FF9D0C    PrimaryCulture    #8cd2c7
E058    Epithelial    #FF9D0C    PrimaryCulture    #8cd2c7
E028    Epithelial    #FF9D0C    PrimaryCulture    #8cd2c7
E027    Epithelial    #FF9D0C    PrimaryCell    #faf7b4
E054    Neurosph    #FFD924    PrimaryCulture    #8cd2c7
E053    Neurosph    #FFD924    PrimaryCulture    #8cd2c7
E112    Thymus    #DAB92E    PrimaryTissue    #f57f73
E093    Thymus    #DAB92E    PrimaryTissue    #f57f73
E071    Brain    #C5912B    PrimaryTissue    #f57f73
E074    Brain    #C5912B    PrimaryTissue    #f57f73
E068    Brain    #C5912B    PrimaryTissue    #f57f73
E069    Brain    #C5912B    PrimaryTissue    #f57f73
E072    Brain    #C5912B    PrimaryTissue    #f57f73
E067    Brain    #C5912B    PrimaryTissue    #f57f73
E073    Brain    #C5912B    PrimaryTissue    #f57f73
E070    Brain    #C5912B    PrimaryTissue    #f57f73
E082    Brain    #C5912B    PrimaryTissue    #f57f73
E081    Brain    #C5912B    PrimaryTissue    #f57f73
E063    Adipose    #AF5B39    PrimaryTissue    #f57f73
E100    Muscle    #C2655D    PrimaryTissue    #f57f73
E108    Muscle    #C2655D    PrimaryTissue    #f57f73
E107    Muscle    #C2655D    PrimaryTissue    #f57f73
E089    Muscle    #C2655D    PrimaryTissue    #f57f73
E090    Muscle    #C2655D    PrimaryTissue    #f57f73
E083    Heart    #D56F80    PrimaryTissue    #f57f73
E104    Heart    #D56F80    PrimaryTissue    #f57f73
E095    Heart    #D56F80    PrimaryTissue    #f57f73
E105    Heart    #D56F80    PrimaryTissue    #f57f73
E065    Heart    #D56F80    PrimaryTissue    #f57f73
E078    Sm._Muscle    #F182BC    PrimaryTissue    #f57f73
E076    Sm._Muscle    #F182BC    PrimaryTissue    #f57f73
E103    Sm._Muscle    #F182BC    PrimaryTissue    #f57f73
E111    Sm._Muscle    #F182BC    PrimaryTissue    #f57f73
E092    Digestive    #C58DAA    PrimaryTissue    #f57f73
E085    Digestive    #C58DAA    PrimaryTissue    #f57f73
E084    Digestive    #C58DAA    PrimaryTissue    #f57f73
E109    Digestive    #C58DAA    PrimaryTissue    #f57f73
E106    Digestive    #C58DAA    PrimaryTissue    #f57f73
E075    Digestive    #C58DAA    PrimaryTissue    #f57f73
E101    Digestive    #C58DAA    PrimaryTissue    #f57f73
E102    Digestive    #C58DAA    PrimaryTissue    #f57f73
E110    Digestive    #C58DAA    PrimaryTissue    #f57f73
E077    Digestive    #C58DAA    PrimaryTissue    #f57f73
E079    Digestive    #C58DAA    PrimaryTissue    #f57f73
E094    Digestive    #C58DAA    PrimaryTissue    #f57f73
E099    Other    #999999    PrimaryTissue    #f57f73
E086    Other    #999999    PrimaryTissue    #f57f73
E088    Other    #999999    PrimaryTissue    #f57f73
E097    Other    #999999    PrimaryTissue    #f57f73
E087    Other    #999999    PrimaryTissue    #f57f73
E080    Other    #999999    PrimaryTissue    #f57f73
E091    Other    #999999    PrimaryTissue    #f57f73
E066    Other    #999999    PrimaryTissue    #f57f73
E098    Other    #999999    PrimaryTissue    #f57f73
E096    Other    #999999    PrimaryTissue    #f57f73
E113    Other    #999999    PrimaryTissue    #f57f73
"), header = TRUE, row.names = 1, stringsAsFactors = FALSE, comment.char = "")

sample_id = intersect(sample_id, rownames(df))
SAMPLE = data.frame(id = sample_id, group = df[sample_id, "group"], sample_type = df[sample_id, "sample_type"], stringsAsFactors = FALSE)
rownames(SAMPLE) = sample_id
df = df[sample_id, ]
COLOR = list(group = structure(unique(df[, "color"]), names = unique(df[, "group"])),
             sample_type = structure(unique(df[, "sample_type_col"]), names = unique(df[, "sample_type"])))


CHROMOSOME = paste0("chr", c(1:22))

GENOME = "hg19"

OUTPUT_DIR = BASE_DIR

##################################################
## transcriptome annotation

cat("load txdb...\n")
TXDB = loadDb(qq("@{BASE_DIR}/data/gen10_long_protein_coding_gene_adjusted.sqlite"))
GTF_FILE = qq("@{BASE_DIR}/data/gen10_long_protein_coding.gtf")
genes = genes(TXDB)
map = structure(names(genes), names = gsub("\\.\\d+$", "", names(genes)))

####################################################
## expression data

cat("load expression...\n")
expression = list()
count = as.matrix(read.table(qq("@{BASE_DIR}/data/expression/57epigenomes.N.pc.gz"), row.names = 1, header = TRUE))
rpkm = as.matrix(read.table(qq("@{BASE_DIR}/data/expression/57epigenomes.RPKM.pc.gz"), row.names = 1, header = TRUE))
rownames(count) = map[rownames(count)]
rownames(rpkm) = map[rownames(rpkm)]
sample_id = intersect(sample_id, colnames(count))


######################################################
## only took samples that have both RNASeq and WGBS
count = count[, sample_id, drop = FALSE]
rpkm = rpkm[, sample_id, drop = FALSE]

SAMPLE = SAMPLE[sample_id, ]
COLOR$group = COLOR$group[unique(SAMPLE$group)]
COLOR$sample_type = COLOR$sample_type[unique(SAMPLE$sample_type)]

######################################################
## genes should have raw count > 0 in at least half samples
l = apply(count, 1, function(x) sum(x > 0) > length(x)/2)
expr = rpkm[l, , drop = FALSE]
EXPR = log2(expr + 1)   # log2(rpkm + 1)
EXPR = EXPR[intersect(rownames(EXPR), names(genes)), , drop = FALSE]

GENE_TYPE = "protein_coding"
df = read.table(pipe(qq("perl @{BASE_DIR}/scripts/perl/extract_field_from_gencode.pl @{GTF_FILE} gene gene_id gene_type")), stringsAsFactors = FALSE)
gt = structure(df[[2]], names = df[[1]])
gt = gt[gt %in% GENE_TYPE]
# only protein coding genes
EXPR = EXPR[intersect(rownames(EXPR), names(gt)), , drop = FALSE]
gene_with_chr = structure(as.vector(seqnames(genes)), names = names(genes))
EXPR = EXPR[gene_with_chr[rownames(EXPR)] %in% CHROMOSOME, ]

###############################################
## available histone marks
MARKS = c("H3K4me1", # enhancer , core
          "H3K4me3",  # promoter, core
          "H3K27ac",  # promoter + enhancer
          "H3K27me3",  # repressed genes, core
          "H3K36me3",  # actively transcribed gene bodies, core
          "H3K9me3",   # repressed genes, core, we removed this mark because the signal is relatively weak which is not sufficient for systematic analysis
          "DNase.macs2"
          )

###############################################
## since not all sample have a certain mark, this hook
## returns samples that have corresponding mark
chipseq_hooks$sample_id = function(mark) {
    sample_id = dir(qq("@{BASE_DIR}/data/narrow_peaks"), pattern = qq("E\\d+-@{mark}.narrowPeak.gz"))
    sample_id = gsub(qq("-@{mark}.narrowPeak.gz"), "", sample_id)
    intersect(sample_id, rownames(SAMPLE))
}


#################################################
## how to get the peak of a certain mark in a sample
chipseq_hooks$peak = function(mark, sid) {
    qqcat("reading peaks: @{sid}, @{mark}\n")
    df = read.table(qq("@{BASE_DIR}/data/narrow_peaks/@{sid}-@{mark}.narrowPeak.gz"), stringsAsFactors = FALSE)
    GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), density = df[[5]])
}

#################################################
## CpG islands and shores (flanking 2kb)
CGI = read.table(qq("@{BASE_DIR}/data/cpgIslandExt.bed"), stringsAsFactors = FALSE)
CGI = GRanges(seqnames = CGI[[1]], ranges = IRanges(CGI[[2]], CGI[[3]]))
CGI_SHORE = setdiff(flank(CGI, 2000, both = TRUE), CGI)
