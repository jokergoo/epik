library(methods)
suppressPackageStartupMessages(library(GetoptLong))
chr = "chr21"
GetoptLong("chr=s", "a single chromosome, should have 'chr' prefix")

suppressPackageStartupMessages(library(GenomicRanges))

PROJECT_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/epik_roadmap/"

# this xls file is actually a plain text file
meta = read.table(qq("@{PROJECT_DIR}/data/EG.mnemonics.name.xls"), sep = "\t", stringsAsFactors = FALSE)
cn = meta[[1]]

qqcat("reading methylation rate for @{chr}...\n")
meth = read.table(qq("@{PROJECT_DIR}/data/FractionalMethylation_Removed_E027_E064_Fixed_E012/@{chr}.fm"))
qqcat("reading CpG coverage for @{chr}...\n")
cov = read.table(qq("@{PROJECT_DIR}/data/ReadCoverage_Removed_E027_E064/@{chr}.rc"))

pos = meth[, 1]
meth = as.matrix(meth[, -1])
cov = as.matrix(cov[, -1])
meth[meth < -0.1] = 0
cov[cov < -0.1] = 0
colnames(meth) = cn
colnames(cov) = cn

# NA rate: E054: 14%, E070:%13, E084:18%, E085:41%
l = cn %in% c("E054", "E070", "E084", "E085")
meth = meth[, !l]
cov = cov[, !l]

source("/home/guz/project/development/epik/R/cpg_dinucleotide_methylation.R")
lt = cpg_dinucleotide_methylation(pos, meth, cov)
pos2 = lt$pos
meth2 = lt$meth
cov2 = lt$cov

# more than 50% of samples have CpG coverage > 2
cpg_cov_cutoff = min(c(5, quantile(cov2[cov2 > 0], 0.1)))
x = apply(cov2, 1, function(x) sum(x >= cpg_cov_cutoff)/length(x))
l = x >= 0.5

meth2 = meth2[l, ]
cov2 = cov2[l, ]
pos2 = pos2[l]

library(bsseq)
bsseq = BSseq(M = round(meth2*cov2), Cov = cov2, pos = pos2, chr = chr, sampleNames = colnames(meth2))
bsseq = BSmooth(bsseq, parallelBy = "sample", mc.cores = 1, verbose = TRUE)
saveRDS(bsseq, file = qq("@{PROJECT_DIR}/rds_methylation/@{chr}_roadmap_merged_bsseq.rds"))

