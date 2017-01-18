library(methods)
suppressPackageStartupMessages(library(GetoptLong))
chr = "chr21"
GetoptLong("chr=s", "a single chromosome, should have 'chr' prefix")

suppressPackageStartupMessages(library(GenomicRanges))

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
OUTPUT_DIR = BASE_DIR

# this xls file is actually a plain text file
meta = read.table(qq("@{BASE_DIR}/data/EG.mnemonics.name.xls"), sep = "\t", stringsAsFactors = FALSE)
cn = meta[[1]]

qqcat("reading methylation rate for @{chr}...\n")
meth = read.table(qq("@{BASE_DIR}/data/FractionalMethylation_Removed_E027_E064_Fixed_E012/@{chr}.fm"))
qqcat("reading CpG coverage for @{chr}...\n")
cov = read.table(qq("@{BASE_DIR}/data/ReadCoverage_Removed_E027_E064/@{chr}.rc"))

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


# check the methylation is on CpG level or C level
# if on C level, calcualte the mean methylation of two Cs weighted by coverage
j = 0
meth2 = matrix(nrow = nrow(meth), ncol = ncol(meth))
cov2 = matrix(nrow = nrow(cov), ncol = ncol(cov))
pos2 = numeric(length(pos))
two_rows_are_same = function(x1, x2) {
	sum(abs(x1 - x2) < 0.05)/length(x1) > 0.8
}

two_row_means = function(x1, x2, w1, w2) {
	y = (x1*w1 + x2*w2)/(w1 + w2)
	y[is.na(y)] = 0
	y
}

next_used = FALSE
for(i in seq_len(nrow(meth))) {
	if(next_used) {
		next_used = FALSE
		next
	}

	if(pos[i] + 1 == pos[i + 1]) {
		# check whether the methylation values are same in all samples
		if(two_rows_are_same(meth[i, ], meth[i+1, ])) {
			j = j + 1
			meth2[j, ] = two_row_means(meth[i, ], meth[i+1, ], cov[i, ], cov[i+1, ])
			cov2[j, ] = round(two_row_means(cov[i, ], cov[i+1, ], cov[i, ], cov[i+1, ]))
			pos2[j] = pos[i]
			next_used = TRUE
		}
	} else { 
		j = j + 1
		meth2[j, ] = meth[i, ]
		cov2[j, ] = round(cov[i, ])
		pos2[j] = pos[i]
		cat("find an alone C", i, "\n")
	}
}
l = apply(meth2, 1, function(x) all(is.na(x)))
meth2 = meth2[!l, , drop = FALSE]
cov2 = cov2[!l, , drop = FALSE]
pos2 = pos2[!l]

colnames(meth2) = colnames(meth)
colnames(cov2) = colnames(cov)

# more than 50% of samples have CpG coverage > 2
cpg_cov_cutoff = min(c(5, quantile(cov2[cov2 > 0], 0.1)))
x = apply(cov2, 1, function(x) sum(x >= cpg_cov_cutoff)/length(x))
l = x >= 0.5

meth2 = meth2[l, ]
cov2 = cov2[l, ]
pos2 = pos2[l]

library(bsseq)
bsseq = BSseq(M = round(meth2*cov2), Cov = cov2, pos = pos2, chr = chr, sampleNames = colnames(meth2))
bsseq = BSmooth(bsseq, parallelBy = "sample", mc.cores = 2, verbose = TRUE)
saveRDS(bsseq, file = qq("@{OUTPUT_DIR}/rds_methylation/@{chr}_roadmap_merged_bsseq.rds"))


# for(chr in paste0("chr", c(1:22, "X", "Y"))) {
#     cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/00.make_rds_for_methylation_data.R --chr @{chr} --output /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/rds_methylation")
#     cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=20G,nodes=1:ppn=2 -N cr_smooth_@{chr}' '@{cmd}'")
#     system(cmd)
# }
