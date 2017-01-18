library(methods)
library(GetoptLong)

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

########################################
cr = GRanges()
for(chr in CHROMOSOME) {
	qqcat("loading @{chr}...\n")
	cr_tmp = readRDS_or_readRData(qq(paste0(OUTPUT_DIR, "/rds/@{chr}_smoothed_cr_w6s3.rds")))
	l = !is.na(cr_tmp$corr)
	cr = c(cr, cr_tmp[l])
}

cr = copy_cr_attribute(cr_tmp, cr)

# separate neg and pos
neg_cr = cr[cr$corr < 0]
pos_cr = cr[cr$corr > 0]

# add corr_fdr column
neg_cr$corr_fdr = p.adjust(neg_cr$corr_p, "BH")
pos_cr$corr_fdr = p.adjust(pos_cr$corr_p, "BH")

if(!is.null(neg_cr$meth_anova)) {
	neg_cr$meth_anova_fdr = p.adjust(neg_cr$meth_anova, "BH")
	pos_cr$meth_anova_fdr = p.adjust(pos_cr$meth_anova, "BH")
}

saveRDS(neg_cr, file = qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3.rds")))
saveRDS(pos_cr, file = qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3.rds")))


# cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/02.correlated_regions_concatenate.R")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=30:00:00,mem=30G -N correlated_regions_concatenate' '@{cmd}'")
# system(cmd)

