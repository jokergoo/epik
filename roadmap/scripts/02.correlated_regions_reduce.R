library(methods)
suppressPackageStartupMessages(library(GetoptLong))
cutoff = 0.05
meandiff = 0.2
direction = "neg"
type = "fdr"
region = NULL
NCORE = 1
GetoptLong("cutoff=f", "cutoff for filter cr",
	       "type=s", "corr|fdr",
	       "direction=s", "pos|neg",
	       "meandiff=f", "cutoff for meandiff",
	       "region=s", "tss|gene|intergenic",
	       "NCORE=i", "number of cores")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

neg_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3.rds")))
pos_cr = readRDS(qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3.rds")))


if(type == "fdr") {
	if(direction == "neg") {
		if(!is.null(neg_cr$meth_anova)) {
			neg_cr2 = reduce_cr(neg_cr[neg_cr$corr_fdr < cutoff & neg_cr$meth_anova_fdr < cutoff & neg_cr$meth_diameter > meandiff], TXDB, mc.cores = NCORE)
		} else {
			neg_cr2 = reduce_cr(neg_cr[neg_cr$corr_fdr < cutoff & neg_cr$meth_IQR > meandiff], TXDB, mc.cores = NCORE)
		}
		saveRDS(neg_cr2, file = qq(paste0(OUTPUT_DIR, "/rds/all_neg_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds")))
	} else {
		if(!is.null(pos_cr$meth_anova)) {
			pos_cr2 = reduce_cr(pos_cr[pos_cr$corr_fdr < cutoff & pos_cr$meth_anova_fdr < cutoff & pos_cr$meth_diameter > meandiff], TXDB, mc.cores = NCORE)
		} else {
			pos_cr2 = reduce_cr(pos_cr[pos_cr$corr_fdr < cutoff & pos_cr$meth_IQR > meandiff], TXDB, mc.cores = NCORE)
		}
		saveRDS(pos_cr2, file = qq(paste0(OUTPUT_DIR, "/rds/all_pos_cr_w6s3_fdr_less_than_@{cutoff}_methdiff_larger_than_@{meandiff}.rds")))
	}
} 

# for(cutoff in c(0.1, 0.05, 0.01)) {
#     for(meandiff in c(0, 0.1, 0.2, 0.3)) {
#         cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/correlated_regions_reduce.R --type fdr --cutoff @{cutoff} --direction neg --meandiff @{meandiff} --NCORE 2")
#         cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=30:00:00,mem=10G,nodes=1:ppn=2 -N neg_cr_reduce_fdr_@{cutoff}_meandiff_@{meandiff}' '@{cmd}'")
#         system(cmd)

#         cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/correlated_regions_reduce.R --type fdr --cutoff @{cutoff} --direction pos --meandiff @{meandiff} --NCORE 2")
#         cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=30:00:00,mem=10G,nodes=1:ppn=2 -N pos_cr_reduce_fdr_@{cutoff}_meandiff_@{meandiff}' '@{cmd}'")
#         system(cmd)
#     }
# }
