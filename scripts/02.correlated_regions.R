library(methods)
suppressPackageStartupMessages(library(GetoptLong))
GetoptLong(
	"chr=s", "chromosomes"
)

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

sample_id = rownames(SAMPLE)

if(!chr %in% CHROMOSOME) stop("'chr' should be included in 'CHROMOSOME'.")

# also apply t-test on subgroup1 and subgroup2
cr = correlated_regions(sample_id, EXPR, TXDB, chr = chr, window_size = 6, window_step = 3, 
	factor = SAMPLE$subgroup, col = COLOR$subgroup)
saveRDS(cr, file = qq("@{OUTPUT_DIR}/rds/@{chr}_smoothed_cr_w6s3.rds"))


# for(chr in paste0("chr", 1:22)) {
# 	cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/02.correlated_regions.R --chr @{chr}")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=30:00:00,mem=12G -N roadmap_cr_@{chr}' '@{cmd}'")
# 	system(cmd)
# }
