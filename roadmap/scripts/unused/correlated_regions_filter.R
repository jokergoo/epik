
suppressPackageStartupMessages(library(GetoptLong))

cutoff = 0.05
GetoptLong("cutoff=f", "cutoff for filter cr")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

cr_filtered = filter_correlated_regions(chromosome = CHROMOSOME, template = qq("`OUTPUT_DIR`/rds/@{chr}_cr.rds", code.pattern = "`CODE`"), cutoff = cutoff)
saveRDS(cr_filtered, file = qq("@{OUTPUT_DIR}/rds/cr_filtered_fdr_@{cutoff}.rds"))

cr_filtered = filter_correlated_regions(chromosome = CHROMOSOME, template = qq("`OUTPUT_DIR`/rds/@{chr}_cr.rds", code.pattern = "`CODE`"), 
	cutoff = cutoff, type = "neg")
saveRDS(cr_filtered, file = qq("@{OUTPUT_DIR}/rds/cr_filtered_fdr_@{cutoff}_neg.rds"))

cr_filtered = filter_correlated_regions(chromosome = CHROMOSOME, template = qq("`OUTPUT_DIR`/rds/@{chr}_cr.rds", code.pattern = "`CODE`"), 
	cutoff = cutoff, type = "pos")
saveRDS(cr_filtered, file = qq("@{OUTPUT_DIR}/rds/cr_filtered_fdr_@{cutoff}_pos.rds"))
