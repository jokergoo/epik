
library(methods)
suppressPackageStartupMessages(library(GetoptLong))

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

sample_id = rownames(SAMPLE)

ha = HeatmapAnnotation(group = SAMPLE$group, 
                       sample_type = SAMPLE$sample_type,
                       subgroup = SAMPLE$subgroup,
                       col = list(group = COLOR$group, 
                                  sample_type = COLOR$sample_type, 
                                  subgroup = COLOR$subgroup))

pdf(qq("@{OUTPUT_DIR}/plots/global_methylation_distribution.pdf"))
global_methylation_distribution(sample_id, annotation = SAMPLE$subgroup, annotation_color = COLOR$subgroup,
        chromosome = CHROMOSOME, ha = ha)
dev.off()


chromInfo = getChromInfoFromUCSC(GENOME)
chromInfo = chromInfo[chromInfo$chrom %in% CHROMOSOME, ]
chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
complement = setdiff(chromGr, union(CGI, CGI_SHORE))

pdf(qq("@{OUTPUT_DIR}/plots/global_methylation_distribution_cgi.pdf"))
global_methylation_distribution(sample_id, annotation = SAMPLE$subgroup, annotation_color = COLOR$subgroup,
        chromosome = CHROMOSOME, ha = ha, background = CGI, p = 0.01)
dev.off()

pdf(qq("@{OUTPUT_DIR}/plots/global_methylation_distribution_shore.pdf"))
global_methylation_distribution(sample_id, annotation = SAMPLE$subgroup, annotation_color = COLOR$subgroup,
        chromosome = CHROMOSOME, ha = ha, background = CGI_SHORE, p = 0.01)
dev.off()


pdf(qq("@{OUTPUT_DIR}/plots/global_methylation_distribution_neither_cgi_nor_shore.pdf"))
global_methylation_distribution(sample_id, annotation = SAMPLE$subgroup, annotation_color = COLOR$subgroup,
        chromosome = CHROMOSOME, ha = ha, background = complement)
dev.off()


# cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/01.global_methylation.R")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=10:00:00,mem=10G -N global_methylation' '@{cmd}'")
# system(cmd)

