library(methods)
suppressPackageStartupMessages(library(GetoptLong))

chr = "chr18"
peak = NULL
GetoptLong("chr=s", "chromosome",
	       "peak=s@", "peak name")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

if(is.null(peak)) peak = MARKS

neg_cr = readRDS(qq("@{OUTPUT_DIR}/rds/all_neg_cr_w6s3.rds"))
pos_cr = readRDS(qq("@{OUTPUT_DIR}/rds/all_pos_cr_w6s3.rds"))
cr = c(neg_cr, pos_cr)
cr = cr[cr$corr_fdr < 0.05 & cr$meth_anova_fdr < 0.05 & cr$meth_diameter > 0.1]

cr = copy_cr_attribute(neg_cr, cr)

if(!chr %in% unique(as.character(seqnames(cr)))) {
	stop("'chr' does not exist in 'cr'.")
}

dir.create(qq("@{OUTPUT_DIR}/gviz/@{chr}"), showWarnings = FALSE)

if(length(peak) > 0) {
	if(!all(peak %in% MARKS)) {
		stop("'peak' should all be in 'MARKS'.")
	}
	peak_list = lapply(peak, get_peak_list)
	names(peak_list) = peak
}

cr_subset = cr[seqnames(cr) == chr]
gn = epic::extract_field_from_gencode(GTF_FILE, level = "gene", primary_key = "gene_id", field = "gene_name")
for(gi in unique(cr_subset$gene_id)) {
    pdf(qq("@{OUTPUT_DIR}/gviz/@{chr}/gviz_@{chr}_@{gi}_@{gn[gi]}.pdf"), width = 12, height = 6 + 1*length(peak))
    # try(
	    cr_gviz(cr_subset, gi, EXPR, TXDB, 
	    	gf_list = list(CGI = CGI), 
	    	hm_list = peak_list, symbol = gn[gi])
   	# )
    dev.off()
}


# for(chr in paste0("chr", 1:22)) {
#     cmd = qq("Rscript-3.3.0 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/04.correlated_regions_gviz.R --chr @{chr}") 
#     cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=15G -N cr_gviz_@{chr}' '@{cmd}'")
#     system(cmd)
# }
