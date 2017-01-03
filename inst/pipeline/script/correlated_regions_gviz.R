
suppressPackageStartupMessages(library(GetoptLong))

cutoff = NULL
peak = NULL
GetoptLong("config=s", "configuration R script",
	       "cutoff=f", "cutoff for filter cr",
	       "chr=s", "chromosome",
	       "peak=s", "peak name")

library(epic)
load_config(config)

if(is.null(cutoff)) {
	cutoff = CR_CUTOFF
}

cr_filtered = readRDS(qq("@{OUTPUT_DIR}/rds/cr_filtered_fdr_@{cutoff}.rds"))

if(!chr %in% unique(as.character(seqnames(cr_filtered)))) {
	stop("'chr' does not exist in 'cr_filtered'.")
}

if(is.null(peak)) {
	peak_list = NULL
	if(!is.null(MARKS)) {
		cat("You can add peak tracks to Gviz plots by manually specifying '--peak' options.\n")
	}
} else {
	if(!peak %in% MARKS) {
		stop("'peak' should be in 'MARKS'.")
	}
	peak_list = get_peak_list(peak)
}

cr_subset = cr_filtered[seqnames(cr_filtered) == chr]
gn = extract_field_from_gencode(GTF_FILE, level = "gene", primary_key = "gene_id", field = "gene_name")
for(gi in unique(cr_subset$gene_id)) {
    pdf(qq("@{OUTPUT_DIR}/gviz/gviz_@{chr}_@{gi}_@{gn[gi]}.pdf"), width = 16, height = 12)
    cr_gviz(cr_subset, gi, EXPR, TXDB, 
    	gf_list = GENOMIC_FEATURE_LIST[intersect(c("cgi", "tfbs", "enhancer"), names(GENOMIC_FEATURE_LIST))], 
    	hm_list = peak_list, symbol = gn[gi])
    dev.off()
}


