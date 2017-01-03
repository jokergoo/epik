
suppressPackageStartupMessages(library(GetoptLong))

cutoff = NULL
peak = NULL
GetoptLong("config=s", "configuration R script",
	       "cutoff=f", "cutoff for filter cr",
	       "peak=s", "name of peak",
	       "which=s", "neg|pos")

library(epic)
load_config(config)

if(!peak %in% MARKS) {
	stop("'peak' should be in 'MARKS'.")
}
if(!which %in% c("neg", "pos")) {
	stop("'which' should be in c('neg', 'pos').")
}


if(is.null(cutoff)) {
	cutoff = CR_CUTOFF
}

cr_filtered = readRDS(qq("@{OUTPUT_DIR}/rds/cr_filtered_fdr_@{cutoff}.rds"))

make_enriched_heatmap = function(histone_mark, by, on, which, type = 1:3) {
	
	hm_list = get_peak_list(histone_mark)

	if(which == "neg") {
	    cr = cr_filtered[cr_filtered$corr < 0]
	} else {
	    cr = cr_filtered[cr_filtered$corr > 0]
	}

	n_factor = length(unique(attr(cr, "factor")[attr(cr, "sample_id") %in% names(hm_list)]))
	if(n_factor == 0) n_factor = 1

	ht_global_opt(heatmap_column_title_gp = gpar(fontsize = 10))
	if(any(type %in% c(1, 2, 3))) {
		qqcat("heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}\n")
		pdf(qq("@{OUTPUT_DIR}/enriched_heatmap/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf"), width = 8 + n_factor*2, height = 10)
		enriched_heatmap_list_on_gene(cr, GENOMIC_FEATURE_LIST$cgi, TXDB, EXPR, hm_list, 
		    hm_name = histone_mark, on = on, by = by)
		dev.off()

		system(qq("convert @{OUTPUT_DIR}/enriched_heatmap/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf @{OUTPUT_DIR}/enriched_heatmap/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.png"))
	}

	if(any(type %in% c(3, 2))) {
		qqcat("heatmap_cgi_@{by}_@{which}_@{histone_mark}\n")
		pdf(qq("@{OUTPUT_DIR}/enriched_heatmap/heatmap_cgi_@{by}_@{which}_@{histone_mark}.pdf"), width = 6 + n_factor*2, height = 10)
		enriched_heatmap_list_on_tss_cgi(cr, GENOMIC_FEATURE_LIST$cgi, TXDB, EXPR, hm_list, 
		    hm_name = histone_mark, by = by)
		dev.off()

		system(qq("convert @{OUTPUT_DIR}/enriched_heatmap/heatmap_cgi_@{by}_@{which}_@{histone_mark}.pdf @{OUTPUT_DIR}/enriched_heatmap/heatmap_cgi_@{by}_@{which}_@{histone_mark}.png"))
	}

	if(any(type %in% 3)) {
		if("tfbs" %in% names(GENOMIC_FEATURE_LIST)) {
			qqcat("heatmap_encode_tfbs_@{which}_@{histone_mark}")
			pdf(qq("@{OUTPUT_DIR}/enriched_heatmap/heatmap_encode_tfbs_@{which}_@{histone_mark}.pdf"), width = 6 + n_factor*2, height = 10)
			enriched_heatmap_list_on_genomic_features(cr, GENOMIC_FEATURE_LIST$tfbs, hm_list, 
			    hm_name = histone_mark)
			dev.off()

			system(qq("convert @{OUTPUT_DIR}/enriched_heatmap/heatmap_encode_tfbs_@{which}_@{histone_mark}.pdf @{OUTPUT_DIR}/enriched_heatmap/heatmap_encode_tfbs_@{which}_@{histone_mark}.png"))
		}

		if("enhancer" %in% names(GENOMIC_FEATURE_LIST)) {
			qqcat("heatmap_encode_strong_enhancer_@{which}_@{histone_mark}")
			pdf(qq("@{OUTPUT}/enriched_heatmap/heatmap_encode_strong_enhancer_@{which}_@{histone_mark}.pdf"), width = 6 + n_factor*2, height = 10)
			enriched_heatmap_list_on_genomic_features(cr, GENOMIC_FEATURE_LIST$enhancer, hm_list, 
			    hm_name = histone_mark)
			dev.off()

			system(qq("convert @{OUTPUT_DIR}/enriched_heatmap/heatmap_encode_strong_enhancer_@{which}_@{histone_mark}.pdf @{OUTPUT_DIR}/enriched_heatmap/heatmap_encode_strong_enhancer_@{which}_@{histone_mark}.png"))
		}
	}
}

make_enriched_heatmap(histone_mark, by = "gene", on = "tss", which = which, type = 1:3)
make_enriched_heatmap(histone_mark, by = "tx", on = "tss", which = which, type = 1:2)
make_enriched_heatmap(histone_mark, by = "gene", on = "body", which = which, type = 1)

