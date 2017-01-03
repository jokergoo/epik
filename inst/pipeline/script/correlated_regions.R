
suppressPackageStartupMessages(library(GetoptLong))
GetoptLong("config=s", "configuration R script",
	       "chr=s", "chromosomes")

library(epic)
load_config(config)

sample_id = rownames(SAMPLE)

if(is.null(EXPR)) stop("'EXPR' should be defined.")
if(is.null(TXDB)) stop("'TXDB' should be defined.")

if(!chr %in% CHROMOSOME) stop("'chr' should be included in 'CHROMOSOME'.")

cr = correlated_regions(sample_id, EXPR, TXDB, chr = chr, factor = SAMPLE$class, col = COLOR$class)
saveRDS(cr, file = qq("@{OUTPUT_DIR}/rds/@{chr}_cr.rds"))

# if(raw == "yes") {
# 	cr_raw = correlated_regions(sample_id, EXPR, TXDB, raw_meth = TRUE, chr = chr, factor = SAMPLE$class, 
# 		col = COLOR$class, cov_cutoff = 5, min_dp = max(c(5, round(length(sample_id)/2))))
# 	saveRDS(cr_raw, file = qq("@{RDS_FOLDER}/cr_raw_meth_@{cr}.rds"))
# }
