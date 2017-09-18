 

DAVID_ALL_ID_TYPES = "AFFYMETRIX_3PRIME_IVT_ID,AFFYMETRIX_EXON_GENE_ID,AFFYMETRIX_SNP_ID,AGILENT_CHIP_ID,AGILENT_ID,AGILENT_OLIGO_ID,ENSEMBL_GENE_ID,ENSEMBL_TRANSCRIPT_ID,ENTREZ_GENE_ID,FLYBASE_GENE_ID,FLYBASE_TRANSCRIPT_ID,GENBANK_ACCESSION,GENOMIC_GI_ACCESSION,GENPEPT_ACCESSION,ILLUMINA_ID,IPI_ID,MGI_ID,PFAM_ID,PIR_ID,PROTEIN_GI_ACCESSION,REFSEQ_GENOMIC,REFSEQ_MRNA,REFSEQ_PROTEIN,REFSEQ_RNA,RGD_ID,SGD_ID,TAIR_ID,UCSC_GENE_ID,UNIGENE,UNIPROT_ACCESSION,UNIPROT_ID,UNIREF100_ID,WORMBASE_GENE_ID,WORMPEP_ID,ZFIN_ID"
DAVID_ALL_ID_TYPES = strsplit(DAVID_ALL_ID_TYPES, ",")[[1]]
DAVID_ALL_CATALOGS = "BBID,BIND,BIOCARTA,BLOCKS,CGAP_EST_QUARTILE,CGAP_SAGE_QUARTILE,CHROMOSOME,COG_NAME,COG_ONTOLOGY,CYTOBAND,DIP,EC_NUMBER,ENSEMBL_GENE_ID,ENTREZ_GENE_ID,ENTREZ_GENE_SUMMARY,GENETIC_ASSOCIATION_DB_DISEASE,GENERIF_SUMMARY,GNF_U133A_QUARTILE,GENETIC_ASSOCIATION_DB_DISEASE_CLASS,GOTERM_BP_2,GOTERM_BP_1,GOTERM_BP_4,GOTERM_BP_3,GOTERM_BP_FAT,GOTERM_BP_5,GOTERM_CC_1,GOTERM_BP_ALL,GOTERM_CC_3,GOTERM_CC_2,GOTERM_CC_5,GOTERM_CC_4,GOTERM_MF_1,GOTERM_MF_2,GOTERM_CC_FAT,GOTERM_CC_ALL,GOTERM_MF_5,GOTERM_MF_FAT,GOTERM_MF_3,GOTERM_MF_4,HIV_INTERACTION_CATEGORY,HOMOLOGOUS_GENE,GOTERM_MF_ALL,HIV_INTERACTION,MINT,NCICB_CAPATHWAY_INTERACTION,INTERPRO,KEGG_PATHWAY,PANTHER_FAMILY,PANTHER_BP_ALL,OMIM_DISEASE,OFFICIAL_GENE_SYMBOL,PANTHER_SUBFAMILY,PANTHER_PATHWAY,PANTHER_MF_ALL,PIR_SUMMARY,PIR_SEQ_FEATURE,PFAM,PRODOM,PRINTS,PIR_TISSUE_SPECIFICITY,PIR_SUPERFAMILY,SMART,SP_COMMENT,SP_COMMENT_TYPE,SP_PIR_KEYWORDS,PROSITE,PUBMED_ID,REACTOME_INTERACTION,REACTOME_PATHWAY,UNIGENE_EST_QUARTILE,UP_SEQ_FEATURE,UP_TISSUE,ZFIN_ANATOMY,SSF,TIGRFAMS,UCSC_TFBS"
DAVID_ALL_CATALOGS = strsplit(DAVID_ALL_CATALOGS, ",")[[1]]

# == title
# Doing DAVID analysis
#
# == param
# -genes a vector of genes
# -email the email user registered on David web service
# -catalog a vector of functional catelogs
# -idtype id types
# -species species if the id type is not unique to a species
# 
submit_to_david = function(genes, email, catalog = c("GOTERM_CC_FAT", "GOTERM_BP_FAT", "GOTERM_MF_FAT", "KEGG_PATHWAY"),
	idtype = "ENSEMBL_GENE_ID", species = "Homo sapiens") {

	if(missing(email)) {
		stop("You need to register to DAVID web service.")
	}

	if(!idtype %in% DAVID_ALL_ID_TYPES) {
		cat("idtype is wrong, it should be in:\n")
		print(DAVID_ALL_ID_TYPES)
		stop("You have an error")
	}
	if(!all(catalog %in% DAVID_ALL_CATALOGS)) {
		cat("catalog is wrong, it should be in:\n")
		print(DAVID_ALL_CATALOGS)
		stop("You have an error")
	}

	if(grepl("ENSEMBL", idtype)) {
		map = structure(genes, names = gsub("\\.\\d+$", "", genes))
		genes = names(map)
	}
	
	DAVID_DWS = "https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint"

	# login
	message(qq("log to DAVID web service with @{email}"))
	response = GET(qq("@{DAVID_DWS}/authenticate"),
		query = list("args0" = email)
	)
	if(xml_text(httr::content(response)) != "true") {
		stop(qq("@{email} has not been registed."))
	}

	# add gene list
	message(qq("add gene list (@{length(genes)} genes)"))
	if(length(genes) > 2500) {
		genes = sample(genes, 2500)
	}
	response = POST(qq("@{DAVID_DWS}/addList"),
		body = list("args0" = paste(genes, collapse = ","),  # inputIds
			         "args1" = idtype,             # idType
			         "args2" = Sys.time(),                    # listName
			         "args3" = 0))                           # listType

	message("set species")
	response = GET(qq("@{DAVID_DWS}/getSpecies"))
	all_species = sapply(xml_children(httr::content(response)), xml_text)
	if(length(all_species) > 1) {
		i = grep(species, all_species)
		if(length(i) != 1) {
			cat("check your species, mapped species are:\n")
			print(all_species)
			stop("you have an error.")
		}
		GET(qq("@{DAVID_DWS}/getSpecies"),
			query = list("arg0" = i))
	}

	message("set catalogs")
	response = GET(qq("@{DAVID_DWS}/setCategories"),
		query = list("args0" = paste(catalog, collapse = ","))
	)

	message(qq("doing enrichment"))
	response = GET(qq("@{DAVID_DWS}/getTermClusterReport"),
		query = list("args0" = 3,         # overlap, int
			         "args1" = 3,         # initialSeed, int
			         "args2" = 3,         # finalSeed, int
			         "args3" = 0.5,         # linkage, double
			         "args4" = 1))        # kappa, int

	message(qq("formatting results"))
	xml = httr::content(response)
	clusters = xml_children(xml)
	lt = lapply(clusters, function(x) {
		terms = xml_children(x)[-(1:2)]
		lt = lapply(terms, function(t) {
			fileds = xml_children(t)
			field_name = sapply(fileds, xml_name)
			field_value = sapply(fileds, xml_text)
			l = !field_name %in% c("scores", "listName")
			field_name = field_name[l]
			field_value = field_value[l]
			names(field_value) = field_name
			return(field_value)
		})
		do.call("rbind", lt)
	})
	for(i in seq_along(lt)) {
		lt[[i]] = cbind(lt[[i]], cluster = i)
	}
	tb = do.call("rbind", lt)
	tb = as.data.frame(tb, stringsAsFactors = FALSE)
	for(i in c(1, 2, 3, 4, 6, 7, 8, 11, 12, 13, 14, 15, 16, 18)) {
		tb[[i]] = as.numeric(tb[[i]])
	}

	gene_ids = lapply(strsplit(tb$geneIds, ", "), function(x) map[x])
	tb$geneIds = gene_ids
	return(tb)
}

# == title
# reduce david results
#
# == param
# -tb object from `submit_to_david`
# -fdr_cutoff cutoff for fdr
# -hit_cutoff cutoff for number of genes in a term
#
reduce_david_results = function(tb, fdr_cutoff = 0.05, hit_cutoff = 5) {
	l = tb$benjamini <= fdr_cutoff & tb$listHits >= hit_cutoff
	if(sum(l) < 1) {
		cat("No record left.\n")
		return(NULL)
	}
	tb = tb[l, , drop = FALSE]

	tb_reduced = do.call("rbind", lapply(split(tb, tb$cluster), function(x) {
		do.call("rbind", lapply(split(x, x$categoryName), function(y) {
			i = which.min(y$benjamini)
			y[i, , drop = FALSE]
		}))
	}))

	tb_reduced = tb_reduced[, names(tb_reduced) != "cluster"]
	attr(tb_reduced, "reduced") = TRUE
	return(unique(tb_reduced))
}

# heatmap_david = function(tb, gene_list, email, group_color = NULL, fdr_cutoff = 0.05,
# 	hit_cutoff = 5) {
# 	if(!is.null(attr(tb, "reduced"))) {
# 		stop("`tb` should be directly from `submit_to_david()`.")
# 	}
# 	if(is.null(names(gene_list))) {
# 		stop("`gene_list` must be a named list.")
# 	}
# 	group_names = names(gene_list)
# 	if(is.null(group_color)) {
# 		group_color = structure(brewer.pal(length(group_names), "Set1"), names = group_names)
# 	}

# 	# all records regardless whether it is significant or not
# 	df_all_list = lapply(gene_list, function(genes) {
# 		Sys.sleep(5)
# 		submit_to_david(genes, email)
# 	})
# 	df_list = lapply(df_all_list, function(x) {
# 		reduce_david_results(x, fdr_cutoff = fdr_cutoff, hit_cutoff = hit_cutoff)
# 	})
# 	for(i in names(gene_list)) {
# 		df_list[[i]] = cbind(df_list[[i]], group = as.numeric(i))
# 	}
# 	term_df = do.call("rbind", df_list)

# 	all_genes = unique(unlist(gene_list))
# 	unique_terms = unique(term_df$termName)
# 	gene_mat = matrix(0, nrow = length(all_genes), ncol = length(unique_terms))
# 	rownames(gene_mat) = all_genes
# 	colnames(gene_mat) = unique_terms
# 	for(i in seq_len(nrow(term_df))) {
# 		g_vector = term_df[, "geneIds"][[i]]
# 		gene_mat[g_vector, term_df[i, "termName"]] = 1
# 	}

# 	ht_list = Heatmap(km4, name = "cluster", col = group_mean_col, show_row_names = FALSE, width = 2*grobHeight(textGrob("1", gp = gpar(fontsize = 10))))
# 	if(!is.null(expr_split)) {
# 		ht_list = ht_list + Heatmap(expr_split, name = "expr_direction", show_row_names = FALSE, width = unit(5, "mm"), 
# 			col = c("high" = "Red", "low" = "darkgreen"))
# 	}
# 	unique_category_name = unique(term_df$categoryName)
# 	category_col = structure(brewer.pal(length(unique_category_name), "Set1"), names = unique_category_name)
# 	ht_list = ht_list + Heatmap(gene_mat, name = "identity", col = c("0" = "white", "1" = "blue"), show_row_names = FALSE, 
# 			cluster_columns = FALSE, combined_name_fun = NULL, column_names_max_height = max_text_width(colnames(gene_mat), gp = gpar(fontsize = 12)),
# 			use_raster = TRUE, raster_quality = 2, show_row_dend = FALSE, split = combined_split, show_heatmap_legend = FALSE,
# 			top_annotation = HeatmapAnnotation(categoryName = term_df$categoryName, col = list(categoryName = category_col), show_annotation_name = TRUE))
# 	draw(ht_list, main_heatmap = "identity", row_order = row_order, cluster_rows = FALSE)

	
# 	for(i in seq_along(df_list))	
# 		for(j in seq_len(ncol(gene_mat))) {
# 			current_term_list = df_list[[i]]$Term
# 			ind = which(gsub("GO:.*~", "", df_list[[i]]$Term) == term_names[j])
# 			if(length(ind)) {
# 				message(qq("@{term_names[j]} is significant in @{names(df_list)[i]}"))
# 				decorate_heatmap_body("identity", slice = i, {
# 					grid.rect(x = j/n_term, width = 1/n_term, default.units = "native", just = "right", gp = gpar(fill = "transparent", lwd = 2))
# 					grid.text(paste(df_list[[i]][ind, "listHits"], sprintf("%.1e", df_list[[i]][ind, "benjamini"]), sep = ", "), 
# 						x = (j-0.5)/n_term, y = unit(1, "npc") - unit(2, "mm"), rot = 90, just = "right", gp = gpar(fontsize = 8))
# 				})
# 			}
# 		}
# 	}
# 	draw(ht_list, main_heatmap = "GO")

# 	return(invisible(res_list))
# }