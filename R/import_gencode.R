

import_gencode_as_txdb = function(gtf, filter = NULL, include_mappings = FALSE) {

	qqcat("import @{gtf} as a GRanges object.\n")
	gr = import(gtf)
	l = eval(substitute(filter), envir = gr@elementMetadata@listData)
	if(!is.null(l)) {
		gr = gr[l]
	}
	cat("Making txdb\n")
	txdb = makeTxDbFromGRanges(gr)

	if(include_mappings) {
		cat("Generate mappings\n")
		l1 = !duplicated(gr$gene_id)
		l2 = !duplicated(gr$transcript_id)

		gene_id2name = NULL
		if(!is.null(gr$gene_name)) gene_id2name = structure(gr$gene_name[l1], names = gr$gene_id[l1])
		
		gene_id2type = NULL
		if(!is.null(gr$gene_type)) gene_id2type = structure(gr$gene_type[l1], names = gr$gene_id[l1])
		
		tx_id2name = NULL
		if(!is.null(gr$transcript_name)) tx_id2name = structure(gr$transcript_name[l2], names = gr$transcript_id[l2])
		
		tx_id2type = NULL
		if(!is.null(gr$transcript_type)) tx_id2type = structure(gr$transcript_type[l2], names = gr$transcript_id[l2])
		
		mapping = list(gene_id2name = gene_id2name,
			           gene_id2type = gene_id2type,
			           tx_id2name = tx_id2name,
			           tx_id2type = tx_id2type)

		return(list(txdb = txdb, mapping = mapping))
	} else {
		return(txdb)
	}
}
