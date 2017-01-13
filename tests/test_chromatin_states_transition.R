source("test_cr_chr21.R")

chipseq_hooks$chromHMM = function(sid, chr) {
	x = qq("@{BASE_DIR}/data/chromatin_states/@{sid}_15_coreMarks_mnemonics.bed.gz")
	qqcat("reading @{x} on @{chr}...\n")
	gr = read.table(pipe(qq("zcat @{x} | grep @{chr}")), sep = "\t")
	gr = gr[gr[[1]] %in% chr, ]
	GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]] + 1, gr[[3]]), states = gr[[4]])
}

states_list = get_chromHMM_list(sample_id, "chr21")

methdiff = 0

gr_list_1 = states_list[subgroup == "group1"]
gr_list_2 = states_list[subgroup == "group2"]

mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2, methylation_diff = methdiff)
mat

sn = state_names(mat)
od = order(as.numeric(gsub("^(\\d+)_.*$", "\\1", sn)))
mat = mat[od, od]
state_names(mat) = gsub("\\d+_", "", state_names(mat))

## color setting and make the plot
state_col = c("TssA" = "#E41A1C",
	          "TssAFlnk" = "#E41A1C",
	          "TxFlnk" = "#E41A1C",
	          "Tx" = "#E41A1C",
	          "TxWk" = "#E41A1C",
	          "EnhG" = "#E41A1C",
	          "Enh" = "#E41A1C",
	          "ZNF/Rpts" = "#E41A1C",
	          "Het" = "#377EB8",
	          "TssBiv" = "#377EB8",
	          "BivFlnk" = "#377EB8",
	          "EnhBiv" = "#377EB8",
	          "ReprPC" = "#377EB8",
	          "ReprPCWk" = "#377EB8",
	          "Quies" = "black")

chromatin_states_transition_chord_diagram(mat, state_col = state_col, group_names = c("g1", "g2"),
	legend_position = c("bottomleft", "bottomright"))



methylation_hooks(RESET = TRUE)
mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2)
mat

sn = state_names(mat)
od = order(as.numeric(gsub("^(\\d+)_.*$", "\\1", sn)))
mat = mat[od, od]
state_names(mat) = gsub("\\d+_", "", state_names(mat))
chromatin_states_transition_chord_diagram(mat, state_col = state_col, group_names = c("g1", "g2"),
	legend_position = c("bottomleft", "bottomright"))

