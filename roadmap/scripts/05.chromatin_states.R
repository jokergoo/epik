library(methods)
suppressPackageStartupMessages(library(GetoptLong))

methdiff = 0
GetoptLong("methdiff=f", "methylation difference")

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis"
source(qq("@{BASE_DIR}/scripts/configure/roadmap_configure.R"))

subgroup = SAMPLE[sample_id, "subgroup"]

fn = dir(qq("@{BASE_DIR}/data/chromatin_states/"))
nm = gsub("^(E\\d+?)_.*$", "\\1", fn)
fn = fn[nm %in% colnames(EXPR)]

all_sample_chromatin_states = lapply(fn, function(x) {
	x = qq("@{BASE_DIR}/data/chromatin_states/@{x}")
	qqcat("reading @{x}...\n")
	gr = read.table(x, sep = "\t")
	gr = gr[gr[[1]] %in% CHROMOSOME, ]
	GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]] + 1, gr[[3]]), states = gr[[4]])
})

names(all_sample_chromatin_states) = gsub("^(E\\d+)_.*$", "\\1", fn)

gr_list_1 = all_sample_chromatin_states[names(all_sample_chromatin_states) %in% colnames(EXPR[, subgroup == "subgroup1"])]
gr_list_2 = all_sample_chromatin_states[names(all_sample_chromatin_states) %in% colnames(EXPR[, subgroup == "subgroup2"])]

mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2, methylation_diff = methdiff)


## reorder rows and columns in `mat`
meth_mean_1 = attr(mat, "meth_mean_1")
meth_mean_2 = attr(mat, "meth_mean_2")

od1 = order(as.numeric(gsub("^(\\d+)_.*$", "\\1", rownames(mat))))
od2 = order(as.numeric(gsub("^(\\d+)_.*$", "\\1", colnames(mat))))
mat = mat[od1, od2]
rn = gsub("^\\d+_", "", rownames(mat))
cn = gsub("^\\d+_", "", colnames(mat))
rownames(mat) = rn
colnames(mat) = cn

if(!is.null(meth_mean_1)) {
	meth_mean_1 = meth_mean_1[od1, od2]
	meth_mean_2 = meth_mean_2[od1, od2]
	rownames(meth_mean_1) = rn
	colnames(meth_mean_1) = cn
	rownames(meth_mean_2) = rn
	colnames(meth_mean_2) = cn

	attr(mat, "meth_mean_1") = meth_mean_1
	attr(mat, "meth_mean_2") = meth_mean_2
}
saveRDS(mat, file = qq("@{OUTPUT_DIR}/rds/chromatin_states_transition_meth_diff_@{methdiff}.rds"))


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

pdf(qq("@{OUTPUT_DIR}/plots/chromatin_states_transites_meth_diff_@{methdiff}.pdf"), width = 8, height = 8)
chromatin_states_transition_chord_diagram(mat, state_col = state_col, group_names = c("g1", "g2"),
	legend_position = c("bottomleft", "bottomright"))
dev.off()


# cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/05.chromatin_states.R --methdiff 0")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=4:00:00,mem=10G -N chromatin_states_methdiff_0' '@{cmd}'")
# system(cmd)
# cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/05.chromatin_states.R --methdiff 0.1")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=4:00:00,mem=10G -N chromatin_states_methdiff_0.1' '@{cmd}'")
# system(cmd)
# cmd = qq("Rscript-3.1.2 /icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/scripts/05.chromatin_states.R --methdiff 0.2")
# cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=4:00:00,mem=10G -N chromatin_states_methdiff_0.2' '@{cmd}'")
# system(cmd)
