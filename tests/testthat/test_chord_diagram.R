context("test chromatin transition diagram")

if(Sys.getenv("IS_PBS") != "") {

files = dir("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/chromatin_states", pattern = "bed.gz$")
gr_list_1 = lapply(files[1:5], function(f) {
	df = read.table(paste0("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/chromatin_states/", f),stringsAsFactors = FALSE)
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), states = df[[4]])
})
gr_list_2 = lapply(files[6:10], function(f) {
	df = read.table(paste0("/icgc/dkfzlsdf/analysis/B080/guz/epic_test/data/chromatin_states/", f),stringsAsFactors = FALSE)
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), states = df[[4]])
})

mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2)

states = rownames(mat)
numbers = as.numeric(gsub("^(\\d+)_.*$", "\\1", states))
states = gsub("^\\d+_", "", states)
rownames(mat) = states
colnames(mat) = states
mat = mat[order(numbers), order(numbers)]

chromatin_states_transition_chord_diagram(mat, remove_unchanged_transition = FALSE, legend_position = c("bottomleft", "bottomright"))
chromatin_states_transition_chord_diagram(mat, remove_unchanged_transition = TRUE, legend_position = c("bottomleft", "bottomright"))

}
