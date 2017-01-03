
suppressPackageStartupMessages(library(GetoptLong))
min1 = 0.5
min2 = 0.5
window = NULL
GetoptLong("config=s", "configuration R script",
	       "window=i", "window size",
	       "min1=f", ">= floor(min1*n1)",
	       "min2=f", ">= floor(min2*n2)")

library(epic)
load_config(config)

sample_id = rownames(SAMPLE)
chromHMM_list = lapply(sample_id, function(sid) {
	oe = try(gr <- chipseq_hooks$chromHMM(sid))
	if(inherits(oe, "try-error")) {
		return(NULL)
	} else {
		return(gr)
	}
})

l = !sapply(chromHMM_list, is.null)
sample_id = sample_id[l]
chromHMM_list = chromHMM_list[l]
class = SAMPLE$class[l]

if(length(sample_id)) stop("no sample detected.")
if(length(unique(class)) <= 1) stop("less than 2 classes.")

all_classes = unique(class)
nc = length(all_classes)

mat_list = list()
for(i in 1:(nc-1)) {
	for(j in (i+1):nc) {
		gr_list_1 = chromHMM_list[class == all_classes[i]]
		gr_list_2 = chromHMM_list[class == all_classes[j]]
		min1 = min1 * length(gr_list_1)
		min2 = min2 * length(gr_list_2)

		mat = make_transition_matrix_from_chromHMM(gr_list_1, gr_list_2, window = window, min_1 = min1, min_2 = min_2)
		mat_list[[paste0(i, "_", j)]] = mat
	}
}

im = which.max(sapply(mat_list, sum))[1]
max_mat1 = mat_list[[im]]

im = which.max(sapply(mat_list, function(m) {
	cn = intersect(rownames(m), colnames(m))
	for(i in cn) {
		m[i, i] = 0
	}
	sum(m)
}))[1]
max_mat2 = mat_list[[im]]


pdf(qq("@{OUTPUT_DIR}/chromatin_states_transitions.pdf", width = 8, height = 8))
for(i in 1:(nc-1)) {
	for(j in (i+1):nc) {

		if(nrow(mat) > 6) {
			legend_position = c("bottomleft", "bottomright")
		} else {
			legend_position = "bottomleft"
		}

		chromatin_states_transition_chord_diagram(mat, max_mat = max_mat1, remove_unchanged_transition = FALSE, legend_position = legend_position)
		text(1, -1, all_classes[i], adj = c(1, 0))
		text(-1, 1, all_classes[j], adj = c(0, 1))

		chromatin_states_transition_chord_diagram(mat, max_mat = max_mat2, remove_unchanged_transition = TRUE, legend_position = legend_position)
		text(1, -1, all_classes[i], adj = c(1, 0))
		text(-1, 1, all_classes[j], adj = c(0, 1))
	}
}
dev.off()
