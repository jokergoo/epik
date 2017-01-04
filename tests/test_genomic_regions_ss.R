context("test gr_ss")

if(Sys.getenv("IS_PBS") != "") {

BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/"

chipseq_hooks$peak = function(mark, sid) {
    qqcat("reading peaks: @{sid}, @{mark}\n")
    df = read.table(qq("@{BASE_DIR}/data/narrow_peaks/@{sid}-@{mark}.narrowPeak.gz"), stringsAsFactors = FALSE)
    GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]]+1, df[[3]]), density = df[[5]])
}

sample_id = c("E003", "E004", "E005", "E006", "E007", "E011", "E012", "E013", "E016",
	          "E024", "E050", "E065", "E066", "E071", "E079", "E094", "E095", "E096",
	          "E097", "E098", "E100", "E104", "E105", "E106", "E109", "E112", "E113")
peak_list = get_peak_list("H3K4me1", sample_id)
subgroup = c(rep("group1", 10), rep("group2", 17))

cr = common_regions(peak_list)
cr2 = common_regions(peak_list, gap = 0.4, min_width = 10000, min_coverage = 7)

res = subgroup_specific_genomic_regions(cr, subgroup = subgroup)
res = subgroup_specific_genomic_regions(cr, subgroup = subgroup, type = c("10", "01"))
res = subgroup_specific_genomic_regions(cr, subgroup = subgroup, present = function(x) sum(x > 0)/length(x) > 0.5, absent = function(x) sum(x > 0)/length(x) < 0.3)
res = subgroup_specific_genomic_regions(cr, subgroup = subgroup, present = 0.6, absent = 0.2)

subgroup_specificity_heatmap(res)

library(GenomicFeatures)
txdb = loadDb(qq("@{BASE_DIR}/data/gen10.long.sqlite"))
genomic_features = list(
	gene = genes(txdb),
	exon = exons(txdb)
)
strand(genomic_features[[1]]) = "*"
strand(genomic_features[[2]]) = "*"

heatmap_subgroup_specificity(res, genomic_features = genomic_features)

}
