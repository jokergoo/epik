

library(GetoptLong)
library(GlobalOptions)
library(parallel)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(matrixStats)
library(GenomeInfoDb)
library(epik.Gviz)
library(HilbertCurve)
library(gtrellis)
library(numbers)
library(ggplot2)
library(Rcpp)
library(memoise) 
library(RColorBrewer)
library(grid)
library(gridBase)
library(rtracklayer)
library(circlize)


if(grepl("tbi", Sys.info()["nodename"]) & Sys.info()["user"] == "guz") {
	Rfiles = list.files("~/project/development/epik/R", full.names = TRUE)
	cpp_files = list.files("~/project/development/epik/src", pattern = "\\.cpp$", full.names = TRUE)
} else {
	Rfiles = list.files("~/project/epik/R", full.names = TRUE)
	cpp_files = list.files("~/project/epik/src", pattern = "\\.cpp$", full.names = TRUE)
}
for(rf in Rfiles) {
	source(rf)
}

for(cf in cpp_files) {
	sourceCpp(cf)
}
