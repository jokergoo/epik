library(GetoptLong)
source('/home/guz/project/development/epik/R/cpg_dinucleotide_methylation.R')

pos = c(1, 2, 6, 7, 8, 12, 15, 18, 19, 20, 21)

meth = cbind(c(0.9, 0.9, 0.4, 0.5, 0.5, 0.7, 0.8, 0.3, 0.5, 0.7, 0.3),
	         c(0.8, 0.8, 0.5, 0.7, 0.7, 0.2, 0.4, 0.3, 0.5, 0.7, 0.4))

cov = cbind(sample(20, 11), sample(20, 11))
meth = round(meth*cov)

lt = cpg_dinucleotide_methylation(pos, meth, cov)


load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/wgbs_gbm_smoothed/smoothed_object_list_chr21.RData")
lt = cpg_dinucleotide_methylation(bs.fit$site, bs.fit$raw, bs.fit$coverage)
