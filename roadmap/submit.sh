for rmd in `ls ~/project/development/epik/roadmap/*.Rmd`; do
	rmd_base=`basename $rmd`
	qsub_single_line "-l walltime=20:00:00,mem=10G -N knit_$rmd_base" "Rscript-3.3.1 -e \"setwd('~/project/development/epik/roadmap');knitr::knit('$rmd_base')\""
done

for chr in `seq 1 22`; do
	qsub_single_line "-l walltime=40:00:00,mem=10G -N generate_meth_dataset_chr$chr" "Rscript-3.3.1 ~/project/development/epik/roadmap/generate_meth_dataset.R --chr chr$chr"
done


for chr in `seq 1 22`; do
	qsub_single_line "-l walltime=10:00:00,mem=10G -N cr_chr$chr" "Rscript-3.3.1 ~/project/development/epik/roadmap/cr.R --chr chr$chr"
done
