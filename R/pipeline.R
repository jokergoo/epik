 

# == title
# run pipeline in HPC through job scheduling system
#
# == param
# -config_file path of configuration script, check `load_config` for all configurations.
# -prefix prefix of the job name
# -email email address if you want to be notified of your jobs
# -enforce enforce to run all the steps even for those successfully finished jobs.
# -Rscript_binary path of Rscript binary
# -submit_by which job scheduling you are using
#
# == details
# A workflow will be submitted to HPC.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
run_pipeline = function(config_file, prefix = "", email = NULL, enforce = FALSE, Rscript_binary = "Rscript", submit_by = "qsub") {

	OUTPUT_DIR = NULL
	SAMPLE = NULL
	EXPR = NULL
	TXDB = NULL
	CHROMOSOME = NULL
	MARKS = NULL

	load_config(config_file, validate = FALSE)

	tmpdir = paste0(OUTPUT_DIR, "/temp")
	dir.create(tmpdir, showWarnings = FALSE)

	if(length(unique(SAMPLE$class)) > 1) {
		x = pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' differential_methylation_in_cgi_and_shore --config '@{config_file}'"), 
			output = c(qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi.rds}"),
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi_shore.rds}"),
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_neither_cgi_nor_shore.rds}"),
				       qq("@{OUTPUT_DIR}/heatmap_diff_methylation_1kb_window.pdf"),
				       qq("@{OUTPUT}/genome_diff_1kb_window_correlation.pdf")),
			name = qq("@{prefix}differential_methylation_in_cgi_and_shore"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			tmpdir = tmpdir,
			submit_by = submit_by)

		pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' differential_methylation_in_genomic_features --config '@{config_file}'"),
			output = c(qq("@{OUTPUT_DIR}/heatmap_diff_methylation_in_genomic_features.pdf")),
			name = qq("@{prefix}differential_methylation_in_genomic_features"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			tmpdir = tmpdir,
			submit_by = submit_by)

		pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' methylation_subtype_classification_in_cgi_and_shore --config '@{config_file}'"),
			output = c(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi_shore.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_neither_cgi_nor_shores.pdf")),
			name = qq("@{prefix}methylation_subtype_classification_in_cgi_and_shore"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			dependency = x,
			tmpdir = tmpdir,
			submit_by = submit_by)
	} else {
		pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' methylation_subtype_classification_in_cgi_and_shore --config '@{config_file}'"),
			output = c(qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_cgi_shore.pdf"),
				       qq("@{OUTPUT_DIR}/methylation_classification_wgbs_neither_cgi_nor_shores.pdf"),
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi.rds}"), 
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_cgi_shore.rds}"), 
				       qq("@{OUTPUT_DIR/rds/mean_meth_1kb_neither_cgi_nor_shore.rds}")),
			name = qq("@{prefix}methylation_subtype_classification_in_cgi_and_shore"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			tmpdir = tmpdir,
			submit_by = submit_by)
	}

	pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' general_methylation_distribution --config '@{config_file}'"),
		output = c(qq("@{OUTPUT_DIR}/general_methylation_distribution.pdf"),
			       qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi.pdf"),
			       qq("@{OUTPUT_DIR}/general_methylation_distribution_cgi_shores.pdf"),
			       qq("@{OUTPUT_DIR}/general_methylation_distribution_neither_cgi_nor_shores.pdf")),
		name = qq("@{prefix}general_methylation_distribution"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		tmpdir = tmpdir,
		submit_by = submit_by)

	#####################################
	### CR
	#####################################

	if(is.null(EXPR) || is.null(TXDB)) {
		return(NULL)
	}

	dependency = NULL
	for(chr in CHROMOSOME) {
		x = pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions --config '@{config_file}' --chr @{chr}"),
			output = qq("@{OUTPUT_DIR}/rds/@{chr}_cr.rds"),
			name = qq("@{prefix}correlated_regions_@{chr}"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			tmpdir = tmpdir,
			submit_by = submit_by)
		dependency = c(dependency, x)
	}

	pid_cr_filter = pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_filter --config '@{config_file}'"),
		output = qq("@{OUTPUT_FOLDER}/rds/cr_filtered_fdr_@{cutoff}.rds"),
		name = qq("@{prefix}correlated_regions_filter"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		dependency = dependency,
		tmpdir = tmpdir,
		submit_by = submit_by)

	pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_reduce --config '@{config_file}'"),
		output = qq("@{OUTPUT_DIR}/rds/cr_reduced_fdr_@{cutoff}.rds"),
		name = qq("@{prefix}correlated_regions_reduce"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		dependency = pid_cr_filter,
		tmpdir = tmpdir,
		submit_by = submit_by)

	pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_downstream --config '@{config_file}'"),
		output = c(qq("@{OUTPUT_DIR}/cr_number.pdf"),
			       qq("@{OUTPUT_DIR}/hilbert_sig.pdf"),
			       qq("@{OUTPUT_DIR}/hilbert_all.pdf"),
			       qq("@{OUTPUT_DIR}/cr_overlap.pdf"),
			       qq("@{OUTPUT_DIR}/cr_tss.pdf")),
		name = qq("@{prefix}correlated_regions_downstream"),
		walltime = "10:00:00",
		mem = "10G",
		nodes = 1,
		email = email,
		enforce = enforce,
		dependency = pid_cr_filter,
		tmpdir = tmpdir,
		submit_by = submit_by)

	for(chr in CHROMOSOME) {
		pipeline_step(qq("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_gviz --config '@{config_file}' --chr @{chr}"),
			output = NULL,
			name = qq("@{prefix}correlated_regions_gviz"),
			walltime = "10:00:00",
			mem = "10G",
			nodes = 1,
			email = email,
			enforce = enforce,
			dependency = pid_cr_filter,
			tmpdir = tmpdir,
			submit_by = submit_by)
	}

	if(is.null(MARKS)) {
		return(NULL)
	}

	for(mk in MARKS) {
		for(which in c("pos", "neg")) {
			pipeline_step("'@{Rscript_binary}' -e 'epic::epic()' correlated_regions_enriched --config '@{config_file}' --peak @{mk} --which @{which}",
				output = NULL,
				name = qq("@{prefix}correlated_regions_enriched_@{mk}_@{which}"),
				walltime = "10:00:00",
				mem = "10G",
				nodes = 1,
				email = email,
				enforce = enforce,
				dependency = pid_cr_filter,
				tmpdir = tmpdir,
				submit_by = submit_by)
		}
	}
}


pipeline_step = function(..., output = NULL, name, walltime = "1:00:00", mem = "1G", nodes = 1, 
	email = NULL, dependency = NULL, enforce = FALSE, tmpdir, submit_by = "qsub") {
	
	if(is.null(output)) {
		enforce = TRUE
	}
	if(!enforce && all(file.exists(output))) {
		return(NULL)
	}
	cmd = unlist(list(...))
	cmd = paste(cmd, collapse = "\n")
	if(submit_by == "qsub") {
		if(length(dependency) > 1) dependency = paste(dependency, collapse = ",")
		script = qq("#!/bin/sh
#PBS -j oe
#PBS -o @{OUTPUT_DIR}/temp/
@{ifelse(is.null(dependency), '', paste0('#PBS -W depend=afterok:', dependency))}
#PBS -N @{name}
#PBS -l walltime=@{walltime}
#PBS -l mem=@{mem}
#PBS -l nodes=@{nodes}
@{ifelse(is.null(email), '', paste0('#PBS -M', email))}

@{cmd}
")
	} else if(submit_by == "bsub") {
		if(grepl("M", mem)) {
			mem = as.numeric(gsub("M|MB", "", mem, ignore.case = TRUE))
		} else if(grepl("G", mem)) {
			mem = as.numeric(gsub("G|GB", "", mem, ignore.case = TRUE))*1024
		}

		if(!is.null(dependency)) {
			if(length(dependency) > 1) {
				dependency = strsplit(dependency, ",")[[1]]	
			}
			dependency = paste("done(", dependency, ")", sep = "", collapse = " && ")
		}
				script = qq("#!/bin/sh
#BSUB -oo @{OUTPUT_DIR}/temp/@{name}.out
#BSUB -eo @{OUTPUT_DIR}/temp/@{name}.err
#BSUB -J @{name}
#BSUB -W @{walltime}
#BSUB -R rusage[mem=@{mem}]
#BSUB -n @{nodes}
@{ifelse(is.null(dependency), '', paste0('#BSUB -w ', dependency))}
@{ifelse(is.null(email), '', paste0('#BSUB -u', email))}

@{cmd}
")
	}

	temp_file = tempfile(tmpdir = tmpdir, fileext = ".sh")
	writeLines(script, temp_file)

	if(submit_by == "qsub") {
		x = system(qq("qsub '@{temp_file}'"), intern = TRUE)
	} else if(submit_by == "bsub") {
		x = system(qq("bsub < '@{temp_file}'"), intern = TRUE)
	}
	qqcat("submit job: @{name}, id: @{x[1]}\n")
	return(x[1])
}


# x = pipeline_step("ls", output = "a", name = "test", walltime="1:00:00", mem="1G", nodes = 1, email = "z.gu@dkfz.de")
# x = pipeline_step("ls -l", output = "a", name = "test", walltime="1:00:00", mem="1G", nodes = 1, email = "z.gu@dkfz.de", dependency = x)

