library(GenomicFeatures)
library(rtracklayer)
library(GetoptLong)

source("/home/guz/project/development/epik/R/import_gencode.R")

PROJECT_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/epik_roadmap"
g10 = qq("@{PROJECT_DIR}/data//gen10.long.chr21.gtf")

import_gencode_as_txdb(g10)
import_gencode_as_txdb(g10, transcript_type == "protein_coding")

g19 = qq("@{PROJECT_DIR}/data/gencode.v19.annotation.chr21.gtf")
match_by_gencode(g10, g19, transcript_type == "protein_coding")


available_gencode_fields(g10, "gene")
available_gencode_fields(g10, "transcript")
extract_field_from_gencode(g10, level = "gene", primary_key = "gene_id", field = "gene_name")
extract_field_from_gencode(g10, level = "transcript", primary_key = "transcript_id", field = "transcript_type")
