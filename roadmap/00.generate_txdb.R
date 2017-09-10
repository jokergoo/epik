library(GenomicFeatures)
library(GetoptLong)
library(epik)
library(epik.cmd)

PROJECT_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/epik_roadmap/"
initialize_project_directory(PROJECT_DIR)

GTF_FILE = qq("@{PROJECT_DIR}/data/gen10.long.gtf")
g19 = qq("@{PROJECT_DIR}/data/gencode.v19.annotation.gtf")
txdb = match_by_gencode(GTF_FILE, g19, transcript_type == "protein_coding")

saveDb(txdb, file = qq("@{PROJECT_DIR}/txdb/gen10_long_protein_coding_gene_adjusted.sqlite"))
