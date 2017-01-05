
BASE_DIR = "/icgc/dkfzlsdf/analysis/B080/guz/roadmap_analysis/re_analysis/"
GTF_FILE = qq("@{BASE_DIR}/data//gen10.long.chr21.gtf")

import_gencode_as_txdb(GTF_FILE)
import_gencode_as_txdb(GTF_FILE, gene_type == "protein_coding")
import_gencode_as_txdb(GTF_FILE, include_mappings = TRUE)
lt = import_gencode_as_txdb(GTF_FILE, transcript_type == "protein_coding", include_mappings = TRUE)
