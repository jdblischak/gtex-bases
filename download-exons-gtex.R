#!/usr/bin/env Rscript

dir.create("data/", showWarnings = FALSE)

# Download GENCODE 26 annotation file
url_gencode <- "https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf"
file_gencode <- "data/gencode.v26.GRCh38.genes.gtf"
if (!file.exists(file_gencode)) {
  message("Downloading GENCODE annotation")
  download.file(url = url_gencode, dest = file_gencode)
}

