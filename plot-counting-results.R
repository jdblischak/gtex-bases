#!/usr/bin/env Rscript

library(data.table)
library(GenomicRanges)
library(rtracklayer)

gtf <- "data/gencode.v26.GRCh38.genes.gtf"
gr <- import(gtf, format = "gtf")

exons <- data.table(
  gene = mcols(gr)[, "gene_id"],
  exon = mcols(gr)[, "exon_id"],
  chr = as.character(seqnames(gr)),
  start = start(gr),
  end = end(gr),
  strand = as.character(strand(gr)),
  name = mcols(gr)[, "gene_name"]
)
exons[, gene := stringr::str_replace(gene, "(ENSG\\d+)\\.\\d+", "\\1")]
setkey(exons, gene)

dir_counts <- "data/counts/"
dir_plots <- "data/plots"
dir.create(dir_plots, showWarnings = FALSE)
files <- list.files(path = dir_counts, pattern = "txt$", full.names = TRUE)

for (f in files) {
  print(f)
  counts <- fread(f)
  gene_id <- unique(counts$GeneID)

  plotfile <- basename(f)
  plotfile <- tools::file_path_sans_ext(plotfile)
  plotfile <- paste0(plotfile, ".png")
  plotfile <- file.path(dir_plots, plotfile)

  totals <- counts[, list(total = base::sum(.SD)),
                   .SDcols = grep("bam$", colnames(counts), value = TRUE),
                   by = Start]

  gene_exons <- exons[gene_id]
  gene_exons <- gene_exons[!is.na(exon), ]

  png(plotfile)
  plot(totals$Start, totals$total, main = gene_id, ylim = c(-50, 1000))
  for (i in seq_len(nrow(gene_exons))) {
    lines(c(gene_exons$start[i], gene_exons$end[i]), c(-50, -50), col = "red", lwd = 5)
  }
  dev.off()
}
