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
setkey(exons, gene)

dir_counts <- "data/counts/"
dir_plots <- "data/plots"
dir.create(dir_plots, showWarnings = FALSE)
files <- list.files(path = dir_counts, pattern = "txt$", full.names = TRUE)

for (f in files) {
  print(f)
  counts <- fread(f)
  plotfile <- basename(f)
  plotfile <- tools::file_path_sans_ext(plotfile)
  plotfile <- paste0(plotfile, ".png")
  plotfile <- file.path(dir_plots, plotfile)
  totals <- counts[, list(total = base::sum(.SD)),
                   .SDcols = grep("bam$", colnames(counts), value = TRUE),
                   by = Start]
  png(plotfile)
  plot(totals$Start, totals$total, main = unique(counts$GeneID))
  dev.off()
}
