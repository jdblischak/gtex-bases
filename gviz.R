#!/usr/bin/env Rscript

# Experiment with Gviz package for visualizing genomic data
#
# https://bioconductor.org/packages/release/bioc/html/Gviz.html

suppressPackageStartupMessages({
  library(data.table)
  library(Gviz)
})

exons <- fread("data/exons.txt")
setkey(exons, ensembl_gene_id)

track_genes <- GeneRegionTrack(
  rstarts = exons$exon_chrom_start,
  rends = exons$exon_chrom_end,
  exon = exons$ensembl_exon_id,
  strand = exons$strand,
  transcript = exons$ensembl_transcript_id,
  gene = exons$ensembl_exon_id,
  symbol = exons$external_gene_name,
  chromosome = exons$chromosome_name,
  genome = "hg38",
  name = "Gene models",
  # Plot settings
  transcriptAnnotation = "symbol"
)

exons["ENSG00000104904", list(start = min(exon_chrom_start),
                              end = max(exon_chrom_end))]
plotTracks(
  track_genes,
  from = 2269509 - 1000,
  to = 2273490 + 1000,
  chromosome = "chr19"
)

counts <- fread("data/counts/ENSG00000104904-OAZ1-chr19-2270291-2273490.txt")
counts_only <- counts[, -(GeneID:Strand)]
samples <- colnames(counts_only)
samples <- sub(pattern = "^SRR", replacement = "", samples)
samples <- sub(pattern = "\\.bam$", replacement = "", samples)
ymin <- 0
ymax <- max(counts_only)

list_track_counts <- list()
for (i in 1:11) {
  track_counts <- DataTrack(
    start = counts$Start,
    end = counts$End,
    data = counts_only[, i, with = FALSE],
    strand = counts$Strand,
    chromosome = unique(counts$Chr),
    genome = "hg38",
    name = samples[i],
    # Plot settings
    type = "histogram",
    ylim = c(ymin, ymax),
    title = samples[i],
    # background.title="red",
    cex.title = 0.75,
    rotation.title = 0
  )
  list_track_counts <- c(list_track_counts, track_counts)
}

png("ENSG00000104904.png", width = 10, height = 10, units = "in", res = 120)
plotTracks(
  c(list_track_counts, track_genes),
  from = 2269509 - 1000,
  to = 2273490 + 1000,
  chromosome = "chr19",
  main = "ENSG00000104904"
)
dev.off()
