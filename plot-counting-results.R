#!/usr/bin/env Rscript

# Setup ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(Gviz)
})

exons_file <- "data/exons.txt"
stopifnot(file.exists(exons_file))

dir_counts <- "data/counts/"
dir_plots <- "data/plots/"
dir.create(dir_plots, showWarnings = FALSE)
files <- list.files(path = dir_counts, pattern = "txt$", full.names = TRUE)

# Genes track ------------------------------------------------------------------

exons <- fread(exons_file)
setkey(exons, ensembl_gene_id)

# Convert chromosome names from Ensembl to UCSC (required by Gviz)
exons[, chromosome_name := paste0("chr", chromosome_name)]
exons[, chromosome_name := sub("chrMT", "chrM", chromosome_name)]

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

# Counts tracks ----------------------------------------------------------------

for (f in files) {
  print(f)
  counts <- fread(f)
  gene_id <- unique(counts$GeneID)
  # Convert chromosome names from Ensembl to UCSC (required by Gviz)
  counts[, Chr := paste0("chr", Chr)]
  counts[, Chr := sub("chrMT", "chrM", Chr)]
  chromosome = unique(counts$Chr)
  start <- min(counts$Start)
  end <- max(counts$End)
  width <- end - start

  counts_only <- counts[, -(GeneID:Strand)]
  samples <- colnames(counts_only)
  samples <- sub(pattern = "^SRR", replacement = "", samples)
  samples <- sub(pattern = "\\.bam$", replacement = "", samples)
  ymin <- 0
  ymax <- max(counts_only)

  list_track_counts <- list()
  for (i in seq_along(samples)) {
    track_counts <- DataTrack(
      start = counts$Start,
      end = counts$End,
      data = counts_only[, i, with = FALSE],
      strand = counts$Strand,
      chromosome = chromosome,
      genome = "hg38",
      name = samples[i],
      # Plot settings
      type = "histogram",
      ylim = c(ymin, ymax),
      title = samples[i],
      # background.title="red",
      cex.title = 0.75,
      rotation.title = 90
    )
    list_track_counts <- c(list_track_counts, track_counts)
  }

  plotfile <- basename(f)
  plotfile <- tools::file_path_sans_ext(plotfile)
  plotfile <- paste0(plotfile, ".png")
  plotfile <- file.path(dir_plots, plotfile)

  png(plotfile, width = 10, height = 10, units = "in", res = 120)
  plotTracks(
    c(list_track_counts, track_genes),
    from = start,
    to = end,
    extend.left = width * 0.1,
    extend.right = width * 0.1,
    chromosome = chromosome,
    main = gene_id
  )
  dev.off()
}
