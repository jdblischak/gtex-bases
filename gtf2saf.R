#!/usr/bin/env Rscript

# Convert the GTF of all gene models to per-gene SAF files with one entry for
# each base pair.
#
# Usage: Rscript gtf2saf.R
#
# Input:
#   GTF: data/gencode.v26.GRCh38.genes.gtf
#   Target genes: data/target-genes.txt
#
# Output:
#   SAF: data/saf/<ensembl-gene-id>.saf

# Setup ------------------------------------------------------------------------

suppressPackageStartupMessages({
library(data.table)
library(GenomicRanges)
library(rtracklayer)
})

gtf <- "data/gencode.v26.GRCh38.genes.gtf"
target_genes_file <- "data/target-genes.txt"
stopifnot(file.exists(gtf), file.exists(target_genes_file))

dir.create("data/saf/", showWarnings = FALSE)

# Import and format GTF file ---------------------------------------------------

gr <- import(gtf, format = "gtf")

# SAF columns:
# GeneID Chr Start End Strand Name

saf <- data.table(
  GeneID = mcols(gr)[, "gene_id"],
  Chr = as.character(seqnames(gr)),
  Start = start(gr),
  End = end(gr),
  Strand = as.character(strand(gr)),
  Name = mcols(gr)[, "gene_name"]
)

# Remove "chr" from chromosome
saf$Chr <- substr(saf$Chr, 4, nchar(saf$Chr))
# Only keep canonical ID
saf$GeneID <- substr(saf$GeneID, 1, 15)
# Remove Y chromosome. Duplicates entries from X
saf <- saf[Chr != "Y", ]

stopifnot(saf$Start <= saf$End)

str(saf)

# Reduce to one entry per gene
setkey(saf, GeneID)
saf_gene <- saf[, list(Start = min(Start), End = max(End)),
                by = list(GeneID, Chr, Strand, Name)]
setkey(saf_gene, GeneID)
stopifnot(nrow(saf_gene) == length(unique(saf_gene$GeneID)))

str(saf_gene)

# Subset to target genes -------------------------------------------------------

target_genes <- scan(target_genes_file, what = "character")
saf_target_genes <- saf_gene[target_genes]
stopifnot(nrow(saf_target_genes) == length(target_genes))
if (any(is.na(saf_target_genes)))
  stop("Unavailable Ensembl gene IDs requested")

for (i in seq_len(nrow(saf_target_genes))) {
  saf_sub <- saf_target_genes[i, ]
  gene_id <- saf_sub$GeneID

  # Expand to one entry per base pair
  saf_base <- saf_sub[, list(Start = seq(Start, End)),
                       by = list(GeneID, Chr, Strand, Name)]
  saf_base[, End := Start]

  str(saf_base)
  stopifnot(saf_base$End == saf_base$Start)

  saf_base <- saf_base[, list(GeneID, Chr, Start, End, Strand, Name)]
  saf_outfile <- paste0("data/saf/", gene_id, ".saf")
  write.table(saf_base, file = saf_outfile,
              sep = "\t", quote = FALSE, row.names = FALSE)
}
