#!/usr/bin/env Rscript

# Create per-gene SAF file with one entry for each base pair.
#
# Usage: Rscript create-saf-files.R <exons-file>
#
# Input:
#   Ensembl exons: data/exons/<ensembl-gene-id>.txt
#
# Output (writes to stdout):
#   SAF: data/saf/<ensembl-gene-id>.saf

# Setup ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  exons_file <- args[1]
} else {
  stop("Usage: Rscript create-saf-file.R <ensembl-gene-id>.txt", call. = FALSE)
}

stopifnot(file.exists(exons_file))

# Import and format exons ------------------------------------------------------

exons <- fread(exons_file)

# SAF columns:
# GeneID Chr Start End Strand Name

saf <- exons[, list(GeneID = ensembl_gene_id,
                    Chr = chromosome_name,
                    Start = exon_chrom_start,
                    End = exon_chrom_end,
                    Strand = strand,
                    Name = external_gene_name)]

# Confirm that the file only contains one gene on one chromosome
stopifnot(
  length(unique(saf$Chr)) == 1,
  length(unique(saf$GeneID)) == 1
)

stopifnot(saf$Start <= saf$End)

# Reduce to a single entry
saf_gene <- saf[, list(Start = min(Start), End = max(End)),
                by = list(GeneID, Chr, Strand, Name)]

# Create per-base features -----------------------------------------------------

# Expand to one entry per base pair
saf_base <- saf_gene[, list(Start = seq(Start, End)),
                     by = list(GeneID, Chr, Strand, Name)]
saf_base[, End := Start]
stopifnot(saf_base$End == saf_base$Start)

# Append base number to each feature, for interpreting assignments from "-R CORE"
saf_base[, GeneID := paste(GeneID, seq_len(nrow(saf_base)), sep = "_")]


saf_base <- saf_base[, list(GeneID, Chr, Start, End, Strand, Name)]
write.table(saf_base, sep = "\t", quote = FALSE, row.names = FALSE)
