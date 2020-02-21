#!/usr/bin/env Rscript

# Create per-gene SAF files with one entry for each base pair.
#
# Usage: Rscript create-saf-files.R
#
# Input:
#   Ensembl exons: data/exons.txt
#   Target genes: data/target-genes.txt
#
# Output:
#   SAF: data/saf/<ensembl-gene-id>.saf

# Setup ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

exons_file <- "data/exons.txt"
target_genes_file <- "data/target-genes.txt"
stopifnot(file.exists(exons_file), file.exists(target_genes_file))

dir.create("data/saf/", showWarnings = FALSE)

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


# Remove "chr" from chromosome
saf[, Chr := substr(Chr, 4, nchar(Chr))]

# Confirm that each gene is only present on one chromosome
stopifnot(
  identical(
    length(unique(paste(saf$GeneID, saf$Chr))),
    length(unique(saf$GeneID))
  )
)

stopifnot(saf$Start <= saf$End)

# Reduce to one entry per gene
setkey(saf, GeneID)
saf_gene <- saf[, list(Start = min(Start), End = max(End)),
                by = list(GeneID, Chr, Strand, Name)]
setkey(saf_gene, GeneID)
stopifnot(nrow(saf_gene) == length(unique(saf_gene$GeneID)))

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

  stopifnot(saf_base$End == saf_base$Start)

  saf_base <- saf_base[, list(GeneID, Chr, Start, End, Strand, Name)]
  saf_outfile <- paste0("data/saf/", gene_id, ".saf")
  write.table(saf_base, file = saf_outfile,
              sep = "\t", quote = FALSE, row.names = FALSE)
}
