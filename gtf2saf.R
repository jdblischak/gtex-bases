#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(data.table)
library(GenomicRanges)
library(rtracklayer)
})

gtf <- "data/gencode.v26.GRCh38.genes.gtf"

gr <- import(gtf, format = "gtf")

# SAF columns:
# GeneID  Chr     Start   End     Strand  Name

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

# Subset to ENSG00000223745
saf_sub <- saf_gene["ENSG00000223745"]

# Expand to one entry per base pair
saf_base <- saf_sub[, list(Start = seq(Start, End), End = seq(Start, End) + 1),
                     by = list(GeneID, Chr, Strand, Name)]

str(saf_base)
stopifnot(saf_base$End - saf_base$Start == 1)

saf_base <- saf_base[, list(GeneID, Chr, Start, End, Strand, Name)]
write.table(saf_base, file = "data/ENSG00000223745.saf",
            sep = "\t", quote = FALSE, row.names = FALSE)
