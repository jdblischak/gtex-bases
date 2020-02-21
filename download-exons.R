#!/usr/bin/env Rscript

# Download exons from Ensembl.
#
# Ensembl 98 (Sep 2019) - http://sep2019.archive.ensembl.org/
# Gencode 32 (Sep 2019) - https://www.gencodegenes.org/human/release_32.html
# GRCh38.p13

library(biomaRt)

archive <- "sep2019.archive.ensembl.org"

ensembl <- useMart(
  host = archive,
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)

exons_all <- getBM(
  attributes = c(
    "ensembl_exon_id",
    "ensembl_transcript_id",
    "ensembl_gene_id",
    "chromosome_name",
    "exon_chrom_start",
    "exon_chrom_end",
    "strand",
    "external_gene_name",
    "gene_biotype",
    "transcript_biotype"
  ),
  mart = ensembl
)

exons_final <- exons_all[exons_all$chromosome_name %in% c(1:22, "X", "Y", "MT"), ]

exons_final <- exons_final[order(exons_final$chromosome_name,
                                 exons_final$exon_chrom_start,
                                 exons_final$exon_chrom_end), ]

exons_final$chromosome_name <- paste0("chr", exons_final$chromosome_name)
exons_final$chromosome_name <- sub("chrMT", "chrM", exons_final$chromosome_name)

exons_final$strand <- ifelse(exons_final$strand == 1, "+", "-")

stopifnot(
  identical(
    length(unique(paste(exons_final$ensembl_transcript_id,
                        exons_final$ensembl_exon_id))),
    nrow(exons_final)
  )
)

write.table(exons_final, "data/exons.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
