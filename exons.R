#!/usr/bin/env Rscript

library(biomaRt)

# GTEx used Gencode release 26 for gene annotations. This corresponds to genome
# build GRCh38.p10. The last updated version of GRCh38.p10 for Ensembl was
# Ensebml 91 from December 2017.
#
# https://www.gencodegenes.org/human/releases.html
# https://www.gencodegenes.org/human/release_26.html
# http://dec2017.archive.ensembl.org/Homo_sapiens/

ensembl <- useMart(host = "dec2017.archive.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")
exons_all <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id",
                                  "chromosome_name", "exon_chrom_start",
                                  "exon_chrom_end", "strand",
                                  "external_gene_name",
                                  "gene_biotype", "transcript_biotype"),
                   filters = "ensembl_gene_id",
                   values = "ENSG00000104904",
                   mart = ensembl)
# exons_final <- exons_all[exons_all$chromosome_name %in% c(1:22, "X", "Y", "MT") &
#                            exons_all$gene_biotype == "protein_coding",
#                          c("ensembl_gene_id", "chromosome_name", "exon_chrom_start",
#                            "exon_chrom_end", "strand", "external_gene_name")]

write.table(exons_all, "data/exons-ENSG00000104904.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
