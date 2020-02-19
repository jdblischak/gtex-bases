#!/usr/bin/env Rscript

library(data.table)
library(stringr)

counts_summary <- fread("data/standard-counts.txt.summary")
counts_summary
(reads_assigned <- counts_summary[Status == "Assigned", -(Status:Status)])
(reads_total <- colSums(counts_summary[, -(Status:Status)]))
(reads_percent <- reads_assigned / reads_total * 100)

counts <- fread("data/standard-counts.txt", skip = 1)
counts
dim(counts)
length(unique(counts$Geneid))

colnames(counts) <- stringr::str_replace(colnames(counts),
                                         ".+(SRR\\d+\\.bam)",
                                         "\\1")
colnames(counts)
counts$gene <- stringr::str_replace(counts$Geneid,
                                    "(ENSG\\d+)\\.\\d+",
                                    "\\1")
stopifnot(length(unique(counts$Geneid)) == length(unique(counts$gene)))

summary(counts[["SRR1069690.bam"]])
hist(counts[["SRR1069690.bam"]])

gene_totals <- apply(counts[, -(Geneid:Length)][, -(gene:gene)], 1, sum)
names(gene_totals) <- counts$gene
summary(gene_totals)
sum(gene_totals > 0)
mean(gene_totals > 0)
gene_totals_cpm <- log2( (gene_totals + 0.25) / sum(gene_totals))
hist(gene_totals_cpm)
sum(gene_totals_cpm > -20)

target_genes_file <- "data/target-genes.txt"
target_genes <- scan(target_genes_file, what = "character")
gene_totals_cpm[target_genes]
sum(gene_totals_cpm[target_genes] > -20)
sort(gene_totals_cpm[target_genes], decreasing = TRUE)[1:12]
gene_totals_cpm["ENSG00000231500"]
gene_totals_cpm["ENSG00000104904"]
