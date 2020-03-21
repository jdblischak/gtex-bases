#!/usr/bin/env Rscript

library(data.table)

dir_counts <- "data/counts/"

files <- list.files(path = dir_counts, pattern = "txt$", full.names = TRUE)

greater_than_0 <- 0
greater_than_1000 <- 0
for (f in files) {
  print(f)
  counts <- fread(f)
  counts <- counts[, -(GeneID:Length)]
  total <- sum(colSums(counts))
  print(total)
  if (total > 0) greater_than_0 <- greater_than_0 + 1
  if (total > 1000) greater_than_1000 <- greater_than_1000 + 1
}
cat(sprintf("Genes with more than 0 counts: %d\n", greater_than_0))
cat(sprintf("Genes with more than 1000 counts: %d\n", greater_than_1000))
