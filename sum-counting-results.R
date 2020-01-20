#!/usr/bin/env Rscript

library(data.table)

dir_counts <- "data/counts/"

files <- list.files(path = dir_counts, pattern = "txt$", full.names = TRUE)

for (f in files) {
  print(f)
  counts <- fread(f)
  counts <- counts[, -(GeneID:Strand)]
  total <- sum(colSums(counts))
  print(total)
}
