#!/usr/bin/env Rscript

# Convert the GTF of all gene models to per-gene SAF files with one entry for
# each base pair.
#
# Usage: Rscript count-bases.R
#
# Input:
#   GTF: data/gencode.v26.GRCh38.genes.gtf
#   Target genes: data/target-genes.txt
#
# Output:
#   Counts: data/counts/<ensembl-gene-id>-<gene-name>-<chromosome>-<start>-<end>.txt

# Setup ------------------------------------------------------------------------

suppressPackageStartupMessages({
library(Rsubread)
})

nthreads <- 1
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  nthreads <- as.integer(args)
}
message("Using ", nthreads, " threads to run featureCounts")
if (nthreads > 28) {
  warning("That's a lot of threads! Did you mean to do that?")
}

bam_files <- list.files(
  path = "/project2/mstephens/dongyue/gtex/",
  pattern = "bam$",
  full.names = TRUE
)
bam_corrupted <- c("/project2/mstephens/dongyue/gtex/SRR1071692.bam",
                   "/project2/mstephens/dongyue/gtex/SRR1077115.bam")
bam_files <- setdiff(bam_files, bam_corrupted)
stopifnot(length(bam_files) > 0)

saf_files <- list.files(
  path = "data/saf/",
  pattern = "saf$",
  full.names = TRUE
)
stopifnot(length(saf_files) > 0)

outdir <- "data/counts/"
dir.create(outdir, showWarnings = FALSE)

# Count and export results -----------------------------------------------------

for (i in seq_along(saf_files)) {
  metadata <- read.table(saf_files[i], header = TRUE, nrow = 1)
  counts <- featureCounts(
    files = bam_files,
    annot.ext = saf_files[i],
    useMetaFeatures = FALSE,
    read2pos = 5,
    isPairedEnd = TRUE,
    requireBothEndsMapped = TRUE,
    countChimericFragments = FALSE,
    nthreads = nthreads,
    tmpDir = tempdir()
  )
  outdata <- cbind(counts$annotation, counts$counts)
  outdata[["Length"]] <- NULL

  outfile <- sprintf("%s-%s-chr%d-%d-%d.txt", metadata$GeneID, metadata$Name,
                     metadata$Chr, metadata$Start, metadata$End)
  outfile <- file.path(outdir, outfile)
  write.table(outdata, file = outfile,
              sep = "\t", quote = FALSE, row.names = FALSE)
}