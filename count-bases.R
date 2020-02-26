#!/usr/bin/env Rscript

# Count the number of reads mapping to each base pair of a given gene.
#
# Usage: Rscript count-bases.R
#
# Input:
#   BAM: /project2/mstephens/dongyue/gtex/SRRXXXXXXX.bam
#   SAF: data/saf/<ensembl-gene-id>.saf
#
# Output:
#   Counts: data/counts/<ensembl-gene-id>-<gene-name>-<chromosome>-<start>-<end>.txt

# Setup ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Rsubread)
  library(stringr)
})

nthreads <- 1

message("Using ", nthreads, " threads to run featureCounts")
if (nthreads > 28) {
  warning("That's a lot of threads! Did you mean to do that?")
}

bam_files <- list.files(
  path = "/project2/mstephens/dongyue/gtex",
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
  message("Processing ", saf_files[i])
  saf <- read.table(saf_files[i], header = TRUE, stringsAsFactors = FALSE)
  saf$Chr <- as.character(saf$Chr)
  counts <- featureCounts(
    files = bam_files,
    annot.ext = saf_files[i],
    useMetaFeatures = FALSE,
    read2pos = 5,
    countMultiMappingReads = FALSE,
    isPairedEnd = TRUE,
    requireBothEndsMapped = TRUE,
    countChimericFragments = FALSE,
    nthreads = nthreads,
    tmpDir = tempdir(),
    verbose = TRUE
  )
  colnames(counts$counts) <- stringr::str_extract(colnames(counts$counts),
                                                  "SRR\\d+\\.bam")
  colnames(counts$stat) <- stringr::str_replace(colnames(counts$stat),
                                                ".+(SRR\\d+\\.bam)",
                                                "\\1")
  # Avoid duplicate rowname warning:
  rownames(counts$counts) <- seq_len(nrow(counts$counts))
  outdata <- cbind(counts$annotation, counts$counts)
  outdata[["Length"]] <- NULL

  outfile <- sprintf("%s-%s-chr%s-%d-%d.txt", saf$GeneID[1], saf$Name[1],
                     saf$Chr[1], min(saf$Start), max(saf$End))
  outfile <- file.path(outdir, outfile)
  write.table(outdata, file = outfile,
              sep = "\t", quote = FALSE, row.names = FALSE)
  # Save the counting summary statistics
  outfile_sum <- paste0(outfile, ".summary")
  write.table(counts$stat, file = outfile_sum,
              sep = "\t", quote = FALSE, row.names = FALSE)
}
