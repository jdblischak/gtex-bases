#!/bin/bash
set -eu

# Run a standard gene-level featureCounts run
#
# sbatch --partition=mstephens --account=pi-mstephens -J featurecounts \
#   -o featurecounts-stdout.txt -e featurecounts-stderr.txt \
#   --mem=32G --nodes=1 --tasks-per-node=8 \
#   run-standard-featurecounts.sh /project2/mstephens/dongyue/gtex/*bam

gtf="data/gencode.v26.GRCh38.genes.gtf"
output="data/standard-counts.txt"

featureCounts -a $gtf -o $output -p -B -C --tmpDir /tmp -T 8 $*
