#!/bin/bash
set -eu

# Check the chromosome names used in BAM file(s)
#
# Usage: bash check-bam-chromosomes.sh *.bam

for bamfile in $*
do
  echo $bamfile
  printf "\n"
  samtools view $bamfile | cut -f3 | sort | uniq
  printf "___\n\n"
done
