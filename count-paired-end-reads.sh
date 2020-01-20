#!/bin/bash
set -eu

# Count the paired end reads in BAM file(s)
#
# Usage: bash count-paired-end-reads.sh *.bam

for bamfile in $*
do
  echo $bamfile
  printf "%-20s:\t%10d\n" "Total reads" `samtools view -c $bamfile`
  printf "%-20s:\t%10d\n" "Paired end reads" `samtools view -c -f 1 $bamfile`
  printf "%-20s:\t%10d\n" "Single end reads" `samtools view -c -F 1 $bamfile`
  printf "\n"
done
