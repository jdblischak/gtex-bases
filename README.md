# Per-base counts for GTEx samples

The goal is to calculate counts for each base of a gene for a subset of the GTEx
samples. These will be the input for a statistical model that will cluster the
samples based on topic modelling of the splicing patterns.

**Notes:**
* The reads are paired-end
* Bases in introns are included
* Each read is reduced to its 5' base. In other words, each read only
  contributes a count to a single base of the gene (Subread option `--read2pos
  5`)

## Data

* The Ensembl gene IDs of the target genes are in [data/target-genes]()
* The GTEx BAM files are currently located on Midway2 at
  `/project2/mstephens/dongyue/gtex/`

BAM files that are currently corrupted:
* SRR1071692.bam
* SRR1077115.bam

## Installation

Install conda via Miniconda and then run the following:

```
conda env create --file environment.yml
conda activate gtex-bases
```

## Instructions

1. Download the GTF file with the gene models used by GTEx V8
    ```
    Rscript setup.R
    ```
    This creates the file `data/gencode.v26.GRCh38.genes.gtf`

1. Convert the GTF file to per-gene SAF files with one entry per base pair
    ```
    Rscript gtf2saf.R
    ```
    For each target gene, this creates a file `data/saf/<ensembl-gene-id>.saf`

1. Count the number of reads per base with featureCounts
    ```
    bash count-bases.sh
    ```

## Miscellaneous

* Count the number of paired and single end reads
    ```
    bash count-paired-end-reads.sh /project2/mstephens/dongyue/gtex/*.bam
    ```
