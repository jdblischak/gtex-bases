# Per-base counts for GTEx samples

The goal is to calculate counts for each base of a gene for a subset of the GTEx
samples. These will be the input for a statistical model that will cluster the
samples based on topic modelling of the splicing patterns.

**Notes:**
* The reads are paired-end
* Bases in introns are included
* Each read is reduced to its 5' base. In other words, each read only
  contributes a count to a single base of the gene (featureCounts option
  `--read2pos 5`)
* The GTEx BAM files use Ensembl chromosome names, e.g. 1, 2, MT, etc. Thus
`data/exons.txt` and the SAF files in `data/saf/` use Ensembl chromosome names.
However, [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html)
uses UCSC chromosome names, e.g. chr1, chr2, chrM, etc. Thus the chromosome
names are converted just prior to plotting.

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

1. Download exons from [Ensembl 98](http://sep2019.archive.ensembl.org/) (Sep
2019, [Gencode 32](https://www.gencodegenes.org/human/release_32.html),
GRCh38.p13).
    ```
    Rscript download-exons.R
    ```
    This creates the file `data/exons.txt`

1. Create per-gene SAF files with one entry per base pair
    ```
    Rscript create-saf-files.R
    ```
    For each target gene, this creates a file `data/saf/<ensembl-gene-id>.saf`

1. Count the number of reads per base with featureCounts
    ```
    Rscript count-bases.R
    ```
    Alternatively submit to the cluster, passing the optional argument for the
    number of threads to use for featureCounts:
    ```
    sbatch --partition=broadwl -J count-bases \
      -o count-bases-stdout.txt -e count-bases-stderr.txt \
      --mem=32G --nodes=1 --tasks-per-node=8 \
      count-bases.R 8
    ```

## Miscellaneous

* Count the number of paired and single end reads
    ```
    bash count-paired-end-reads.sh /project2/mstephens/dongyue/gtex/*.bam
    ```

* Run a standard gene-level featureCounts run (using GTEx exons)
    ```
    Rscript download-exons-gtex.R
    bash run-standard-featurecounts.sh /project2/mstephens/dongyue/gtex/*.bam
    ```
