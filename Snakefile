import os

configfile: "config.yml"


dir_data = config["dir_data"]
dir_external = config["dir_external"]
dir_scratch = config["dir_scratch"]
dir_bam = dir_external + "bam/"
dir_bam_pe = dir_external + "bam-pe/"
dir_exons = dir_data + "exons/"
dir_saf = dir_data + "saf/"
dir_counts = dir_data + "counts/"
dir_counts_gene = dir_data + "counts-gene/"

assert os.path.exists(dir_data), "Local data directory exists"
assert os.path.exists(dir_external), "External data directory exists"

# Directory to send log files. Needs to be created manually since it
# is not a file created by a Snakemake rule.
dir_log = config["dir_log"]
if not os.path.isdir(dir_log):
    os.mkdir(dir_log)

samples = config["samples"]
genes = config["genes"]

localrules: download_exons

rule all:
    input: expand(dir_counts + "{gene}.txt", gene = genes),
           expand(dir_counts_gene + "{gene}.txt", gene = genes)

rule target_genes:
    input: expand(dir_saf + "{gene}.saf", gene = genes)

rule target_samples:
    input: expand(dir_bam_pe + "{sample}.bam", sample = samples)

# The experiment was paired-end, but the BAM files contain some single end
# reads. featureCounts especially does not like this, so remove them.
rule remove_se_reads:
    input: bam = dir_bam + "{sample}.bam"
    output: dir_bam_pe + "{sample}.bam"
    conda: "envs/samtools.yml"
    shell: "samtools view -b -f 1 {input} > {output}"

rule download_exons:
    output: dir_exons + "{gene}.txt"
    params: gene = "{gene}"
    conda: "envs/biomart.yml"
    shell: "Rscript download-exons.R {params.gene} > {output}"

rule create_saf:
    input: dir_exons + "{gene}.txt"
    output: dir_saf + "{gene}.saf"
    params: gene = "{gene}"
    conda: "envs/dt.yml"
    shell: "Rscript create-saf-file.R {input} > {output}"

rule count_bases:
    input: saf = dir_saf + "{gene}.saf",
           bam = expand(dir_bam_pe + "{sample}.bam", sample = samples)
    output: dir_counts + "{gene}.txt"
    threads: 8
    conda: "envs/subread.yml"
    shell:
      "featureCounts -a {input.saf} -o {output} -F SAF -f --read2pos 5 -p -B -C -T {threads} --tmpDir {dir_scratch} -O -R CORE --Rpath {dir_external} {input.bam} ;\n"
      # Format the output file
      "sed -i '1d' {output} ;\n"
      "sed -i 's/Geneid/GeneID/' {output} ;\n"

# Instead of counting per base, count for the entire gene (i.e. meta-feature)
rule count_gene:
    input: saf = dir_saf + "{gene}.saf",
           bam = expand(dir_bam_pe + "{sample}.bam", sample = samples)
    output: dir_counts_gene + "{gene}.txt"
    threads: 8
    conda: "envs/subread.yml"
    shell:
      "featureCounts -a {input.saf} -o {output} -F SAF --read2pos 5 -p -B -C -T {threads} --tmpDir {dir_scratch} {input.bam} ;\n"
      # Format the output file
      "sed -i '1d' {output} ;\n"
      "sed -i 's/Geneid/GeneID/' {output} ;\n"
