#!/bin/bash

# sbatch submission script to run main snakemake process. It then submits
# individual jobs from the compute node.

#SBATCH --job-name=snakemake
#SBATCH --output=snakelog.out
#SBATCH --time=24:00:00
#SBATCH --partition=mstephens
#SBATCH --account=pi-mstephens
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4

snakemake \
    -kp \
    --ri \
    -j 100 \
    --max-jobs-per-second 5 \
    --use-conda \
    --cluster-config cluster.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition=broadwl \
        --account=pi-mstephens \
        --job-name={cluster.name} \
	    --output={cluster.logfile}" \
    $*
