#!/bin/bash
#SBATCH --job-name=filter.rdata.downsampled
#SBATCH --time=04:00:00
#SBATCH --mem=300GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem
#SBATCH --out=out.filter.rdata.downsampled

ml R

Script=${1}

Rscript ${Script}
