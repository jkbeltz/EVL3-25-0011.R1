#!/bin/bash
#SBATCH --job-name=ComputeEEC
#SBATCH --time=03:00:00
#SBATCH --mem=300GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem
#SBATCH --out=out.ComputeEEC

ml R

Rscript ComputeEEC.R

