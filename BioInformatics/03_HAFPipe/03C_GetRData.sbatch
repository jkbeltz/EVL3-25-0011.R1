#!/bin/bash
#SBATCH --job-name=Get_tp_RData
#SBATCH --time=08:00:00
#SBATCH --mem=500GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem
#SBATCH --out=out.Get_tp_RData

ml R

Rscript GetRData.R
