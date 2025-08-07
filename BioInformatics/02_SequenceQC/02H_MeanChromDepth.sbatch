#!/bin/bash
#SBATCH --job-name=02H_MeanChromDepth
#SBATCH --time=23:00:00
#SBATCH --mem=300GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem
#SBATCH --out=out.02H_MeanChromDepth

#cd /scratch/groups/dpetrov/MarkB/Orchard2021Data/03_bams/
Dir=${1}

cd ${Dir}

ml R

Rscript /home/users/mcbitter/OrchardProject/Code/Orchard2021_Scripts/BioInformatics/02_SequenceQC/MeanChromDP.R  
