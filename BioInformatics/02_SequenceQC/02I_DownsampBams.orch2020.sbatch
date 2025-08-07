#!/bin/bash
#SBATCH --job-name=DownSampBams
#SBATCH --time=02:00:00
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=hns
#SBATCH --out=out.DownSampBams

ml biology
ml samtools

cd /home/users/mcbitter/dpetrov/MarkB/Orchard2020Data/03_bams

File=${1}
AmountDown=${2}


samtools view -bs 42.${AmountDown} ${File} > ./DownsampledBams/${File/.bam/_downsamped.bam}
