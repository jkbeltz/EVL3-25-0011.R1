#!/bin/bash
#SBATCH --job-name=02F_AlignStats
#SBATCH --time=08:00:00
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=4
#SBATCH --partition=hns
#SBATCH --out=out.02F.AlignStats.sbatch

module load biology
module load bamtools

#cd /scratch/groups/dpetrov/MarkB/Orchard2021Data/03_bams/
Dir=${1}

cd ${Dir}

for i in ./*.bam;do
    bamtools stats -in $i > ${i/.bam/.AlignStats.txt}
done
