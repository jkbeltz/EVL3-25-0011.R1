#!/bin/bash
#SBATCH --job-name=02G_MeanDepth
#SBATCH --time=23:00:00
#SBATCH --mem=400GB
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem
#SBATCH --out=out.02G_MeanDepth

#cd /scratch/groups/dpetrov/MarkB/Orchard2021Data/03_bams/
Dir=${1}

cd ${Dir}

for i in ./*.dp.txt;do
    awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' $i > ${i/.dp.txt/_MeanDP.txt}
done

