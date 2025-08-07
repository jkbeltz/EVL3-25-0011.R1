#!/bin/bash
#SBATCH --job-name=01.orch.2020.alignment
#SBATCH --time=48:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=6
#SBATCH --partition=hns
#SBATCH --out=out.01.orch.2020.alignment

module load biology
module load bwa

Dir=${1}

cd ${Dir}

bwa index ../../../Dmel_RefGenome/all_dmel.fasta


awk 'NR>0' orch2020.samps.csv > orch2020.samps.txt

while read samp
do
    bwa mem -t 12 ../../../Dmel_RefGenome/all_dmel.fasta "$samp"_R1_001.fastq.gz "$samp"_R2_001.fastq.gz > "$samp".sam
done < orch2020.samps.txt
