#!/bin/bash                                                                                                                                                                                                                           
#SBATCH --job-name=02.sequenceQC.orch2021                                                                                                                                                                                                
#SBATCH --time=08:00:00                                                                                                                                                                                                               
#SBATCH --mem=50                                                                                                                                                                                                                      
#SBATCH --cpus-per-task=4                                                                                                                                                                                                             
#SBATCH --partition=hns                                                                                                                                                                                                               
#SBATCH --out=02.sequenceQC.orch2021.out 

module load biology
module load samtools
module load bamtools
module load java


cd /scratch/groups/dpetrov/MarkB/Orchard2021Data/02_fastqs/raw/
scp ./*.sam ../../03_bams/


cd /scratch/groups/dpetrov/MarkB/Orchard2021Data/03_bams/
for i in ./*.sam;do
    samtools view -S -b $i > ${i/.sam/.bam}
done



for i in ./*.bam;do
    samtools view -bq 20 $i > ${i/.bam/.filtered.bam}
done

for i in ./*.filtered.bam;do
    samtools sort $i > ${i/.filtered.bam/.filtered.sorted.bam}
done

for i in ./*.filtered.sorted.bam;do
    samtools index -b $i 
done

for i in ./*.filtered.sorted.bam;do
    java -jar /home/users/mcbitter/dpetrov/SOFTWARE/picard.jar MarkDuplicates \
        I=./$i \
        O=./${i/.bam/.deduped.bam} \
        M=/scratch/groups/dpetrov/MarkB/SequencingPractice/${i/.bam/marked_dup_metrics.txt}
        
done

rm *.filtered.bam
rm *.filtered.sorted.bam

for i in ./*.sorted.deduped.bam;do
    samtools depth -a $i > ${i/.sorted.deduped.bam/.dp.txt}
done

for i in ./*.sorted.deduped.bam;do
    bamtools stats -in $i > ${i/.sorted.deduped.bam/.AlignStats.txt}
done

for i in ./*.dp.txt;do
    awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' $i > ${i/.dp.txt/mean.dp.txt}
done

mv *dp* ./DP_AlignStats
mv *AlignStats.txt ./DP_AlignQual
