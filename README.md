# EVL3-25-0011.R1
Code which accompanies for manuscript published in Evolution Letters (accepted 8/25)

Variation in the resource environment affects patterns of seasonal adaptation at phenotypic and genomic levels in Drosophila melanogaster, Jack K. Beltz1, Mark Christopher Bitter2, August Goldfischer1, Dmitri A. Petrov2, Paul Schmidt1, currently in review at Evolution Letters

Phenotypic Analysis - Raw data files found in raw pheno folder, and analysis scripts found in .R file. Data files contained in "Raw", code is designed to be run O20_AB_fig1.R, then O20_AB_fig2.R, followed by the Haldanes, and PCAPheno files. The final version of Figure 2 is generated in O20-Haldanes, and figure 3 in ClusterHeatMap.

16s Analysis - Microbial analysis conducted using R and QIIME2 pipeline, following (https://github.com/jbisanz/16Spipelines.git), 16s Fasta files found on dryad (LINK).

WGS analysis -All samples sequenced on Illumina Novaseq 6000 using 150 bp, paired-end reads -Samples with coverage < 5x were re-sequenced on a second round/date

RD1Sequencing_SampleKey.csv and RD2Sequencing_SampleKey.csv contain information regarding .fastq file names
-BioinformaticsScripts: contains scripts used for read trimming, quality control, alignment, and haplotype-derived allele frequency estimation.

-RData: consists of 4 data frames: -samps: sample information whereby row order corresponds to column order of the allele frequency matrix and effective coverage matrix. -afmat: numeric data frame containing all haplotype-derived allele freqeuncies - eec: numeric data frame containing estimated effective coverage for each sample/site -sites: dataframe containing chromosome and site information (corresponding to rows of afmat and eec_