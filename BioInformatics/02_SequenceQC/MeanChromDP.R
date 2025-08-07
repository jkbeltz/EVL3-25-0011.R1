library(tidyverse)
library(data.table)
setwd('/home/users/mcbitter/dpetrov/MarkB/Orchard2020Data/03_bams/DownsampledBams')
files = list.files(path = "./", pattern = '.dp.txt')

for (file in files){
    df = fread(file)
    samp = strsplit(file, ".dp.txt")[[1]][1]
    colnames(df) = c('chrom', 'pos', 'dp')
    df = df %>%
        group_by(chrom) %>%
        summarise(mean_dp = mean(dp))
    df = df %>%
        filter(chrom %in% c("2L", "2R", "3L", "3R", "X"))
    write.csv(df, paste0(samp, ".MeanChromDP.txt"), row.names = FALSE)
}