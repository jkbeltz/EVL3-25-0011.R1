##joining TP .RData
library(tidyverse)
library(data.table)

setwd("/home/users/mcbitter/dpetrov/MarkB/Orchard2021Data/PigTraitMapping/RData/")
files <- c('orch2021.PigTraitMap.8.RData', 'orch2021.PigTraitMap.2.RData')


chrom = c("2L", "2R", "3L", "3R", "X")
pos = NA
af.meta = as.data.frame(cbind(chrom, pos))
af.meta$pos = as.integer(af.meta$pos)
samps.meta = data.frame()

for (file in files){
    load(file)
    df.samps = samps
    df.af = cbind(sites, afmat)
    af.meta = right_join(af.meta, df.af)
    samps.meta = rbind(samps.meta, df.samps)
}


afmat = na.omit(af.meta)
sites = afmat[,c(1:2)]
afmat = afmat[,-c(1:2)]
samps = samps.meta

save(afmat, sites, samps, file = "./orch2021.PigTraitMap.RData")