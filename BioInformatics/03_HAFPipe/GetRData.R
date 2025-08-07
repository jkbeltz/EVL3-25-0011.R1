##Get unfiltered RData -- Run on Cluster

library(tidyverse)
library(data.table)
setwd('~/dpetrov/MarkB/Orchard2020Data/04_HAFs/')

tpts = c(0, 1, 2, 3, 4, 5)


for (tp in tpts){
    data = data.frame()
    files = list.files(pattern=paste0("TP",tp,"_"), full.names=TRUE, recursive=FALSE)
    data = data.frame()
    for (file in files){
    clean.file <- strsplit(file, "./")[[1]][2]
    full.sample.name <- strsplit(clean.file, "_downsampled")[[1]][1]
    chr <- strsplit(clean.file, "[.]")[[1]][3] #extract chrom info
    freqs <- fread(file, header = TRUE)
    freqs = freqs %>%
        mutate(sample.name = full.sample.name) %>%
        mutate(chrom = chr)
    freqs = freqs %>%
        dplyr::select(sample.name, chrom, pos, af)
    data = rbind(data, freqs)
    }
    sp.data = data %>%
        spread(sample.name, af)
    sp.data = na.omit(sp.data)

    sites = as.data.frame(sp.data %>%
        dplyr::select(chrom, pos))
    colnames(sites) = c('chrom', 'pos')

    samps = as.data.frame(colnames(sp.data)[3:ncol(sp.data)])
    colnames(samps) = c('sample.name')

    afmat = as.matrix(sp.data[,3:ncol(sp.data)])
    samps = samps %>% rowwise() %>%
        rename(full.sample.name = sample.name) %>%
        mutate(tpt = strsplit(strsplit(as.character(full.sample.name), '[_]')[[1]][1], "TP")[[1]][2]) %>%
        mutate(treatment = strsplit(as.character(full.sample.name), '[_]')[[1]][2]) %>%
        mutate(cage = strsplit(as.character(full.sample.name), '[_]')[[1]][3])
    save(sites, samps, afmat, file = paste0("../RData/orch2020.",tp , ".RData"))
}


setwd('~/dpetrov/MarkB/Orchard2020Data/RData/')
files <- c('orch2020.0.RData', 'orch2020.1.RData', 'orch2020.2.RData', 'orch2020.3.RData', 'orch2020.4.RData', 'orch2020.5.RData')


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

save(afmat, sites, samps, file = "./orch2020_Unfiltered_META.RData")

