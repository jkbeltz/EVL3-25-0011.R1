###Mean baseline > 0.02 and across evolved samples one has MAF > 0.01
setwd('/scratch/groups/dpetrov/MarkB/Orchard2021Data/RData_V2/')
library(tidyverse)


load('./orch2021_Baseline_Downsampled_V2_RAW.RData', verbose = TRUE)
load('./orch2021_Downsampled_META_RAW_V2.RData', verbose = TRUE)


#Combine baseline and time point samples
##Filter so that across all samples, one sample has a HAF > 0.01
df.base = cbind(sites.base, afmat.base)
df.base.eec = cbind(sites.base, eec.base)
df.tps = cbind(sites, afmat)
df.eec = cbind(sites, eec)

df.meta = left_join(df.tps, df.base)
df.meta = na.omit(df.meta)

df.meta = df.meta %>%
    select(-contains("rep")) %>%#exclude biol rep cages
    select(-contains("Rd")) %>% #exclude tech rep cages
    mutate(af.base.mean = (rowMeans(select(., starts_with("Orch")), na.rm = TRUE))) %>%
    filter(af.base.mean > 0.02 & af.base.mean < 0.98) %>%
    filter(if_any(contains("_E"), ~ . > 0.01), if_any(contains("_E"), ~ . < 0.99)) %>%
   filter(if_any(contains("_I"), ~ . > 0.01), if_any(contains("_I"), ~ . < 0.99)) %>%
   filter(if_any(contains("_P"), ~ . > 0.01), if_any(contains("_P"), ~ . < 0.99))



meta.sites = df.meta %>% dplyr::select(chrom, pos) %>%
    mutate(snp = paste0(chrom, "_", pos))


eec.base = df.base.eec %>%
    mutate(snp = paste0(chrom, "_", pos)) %>%
    filter(snp %in% meta.sites$snp)%>%
    select(-c(chrom, pos, snp))

afmat.base =  df.base %>%
    mutate(snp = paste0(chrom, "_", pos))%>%
    filter(snp %in% meta.sites$snp)%>%
    select(-c(chrom, pos, snp))

eec = df.eec %>%
    mutate(snp = paste0(chrom, "_", pos)) %>%
    filter(snp %in% meta.sites$snp)%>%
    select(-c(chrom, pos, snp))

afmat = df.tps %>%
    mutate(snp = paste0(chrom, "_", pos))%>%
    filter(snp %in% meta.sites$snp) %>%
    select(-c(chrom, pos, snp))

sites.base = meta.sites %>% select(-snp)
sites = meta.sites %>% select(-snp)

save(samps.base, sites.base, afmat.base, eec.base, file = './orch2021_Baseline_Downsampled_V2_Filtered_MAF02.RData ')
save(samps, sites, afmat, eec, file = './orch2021_Downsampled_META_Filtered_MAF02_V2.RData')


df = cbind(samps, t(afmat))
df.eec = cbind(samps, t(eec))
df = df %>%
    filter(biol.rep == "No" & tech.rep == "No" & treatment == "E")
samps = df[,1:ncol(samps)]
afmat = df[,-(1:ncol(samps))]
afmat = t(afmat)

df.eec = df.eec %>%
    filter(biol.rep == "No" & tech.rep == "No" & treatment == "E")
eec = df.eec[,-(1:ncol(samps))]
eec = t(eec)

save(sites, samps, afmat, eec, file = "orch2021_Downsampled_ECage_Filtered_MAF02_V2.RData")