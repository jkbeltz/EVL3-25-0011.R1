##Generate EEC and add to RData
##Eec for time point samples   
source("/home/users/mcbitter/OrchardProject/Code/workflow_functions.R")
source("/home/users/mcbitter/OrchardProject/Code/general_cage_functions.R")
setwd('/scratch/groups/dpetrov/MarkB/Orchard2021Data/PigTraitMapping/AlignStats_Coverage/')

load("../RData/orch2021.PigTraitMap.RData")

pct.missing = fread('../../founders/inbredv2.filtered.Orch21.PctMissing.csv')
gens =  as.data.frame(cbind(c(2, 8), c(6, 10))) 
names(gens) = c("tpt", "gens")
gens = gens %>%
    mutate(tpt = as.character(tpt))
chroms = c("2L", "2R", "3L", "3R", "X")
n.snps = c(nrow(filter(sites, chrom == "2L")), nrow(filter(sites, chrom == "2R")), nrow(filter(sites, chrom == "3L")), nrow(filter(sites, chrom == "3R")), nrow(filter(sites, chrom == "X")))
chrom_length <- 23000000
recomb_rate <- 0.0000000239

files <- list.files(pattern="*MeanChromDP*", full.names = TRUE, recursive=FALSE)

##EEC matrix generation
#files <- list.files(pattern="*MeanChromDP*", full.names = TRUE, recursive=FALSE)
data <- data.frame()
for (file in files){
    clean.file <- strsplit(file, "./")[[1]][2]
    full.sample.name <- strsplit(clean.file, ".MeanChromDP.txt")[[1]][1]
    chrom.dp = read.csv(file)
    chrom.dp$sample = full.sample.name
    chrom.dp = chrom.dp %>%
        dplyr::select(sample, chrom, mean_dp)
    data = rbind(data, chrom.dp)
}

samps.ext = as.data.frame(samps)
samps.ext$n.times = 5
samps.ext = as.data.frame(lapply(samps.ext, rep, samps.ext$n.times))
samps.ext = as.data.frame(samps.ext[,-5])
samps.ext = samps.ext %>% rename(sample = sample.name)

samps.ext = samps.ext %>%
    mutate(sample = as.character(sample))

data = data[order(match(data[,1],samps.ext[,1])),]

data = data %>% rowwise() %>% mutate(tpt = strsplit(sample, '[.]')[[1]][1]) %>% separate(tpt, into = c(NA, 'tpt'), 'TP')
data = left_join(data, gens)
data = data %>%
    mutate(num.snps = rep(as.numeric(as.character(n.snps)), length(unique(data$sample))),
           pct.missing = rep(pct.missing$pct, length(unique(data$sample))),
            chrom.length = chrom_length,
           recomb.rate = recomb_rate)

data = data %>% rowwise() %>%
    mutate(eec = calc_expected_ec(mean_dp, gens, pct.missing, num.snps, chrom.length, recomb.rate))

sp.data = data %>%
    dplyr::select(sample, chrom, eec) %>%
    spread(sample, eec)
sp.data = sp.data %>%
    mutate(n.times = as.numeric(as.character(n.snps)))
sp.data = as.data.frame(lapply(sp.data, rep, sp.data$n.times))
eec.sites = sp.data %>%
    dplyr::select(chrom)
sp.data = sp.data %>%
    dplyr::select(-c(chrom, n.times))
eec.samps = as.data.frame(names(sp.data))
eec = sp.data
save(sites, samps, eec, afmat, file = "../RData/orch2021.PigTraitMap_EECIncluded.RData")