setwd("/Users/jackbeltz/Documents/PENN/Dissertation/CH 2:3 (Orchard 20) /AB/RAW/raw geno")
library('tidyverse')
library('dplyr')
library('ggplot2')
library(tibble)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(ggpubr)
library(cowplot)
library(grid)


####PCAS ####
pca_data=read.csv("PCAResults_AppleVBloom.csv")
View(pca_data)


p1 = ggplot(pca_data, aes(x = PC1, y = PC2, color = Timepoint, shape = treatment)) +
  geom_point(size = 3) +
  scale_colour_gradient2(low = "red", mid= "lightblue",high = "blue", midpoint = 3 , limits = c(1, 5), ) +
  scale_shape_manual(name="",breaks =c("A", "B"), values=c(17,19),labels=c("LQ","HQ"))+
  theme_minimal(base_size=12)


p2 = ggplot(pca_data, aes(x = PC1, y = PC3, color = treatment, shape = treatment )) +
  geom_point(size = 3) +
  scale_shape_manual(name="",breaks =c("A", "B"), values=c(17,19),labels=c("Low Quality","High Quality"))+
  #scale_colour_gradient2(low = "red", mid= "white",high = "blue", midpoint = 3 , limits = c(1, 5), ) +
  scale_color_manual(name="",breaks =c("A", "B"),values = c("grey","black"), labels=c("LQ","HQ"))+
  theme_minimal(base_size=12)


pAB<-ggarrange(p1, p2, ncol=2, widths = c(1,1), labels = c("A","B"), common.legend = TRUE, legend=c('right') )
pAB

####PCA STATS####
#Mark Stats
treats = c('A', 'B')
stats = data.frame()
for (treat in treats){
  d = pca_data %>% filter(treatment == treat)
  pval = summary(aov(PC1 ~ Timepoint, data = d))[[1]][1, 5]
  d.s = cbind(treat, pval)
  stats = rbind(stats, d.s)
  
}
stats

geno.pca.treat.manov<-manova(cbind(PC1,PC3)~Timepoint*treatment,pca_data)
summary(geno.pca.treat.manov)

geno.pca.time.manov<-manova(cbind(PC1,PC2)~Timepoint,pca_data)
summary(geno.pca.time.manov)

geno.pca.interaction.manov<-manova(cbind(PC1,PC2)~Timepoint*treatment,pca_data)
summary(geno.pca.interaction.manov)


####GLM####
##GLM of apple and bloom t1->5 (run on cluster): 

library('dplyr')
library('tidyr')
library('data.table')
install.packages("doMC")
library('doMC')

###FUNCTIONS
fit_GLM_ContinuousTime = function (afMatrix, rdMatrix, sampleData, vec, model.vars, poolCt = 100, 
                                   ncores) 
{
  registerDoMC(ncores)
  df = as.data.frame(sampleData[, colnames(sampleData) %in% 
                                  model.vars])
  colnames(df) <- model.vars
  formulaString = paste0(colnames(df), collapse = " + ")
  formulaString = paste0("cts ~ ", formulaString)
  cat("Model Formula is: \n", formulaString, "\n")
  Neff = calc_Neff(rdMatrix, poolCt)
  do.call(rbind, mclapply(1:nrow(afMatrix), function(ix) {
    if (ix%%10000 == 0) {
      cat("working on site ", ix, "\n")
    }
    cts = cbind(t(round(Neff[ix, ] * afMatrix[ix, ])), t(round(Neff[ix, 
    ] * (1 - afMatrix[ix, ]))))
    df$cts = cts
    df = as.data.frame(df)
    model = glm(formulaString, family = "quasibinomial", 
                data = df)
    cp = summary(model)$coefficients[-1, c(1, 4), drop = FALSE]
    results = c(cp[, 1], cp[, 2])
    names(results) = c(paste0("coef.", vec[1], "_", vec[length(vec)]), 
                       paste0("p.", vec[1], "_", vec[length(vec)]))
    return(results)
  }, mc.cores = ncores))
}

#calucluate N effective
calc_Neff = function (rd, poolCt) 
{
  ((poolCt * 2 * rd) - 1)/(poolCt * 2 + rd)
}

###Loop through treatments and run glm - save results as RData file
treats = c('A', 'B')
load('orch2020_filtered.RData')
head(orch2020_filtered)
df.glm.meta = sites
vec = c(1:5)
for (treat in treats){
  load('orch2020_filtered.RData')
  samps = samps %>% mutate(tpt = as.numeric(tpt))
  df.glm = sites
  afmat = t(afmat)
  eec = t(eec)
  afmat = cbind(samps, afmat)
  eec = cbind(samps, eec)
  afmat = afmat %>%
    filter(tpt %in% vec & treatment == treat)
  eec = eec %>%
    filter(tpt %in% vec & treatment == treat)
  samps = afmat[,1:ncol(samps)]
  afmat = afmat[,-c(1:ncol(samps))]
  afmat = t(afmat)
  afmat = as.data.frame(afmat)
  eec = eec[,-c(1:ncol(samps))]
  eec = t(eec)
  eec = as.data.frame(eec)
  res = fit_GLM_ContinuousTime(afMatrix = afmat ,rdMatrix = eec, vec = vec, sampleData = samps, model.vars = 'tpt', poolCt=100, ncores = 16)
  df.glm = cbind(df.glm, res) %>% dplyr::select(-chrom, -pos)
  names(df.glm) = c(paste0('coef.1_5'),paste0('p.1_5'))
  save(df.glm, file = paste0("../05_GLM/glm.", treat, ".TP1_5.RData"))
}
View(df.glm)

##Getting significant sites of each treatment through time:
get.sig.sites = function(glm.file, rdata, comps, fdrThreshs, esThreshs){
  chroms = c('2L', '2R', '3L', '3R', 'X')
  load(rdata)
  af.shifts = get_af_shifts(afmat, samps, cage_set = NULL, comparisons = comps)
  load(glm.file)
  FDR = get_glm_FDR.V2(df.glm = df.glm)
  load(glm.file)
  df.glm = cbind(sites, df.glm)
  df.sig = get_sig_sites(df.glm, comparisons = comps, FDR  , af.shifts, fdrThreshs,
                         esThreshs)
  return(df.sig)
}


df.sig.A = get.sig.sites(glm.file = 'glm.A.TP1_5.RData', rdata = 'orch2020_filtered.ACages.RData',
                         comps = c('1_5'), fdrThreshs= c(.1, .05, 0.01), esThreshs= c(0.005, 0.02, 0.02))
View(df.sig.A)
df.sig.A = read.csv("df.sig.A.1_5.csv")
df.sig.A %>% group_by(chrom) %>% count()
df.sig.A %>% group_by(sigLevel) %>% count()

2203 + 12644
14847
df.sig.B = get.sig.sites(glm.file = './glm.B.TP1_5.RData', rdata = '../RData/orch2020_filtered.BCages.RData',
                         comps = c('1_5'), fdrThreshs= c(.1, .05, 0.01), esThreshs= c(0.005, 0.02, 0.02))

df.sig.B = read.csv("df.sig.B.1_5.csv")
df.sig.B %>% group_by(chrom) %>% count()
df.sig.B %>% group_by(sigLevel) %>% count()
26056 +3247
29303  + 52023
81326
###Plotting significant sites:
df.sig = read.csv('df.sig.A.1_5.csv')

df.sig <- df.sig %>% mutate(sigLabel=as.factor(sigLevel),compLabel = as.factor(comparison))

df.sig = df.sig %>% filter(abs(afShift) > 0.005)
p.A = df.sig %>% ggplot(aes(x=pos/1000000))+ theme_minimal() + lims(y=c(0,3.5))  +
  facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
  theme(legend.position = 'none',strip.text = element_text(size=15),strip.text.y = element_text(size=12),
        axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),
        axis.title = element_text(size = 10)) +
  labs(x="Position (Mb)",y="-log10(FDR-corr. p-value)") +
  suppressWarnings(geom_point(aes(y=-log10(FDR),color=sigLabel),alpha=.8,size=1)) +
  scale_colour_manual(values = c('grey', 'orange', 'red')) +
  theme(strip.text.y = element_blank())
p.A 

#ggsave('../Figures/manhattan.A.1_5.png',p.A, height = 10, width = 10)

df.sig = read.csv('./df.sig.B.1_5.csv')
df.sig <- df.sig %>% mutate(sigLabel=as.factor(sigLevel),compLabel = as.factor(comparison))

df.sig = df.sig %>% filter(abs(afShift) > 0.005)
p.B = df.sig %>% ggplot(aes(x=pos/1000000))+ theme_minimal() + lims(y=c(0,3.5))  +
  facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
  theme(legend.position = 'none',strip.text = element_text(size=15),strip.text.y = element_text(size=12),
        axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),
        axis.title = element_text(size = 10)) +
  labs(x="Position (Mb)",y="-log10(FDR-corr. p-value)") +
  suppressWarnings(geom_point(aes(y=-log10(FDR),color=sigLabel),alpha=.8,size=1)) +
  scale_colour_manual(values = c('grey', 'orange', 'red')) +
  theme(strip.text.y = element_blank())
p.B 



#ggsave('../Figures/manhattan.B.1_5.png',p.A, height = 10, width = 10)
####ANALYSIS / PREP####
#Clusters of enriched glm signal identified in each treatment separately
clusters.a = read.csv('df.clust.Apple.1_5.csv')
clusters.b = read.csv('df.clust.Bloom.1_5.csv')
View(clusters.a)

##I took all SNPs w/ FDR < 0.2 and allele frequency change > 1% w/in each of those clusters and assesed the dynamics of the favored allele across cages in the other treatment
##I did the same with a set of matched control SNPs (matched on chromosmal arm and starting frequency)
##These files correspond to the allele frequency shifts at the target and matched control sites of each treatment,
##of clusters identified in the other treatment
##These files correspond to the allele frequency shifts at the target and matched control sites of each treatment,
##of clusters identified in the other treatment

###Use the 'phsed.shift' for assessing shifts at target sites, and matched.shift for shifts at control sites

shifts.a = read.csv('./BloomClusters.AppleShifts.csv')
shifts.b = read.csv('./AppleClusters.BloomShifts.csv')


##Get statistics of each cluster:
###Do the distributions of phased, target site shifts significantly differ from matched control shifts, for each cluster?
#Stats on Bloomington shifts in Apple clusters

data.meta = data.frame()
for(clust in (as.character(unique(shifts.b$cluster)))){
  d = data.frame()
  df.clust = shifts.b %>% filter(clust == cluster)
  cluster = clust
  median.target = median(df.clust$phased.shift)
  median.matched = median(df.clust$shift.matched)
  pval = as.numeric(t.test(df.clust$phased.shift, df.clust$shift.matched)$p.value)
  d = cbind(cluster, median.target, median.matched, pval)
  data.meta = rbind(data.meta, d)
}
data.meta = data.meta %>% mutate(pval = as.numeric(as.character(pval)))
data.meta = data.meta %>% mutate(median.target = as.numeric(as.character(median.target)))
data.meta = data.meta %>% mutate(FDR = p.adjust(pval, method = 'BH'))
view(data.meta)
write.csv(data.meta, './AppleClusters.BloomShifts.Stats.csv', row.names = FALSE)
#Stats on Apple shifts in Bloomington clusters

data.meta = data.frame()
for(clust in (as.character(unique(shifts.a$cluster)))){
  d = data.frame()
  df.clust = shifts.a %>% filter(clust == cluster)
  cluster = clust
  median.target = median(df.clust$phased.shift)
  median.matched = median(df.clust$shift.matched)
  pval = as.numeric(t.test(df.clust$phased.shift, df.clust$shift.matched)$p.value)
  d = cbind(cluster, median.target, median.matched, pval)
  data.meta = rbind(data.meta, d)
}
data.meta = data.meta %>% mutate(pval = as.numeric(as.character(pval)))
data.meta = data.meta %>% mutate(median.target = as.numeric(as.character(median.target)))
data.meta = data.meta %>% mutate(FDR = p.adjust(pval, method = 'BH'))
write.csv(data.meta, './BloomClusters.AppleShifts.Stats.csv', row.names = FALSE)

#For each treatment, what are the percentage significantly parallel, neutral, and anti-parallel cluster dynamics?
stats.a = read.csv('./BloomClusters.AppleShifts.Stats.csv')
stats.b = read.csv('./AppleClusters.BloomShifts.Stats.csv')
stats.a %>% summarise(n.parallel = sum(FDR < 0.05 & median.target > 0), n.neut = sum(FDR > 0.05 ), n.anti = sum(FDR < 0.05 & median.target < 0))
stats.b %>% summarise(n.parallel = sum(FDR < 0.05 & median.target > 0), n.neut = sum(FDR > 0.05 ), n.anti = sum(FDR < 0.05 & median.target < 0))
#what are the anti-parallel clusters for each treatment?
stats.a %>% filter(FDR < 0.05 & median.target < 0)
stats.b %>% filter(FDR < 0.05 & median.target < 0)

####HEAPMAP Plotting #### 
#heatmap of median target shift at each cluster - arranged by chromosomal arm?
stats.a$Treatment = 'LQ Shifts in HQ clusters'
stats.b$Treatment = 'HQ Shifts in LQ clusters'
View(stats.b)

####Split label column and reorder

stats.b = stats.b%>% separate(cluster, c("cluster", "chromosome", "window"))
stats.b$cluster<-gsub("^c","",as.character(stats.b$cluster))
stats.b$chromcluster <- paste(stats.b$chromosome,stats.b$cluster)


stats.a = stats.a%>% separate(cluster, c("cluster", "chromosome", "window"))
stats.a$cluster<-gsub("^c","",as.character(stats.a$cluster))
stats.a$chromcluster <- paste(stats.a$chromosome,stats.a$cluster)


stats.a.arr=stats.a%>%
  arrange(factor(chromosome, levels = c("2R","2L","3R","3L","X")))

stats.b.arr=stats.b%>%
  arrange(factor(chromosome, levels = c("2R","2L","3R","3L","X")))
range(stats.b.arr$median.target)

b.sig = subset(stats.b.arr, FDR < 0.05 & median.target < 0) ## 3R 11 & 13 apple clusters
##

a.sig = subset(stats.a.arr, FDR < 0.05 & median.target < 0) ## 3L 5 & 13 bloomington clusters

View(b.sig)
B <- ggplot(stats.a.arr, aes(x  = chromcluster, y = Treatment, fill = median.target)) +
  geom_tile() +
  scale_x_discrete(name = "Position", breaks = c("2L 1", "2R 1", "3L 1", "3R 1", "X 1"), labels = c("2L", "2R", "3L", "3R", "X")) +
  scale_fill_gradientn(name="  Allele 
  Frequency 
  Shift",
                       colours = c("blue", "white", "red"),
                       values = scales::rescale(c(-0.006046, -0.0001, 0, 0.02, 0.046)),
                       limits = c(-0.01, 0.07),
                       breaks = c(-0.01, 0.01, 0.03, 0.05)) +
  #scale_fill_gradient2(name = "Shift", high = "red", mid="white",low = "navy", limits=c(-0.006046,.046 )) +
  theme_minimal_hgrid() +
  scale_y_discrete(name = "")+
  theme(legend.text=element_text(size=8),legend.title = element_text(size=8),axis.text.x=element_text(size=8),
        axis.title = element_text(size = 8), axis.text.y = element_blank()) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 7 , direction = "vertical")) +
  #geom_text(data = subset(stats.a.arr, median.target > 0), aes(x = chromcluster, y = Treatment, label = "P"), vjust = -9)+
  geom_text(data = subset(stats.a.arr, FDR < 0.05 & median.target < 0), aes(x = chromcluster, y = Treatment, label = "*"), vjust = 0, size=8)
#geom_text(data = subset(stats.a.arr, FDR > 0.05), aes(x = chromcluster, y = Treatment, label = "N"), vjust = 0) 
#min(stats.b.arr$median.target)
B
A <- ggplot(stats.b.arr, aes(x = chromcluster, y = Treatment, fill = median.target)) +
  geom_tile() +
  scale_x_discrete(name = "Position", breaks = c("2L 1", "2R 1", "3L 1", "3R 1", "X 1"), labels = c("2L", "2R", "3L", "3R", "X")) +
  scale_fill_gradientn(name="",
                       colours = c("blue", "white", "red"),
                       values = scales::rescale(c(-0.006046, -0.0001, 0, 0.02, 0.046)),
                       limits = c(-0.01, 0.07),
                       breaks = c(-0.01, 0.01, 0.03, 0.05)) +
  #scale_fill_gradient2(name = "Shift", high = "blue", mid="white",low = "darkred", limits=c(-0.0091,.0676 )) +
  theme_minimal_hgrid() +
  scale_y_discrete(name = "")+
  theme(legend.text=element_text(size=8),axis.text.x=element_text(size=8),
        axis.title = element_text(size = 8), axis.text.y = element_blank()) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 7, direction = "vertical")) +
  #geom_text(data = subset(stats.a.arr, median.target > 0), aes(x = chromcluster, y = Treatment, label = "P"), vjust = -9)+
  geom_text(data = subset(stats.b.arr, FDR < 0.05 & median.target < 0), aes(x = chromcluster, y = Treatment, label = "*"), vjust = 0, size=8) #+
#geom_text(data = subset(stats.b.arr, FDR > 0.05), aes(x = chromcluster, y = Treatment, label = "N"), vjust = 0)
A
AB<-ggarrange(B + rremove("xlab"), A, ncol=1, heights = c(1,1.2),common.legend = TRUE, legend = c("right") )
AB

####GLM with all samples through time and treatment*time interaction####


chroms=c("2L","2R","3L","3R","X")
comparisons = c("tpt", "treat", "int")
load('orch2020_filtered.RData')


library('dplyr')
library('tidyr')
library('data.table')
library('doMC')


##Functions

fit_GLM_Treatment = function (afMatrix, rdMatrix, sampleData, model.vars = "treatment", 
                              poolCt = 100, ncores) 
{
  registerDoMC(ncores)
  df = as.data.frame(sampleData[, colnames(sampleData) %in% 
                                  model.vars])
  colnames(df) <- model.vars
  formulaString = paste0(colnames(df), collapse = " + ")
  formulaString = paste0("cts ~ ", formulaString)
  cat("Model Formula is: \n", formulaString, "\n")
  Neff = calc_Neff(rdMatrix, poolCt)
  do.call(rbind, mclapply(1:nrow(afMatrix), function(ix) {
    if (ix%%10000 == 0) {
      cat("working on site ", ix, "\n")
    }
    cts = cbind(t(round(Neff[ix, ] * afMatrix[ix, ])), t(round(Neff[ix, 
    ] * (1 - afMatrix[ix, ]))))
    df$cts = cts
    df = as.data.frame(df)
    model = glm(formulaString, family = "quasibinomial", 
                data = df)
    cp = summary(model)$coefficients[-1, c(1, 4), drop = FALSE]
    results = c(cp[, 1], cp[, 2])
    return(results)
  }, mc.cores = ncores))
}


calc_Neff = function (rd, poolCt) 
{
  ((poolCt * 2 * rd) - 1)/(poolCt * 2 + rd)
}



treats = c('A', 'B')
#load('orch2020.RData')

df.glm = sites

#load('./RData/orch2020.RData')
afmat = cbind(samps, t(afmat))
samps = samps %>% mutate(treatment = as.factor(treatment))
eec = cbind(samps, t(eec))
afmat = afmat %>% filter(treatment %in% treats)
eec = eec %>% filter(treatment %in% treats)
samps = afmat[,1:ncol(samps)]
afmat = afmat[,-c(1:ncol(samps))]
afmat = as.data.frame(t(afmat))
eec = eec[,-c(1:ncol(samps))]
eec = as.data.frame(t(eec))
res = fit_GLM_Treatment(afMatrix = afmat, rdMatrix = eec, sampleData = samps, model.vars = 'treatment',
                        poolCt = 100, ncores = 16)
res = as.data.frame(res)
names(res) = c(paste0("coef.", treats[1], "_", treats[2]), paste0("pval.", treats[1], "_", treats[2]))
df.glm = cbind(df.glm, res)



save(df.glm, file = "./05_GLM/glm.AppleVBloom.RData")
#t1->5 
load('./glm.AppleVBloom.T1_5.RData')
head(df.glm)
df.glm=df.glm %>% mutate(fdr.tpt = p.adjust(p.tpt, method = "BH"))

df.glm = df.glm %>% mutate(fdr.int = p.adjust(p.int, method = "BH"))

dim(df.glm %>% filter(fdr.tpt < 0.1))
dim(df.glm %>% filter(fdr.int < 0.1))

fdr.1_5 = as.data.frame(get_glm_FDR.V2(df.glm))
write.csv(fdr.1_5, './fdr.AppleVBloom.T1_5.csv', row.names = FALSE)
fdr.1_5 = read.csv('./fdr.AppleVBloom.T1_5.csv')


fdr.1_5 = read.csv('fdr.AppleVBloom.T1_5.csv')
load('orch2020_filtered.RData')

fdr.1_5 = cbind(sites, fdr.1_5)
head(fdr.1_5)

fdr.1_5 = fdr.1_5 %>% dplyr::select(chrom, pos, p.int) %>% rename(FDR = p.int)
fdr.1_5$comparison = 'Treatment_Interaction'
fdr.1_5 <- fdr.1_5 %>% mutate(compLabel = as.factor(comparison))
fdr.1_5 = fdr.1_5 %>% mutate(sigLabel = if_else(FDR < 0.1, 'sig', 'non-sig')) %>% 
  mutate(sigLabel=as.factor(sigLabel))
fdr.1_5 = fdr.1_5 %>% filter(FDR < 0.9)

rows(fdr.1_5)
470706
p.int = fdr.1_5 %>% ggplot(aes(x=pos/1000000))+ theme_minimal() + lims(y=c(0,2))  +
  facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
  theme(legend.position = 'none',strip.text = element_text(size=15),strip.text.y = element_text(size=12),
        axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),
        axis.title = element_text(size = 10)) +
  labs(x="Position (Mb)",y="-log10(FDR-corr. p-value)") +
  suppressWarnings(geom_point(aes(y=-log10(FDR),color=sigLabel),alpha=.8,size=4)) +
  scale_colour_manual(values = c('dark grey', 'red')) +
  theme(strip.text.y = element_blank())
p.int

chroms = c('3L','3R')

p.int.sub = fdr.1_5 %>%
  filter(chrom %in% chroms) %>%
  ggplot(aes(x=pos/1000000))+ theme_minimal() + lims(y=c(0,2))  +
  facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
  theme(legend.position = 'none',strip.text = element_text(size=15),strip.text.y = element_text(size=12),
        axis.text.y=element_text(size=12), axis.text.x=element_text(size=12),
        axis.title = element_text(size = 10)) +
  labs(x="Position (Mb)",y="-log10(FDR-corr. p-value)") +
  suppressWarnings(geom_point(aes(y=-log10(FDR),color=sigLabel),alpha=.8,size=2)) +
  scale_colour_manual(values = c('dark grey', 'red')) +
  theme(strip.text.y = element_blank())
p.int.sub

Combined2 <-ggarrange(pAB, p.A + rremove("xlab"),p.B, ncol=1, nrow = 3, heights = c(1,1,1), labels = c("","C", "D"))
Combined2

View(fdr.1_5)

SIG_fdr.1_5= fdr.1_5 %>%
  filter(sigLabel == "sig")
View(SIG_fdr.1_5)

###Get trajectories of interaction-significant alleles:
int.alleles = fdr.1_5 %>% filter(FDR < 0.1)
int.alleles = int.alleles %>% mutate(snp = paste0(chrom, pos))
##apple trjaectory generation
load('../RData/orch2020_filtered.ACages.RData', verbose = TRUE)
df = cbind(sites, afmat)
df = df %>% mutate(snp = paste0(chrom, pos))
df = df %>% filter(snp %in% int.alleles$snp) #filter frequencies to sites of interest
df = df %>% dplyr::select(-snp)

shifts.int.apple = get_af_shifts(afmat = df[,-c(1:2)], samps = samps, comparison = c('1_5'))

shifts.int.apple = shifts.int.apple %>% mutate(sign.shift = sign(dAF.1_5))

load('../RData/orch2020_filtered.ACages.RData', verbose = TRUE)
df = cbind(sites, afmat)
df = df %>% mutate(snp = paste0(chrom, pos))
df = df %>% filter(snp %in% int.alleles$snp) #filter frequencies to sites of interest
df = df %>% dplyr::select(-snp)

shifts.int.apple = get_af_shifts(afmat = df[,-c(1:2)], samps = samps, comparison = c('1_5'))

shifts.int.apple = shifts.int.apple %>% mutate(sign.shift = sign(dAF.1_5))


traj.apple = get_af_shifts(df %>% dplyr::select(-chrom, -pos), samps,cage_set = NULL, comparisons = c('1_2', '2_3','3_4','4_5'))

#Generate trajectories:
shifts.apple = matrix(nrow = nrow(df), ncol = 5)
for (row in 1:nrow(traj.apple)){
  shifts = as.vector(as.numeric(traj.apple[row, ]))
  shifts = cumsum(shifts)
  shifts = c(0, shifts)
  shifts.apple = rbind(shifts.apple, shifts)
  
}

shifts.apple = as.data.frame(na.omit(shifts.apple))
names(shifts.apple) = c('t1', 't2', 't3','t4','t5')
shifts.apple = cbind((df %>% dplyr::select(chrom, pos)), shifts.apple)
head(shifts.apple)

##alternative trjaectory generation
load('../RData/orch2020_filtered.BCages.RData', verbose = TRUE)
df = cbind(sites, afmat)
df = df %>% mutate(snp = paste0(chrom, pos))
df = df %>% filter(snp %in% int.alleles$snp) #filter frequencies to sites of interest
df = df %>% dplyr::select(-snp)

shifts.int.bloom = get_af_shifts(afmat = df[,-c(1:2)], samps = samps, comparison = c('1_5'))

shifts.int.bloom = shifts.int.bloom %>% mutate(sign.shift = sign(dAF.1_5))


traj.bloom = get_af_shifts(df %>% dplyr::select(-chrom, -pos), samps,cage_set = NULL, comparisons = c('1_2', '2_3','3_4','4_5'))

#Generate trajectories:
shifts.bloom = matrix(nrow = nrow(df), ncol = 5)
for (row in 1:nrow(traj.bloom)){
  shifts = as.vector(as.numeric(traj.bloom[row, ]))
  shifts = cumsum(shifts)
  shifts = c(0, shifts)
  shifts.bloom = rbind(shifts.bloom, shifts)
  
}

shifts.bloom = as.data.frame(na.omit(shifts.bloom))
names(shifts.bloom) = c('t1', 't2', 't3','t4','t5')
shifts.bloom = cbind((df %>% dplyr::select(chrom, pos)), shifts.bloom)
head(shifts.bloom)
shifts.apple = shifts.apple %>% gather(3:7, key = tpt, value = frequency)
shifts.bloom = shifts.bloom %>% gather(3:7, key = tpt, value = frequency)
shifts.apple$treatment = 'Apple'
shifts.bloom$treatment = 'Bloomington'
shifts.apple = shifts.apple %>% mutate(snp = paste0(chrom, pos))
shifts.bloom = shifts.bloom %>% mutate(snp = paste0(chrom, pos))
shifts.meta = rbind(shifts.apple, shifts.bloom)\
write.csv(shifts.meta, './shifts.InteractionSites.csv', row.names = FALSE)

ggplot(shifts.apple, aes(x = tpt, y = frequency, group = snp, colour = treatment)) +
  geom_line(size = 0.8, alpha = 0.8) +
  theme_bw(base_size = 35) +
  scale_colour_manual(values = c('grey')) +
  theme(legend.position = 'none')

ggplot(shifts.bloom, aes(x = tpt, y = frequency, group = snp, colour = treatment)) +
  geom_line(size = 0.8, alpha = 0.8) +
  theme_bw(base_size = 35) +
  scale_colour_manual(values = c('grey')) +
  theme(legend.position = 'none')

shifts.meta = read.csv('shifts.InteractionSites.csv')
shifts.meta = shifts.meta %>% mutate(tpt = case_when(
  tpt == 't1' ~ 1,
  tpt == 't2' ~ 2,
  tpt == 't3' ~ 3,
  tpt == 't4' ~ 4,
  tpt == 't5' ~ 5,))

shifts.meta = shifts.meta %>% mutate(treatment = case_when(
  treatment == 'Apple' ~ 'LQ',
  treatment == 'Bloomington' ~ 'HQ',))
View(shifts.meta)

p.int.line = ggplot(shifts.meta, aes(x = tpt, y = frequency, color = treatment, group = interaction(snp, treatment))) +
  geom_line(linewidth=.5) +
  labs(
    x = "Sample Time Point",
    y = "Allele Frequency", 
    color = " "
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_colour_manual(values = c('black', 'lightgray'))

p.int.line
CombinedEF <-ggarrange(p.int.sub + rremove("xlab"),p.int.line, ncol=2, nrow = 1, widths = c(1,1), labels = c("E", "F"))
CombinedEF
#plotlist = list( ,Combined2,CombinedEF)
CombinedFULL<-ggarrange(Combined2, CombinedEF, ncol=1, nrow = 2, heights = c(3,1))
CombinedFULL

#####Leave Treatment out, parallelism####
##Testing apple sites in bloom and vice versa:
df.a = read.csv('./df.sig.Matched.A.1_5.csv')
dim(df.a)
5684712
max(df.a$FDR)
0.199981805393494
#apple sites in bloom
df.a = read.csv('./df.sig.Matched.A.1_5.csv')
df.a = df.a %>% mutate(sign.shift = sign(afShift))
df.a = df.a %>% filter(sigLevel > 1)
target.a.sites = df.a %>% dplyr::select(chrom, pos)
matched.a.sites = df.a %>% dplyr::select(chrom, pos.matched) %>% rename(pos = pos.matched)
cages = c(1:6)
for (c in cages){
  load(paste0('../RData/RData_ByCage/orch2020.B',c, '.RData'))
  shifts = get_af_shifts(afmat = afmat.lo, samps = samps.lo, cage_set = NULL, comparison = c('1_5'))
  shifts = cbind(sites, shifts)
  shifts.target = left_join(target.a.sites, shifts)
  shifts.matched = left_join(matched.a.sites, shifts)
  shifts.target = shifts.target %>% mutate(dAF.1_5.target = dAF.1_5 * df.a$sign.shift )
  shifts.matched = shifts.matched %>% rename(dAF.1_5.matched = dAF.1_5) %>% rename(pos.matched = pos)
  shifts = cbind(shifts.target %>% dplyr::select(chrom, pos,dAF.1_5.target ), 
                 shifts.matched %>% dplyr::select(pos.matched, dAF.1_5.matched))
  shifts$cage = c
  shifts$id_treatment = 'Apple'
  shifts$test_treatment = 'Bloom'
  write.csv(shifts, paste0('./Apple.Sites.BloomCage',c,'Shifts.csv'), row.names = FALSE)
}

#bloom sites in apple
df.b = read.csv('./df.sig.Matched.B.1_5.csv')
df.b = df.b %>% mutate(sign.shift = sign(afShift))
df.b = df.b %>% filter(sigLevel > 1)
target.b.sites = df.b %>% dplyr::select(chrom, pos)
matched.b.sites = df.b %>% dplyr::select(chrom, pos.matched) %>% rename(pos = pos.matched)
cages = c(1:6)
for (c in cages){
  load(paste0('../RData/RData_ByCage/orch2020.A',c, '.RData'))
  shifts = get_af_shifts(afmat = afmat.lo, samps = samps.lo, cage_set = NULL, comparison = c('1_5'))
  shifts = cbind(sites, shifts)
  shifts.target = left_join(target.b.sites, shifts)
  shifts.matched = left_join(matched.b.sites, shifts)
  shifts.target = shifts.target %>% mutate(dAF.1_5.target = dAF.1_5 * df.b$sign.shift )
  shifts.matched = shifts.matched %>% rename(dAF.1_5.matched = dAF.1_5) %>% rename(pos.matched = pos)
  shifts = cbind(shifts.target %>% dplyr::select(chrom, pos,dAF.1_5.target ), 
                 shifts.matched %>% dplyr::select(pos.matched, dAF.1_5.matched))
  shifts$cage = c
  shifts$id_treatment = 'Bloom'
  shifts$test_treatment = 'Apple'
  write.csv(shifts, paste0('./Bloom.Sites.AppleCage',c,'Shifts.csv'), row.names = FALSE)
}
files = list.files(pattern = 'Shifts')
head(read.csv('Bloom.Sites.AppleCage5Shifts.csv'))
head(read.csv('Bloom.Sites.AppleCage4Shifts.csv'))

#genome wide stats:
data.stats.genome = data.frame()
files = list.files(pattern = 'Shifts')
for (f in files){
  d = read.csv(f)
  id.treat = strsplit(f, '[.]')[[1]][1]
  test.treat = strsplit(strsplit(f, '[.]')[[1]][3], 'Cage' )[[1]][1]
  test.cage = as.character(d %>% distinct(cage))
  chrom = 'Genome'
  median.target = median(d$dAF.1_5.target)
  median.matched = median(d$dAF.1_5.matched)
  pval = t.test(d$dAF.1_5.target, d$dAF.1_5.matched)$p.val
  stats = as.data.frame(cbind(id.treat, test.treat, test.cage,chrom, median.target, median.matched, pval))
  data.stats.genome = rbind(data.stats.genome, stats)
}
#chromsome-level stats:
data.stats.chr = data.frame()
files = list.files(pattern = 'Shifts')
for (f in files){
  d = read.csv(f)
  id.treat = strsplit(f, '[.]')[[1]][1]
  test.treat = strsplit(strsplit(f, '[.]')[[1]][3], 'Cage' )[[1]][1]
  test.cage = as.character(d %>% distinct(cage))
  stats.cage = data.frame()
  for (chr in as.character(unique(d$chrom))){    
    d.c = d %>% filter(chrom == chr)
    c = chr
    median.target = median(d.c$dAF.1_5.target)
    median.matched = median(d.c$dAF.1_5.matched)
    pval = t.test(d.c$dAF.1_5.target, d.c$dAF.1_5.matched)$p.val
    stats = as.data.frame(cbind(id.treat, test.treat, test.cage, c, median.target, median.matched, pval))
    stats.cage = rbind(stats.cage, stats)    
  }
  data.stats.chr = rbind(data.stats.chr, stats.cage)
}
data.stats.chr = data.stats.chr %>% rename(chrom = c)
stats.meta = rbind(data.stats.genome, data.stats.chr)
stats.meta = unique(stats.meta)
stats.meta = stats.meta %>% mutate(median.target = as.numeric(as.character(median.target)),
                                   median.matched = as.numeric(as.character(median.matched)),
                                   pval = as.numeric(as.character(pval)))

stats.meta = stats.meta %>% rowwise() %>% mutate(FDR = p.adjust(pval))
stats.meta = stats.meta %>% mutate(sig = if_else(FDR < 0.05 & abs(median.target) > abs(median.matched), 'sig', 'non-sig'))
stats.meta = stats.meta %>% rowwise() %>%
  mutate(dynamic = case_when(sig == 'sig' & abs(median.target > 0) ~ 'parallel',
                             sig == 'sig' & abs(median.target < 0) ~ 'fluctuating',
                             sig == 'non-sig' & abs(median.target > 0) ~ 'non-sig',))
write.csv(stats.meta, './GLM.TreatmentSiteTest.Stats.csv', row.names = FALSE)

stats.meta = read.csv('./GLM.TreatmentSiteTest.Stats.csv')
stats.meta = stats.meta %>% gather(median.target, median.matched, key = type, value = median)
head(stats.meta)

stats.meta = stats.meta %>% rowwise() %>%
  mutate(SNPs = strsplit(as.character(type), '[.]')[[1]][2]) %>%
  mutate(Color = case_when(SNPs == 'matched' ~ 'black',
                           SNPs == 'target' & dynamic == 'parallel' ~ 'red',
                           SNPs == 'target'& dynamic == 'fluctuating' ~ 'blue',
                           SNPs == 'target'& dynamic == 'non-sig' ~ 'dark grey'))
stats.meta$chrom = factor(stats.meta$chrom, levels = c('Genome', '2L', '2R', '3L', '3R', 'X'))
stats.meta = stats.meta %>% mutate(SNPs = if_else(SNPs == 'matched', 'Matched Control', 'Target'))
p = ggplot(stats.meta %>% filter(id.treat == 'Apple'), aes(x = SNPs, y = median, shape = SNPs, colour = Color))+
  geom_jitter(width = 0.2, size = 2) +
  scale_colour_identity() + 
  facet_wrap(~chrom, ncol=6) +
  scale_shape_manual(values = c(4, 19)) +
  ylim(c(-0.02, 0.07)) +
  scale_x_discrete(limits = c('Target', 'Matched Control')) +
  theme_bw(base_size = 10) +
  xlab('')+ 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) 
p
ggsave('../Figures/AppleSites.BloomCageTest.pdf', p, height = 10, width = 15)


q = ggplot(stats.meta %>% filter(id.treat == 'Bloom'), aes(x = SNPs, y = median, shape = SNPs, colour = Color))+
  geom_jitter(width = 0.2, size = 2) +
  scale_colour_identity() + 
  facet_wrap(~chrom, ncol=6) +
  ylim(c(-0.02, 0.07)) +
  scale_shape_manual(values = c(4, 19)) +
  scale_x_discrete(limits = c('Target', 'Matched Control')) +
  theme_bw(base_size = 10) +
  xlab('')+ 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
q
ggsave('../Figures/BloomSites.AppleCageTest.pdf', p, height = 10, width = 15)

pq <-ggarrange(p,q, ncol=1, labels = c("E", "F"), common.legend = TRUE)
pq

Combined2

p.int.line

CombinedGH <-ggarrange(p.int.sub + rremove("xlab"),p.int.line, ncol=2, nrow = 1, widths = c(1,1), labels = c("G", "H"))
CombinedGH


CombinedFULL<-ggarrange(Combined2, pq, CombinedEF, ncol=1, nrow = 3, heights = c(3,2,1))
CombinedFULL




