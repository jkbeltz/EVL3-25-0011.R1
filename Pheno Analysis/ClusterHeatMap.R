setwd("/Users/jackbeltz/Documents/PENN/Dissertation/CH 2:3(Orchard 20)/AB/RAW/raw geno")
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

shifts.a = read.csv('BloomClusters.AppleShifts.csv')
shifts.b = read.csv('AppleClusters.BloomShifts.csv')


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
stats.a = read.csv('BloomClusters.AppleShifts.Stats.csv')
stats.b = read.csv('AppleClusters.BloomShifts.Stats.csv')
stats.a %>% summarise(n.parallel = sum(FDR < 0.05 & median.target > 0), n.neut = sum(FDR > 0.05 ), n.anti = sum(FDR < 0.05 & median.target < 0))
stats.b %>% summarise(n.parallel = sum(FDR < 0.05 & median.target > 0), n.neut = sum(FDR > 0.05 ), n.anti = sum(FDR < 0.05 & median.target < 0))

#what are the anti-parallel clusters for each treatment?
stats.a %>% filter(FDR < 0.05 & median.target < 0)
stats.b %>% filter(FDR < 0.05 & median.target < 0)

####HEAPMAP Plotting #### 
#heatmap of median target shift at each cluster - arranged by chromosomal arm?
stats.a$Treatment = 'LQ Shifts in HQ clusters'
stats.b$Treatment = 'HQ Shifts in LQ clusters'
#View(stats.b)

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
View(a.sig)
B <-ggplot(stats.a.arr, aes(x  = chromcluster, y = Treatment, fill = median.target)) +
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
AB<-ggarrange(B + rremove("xlab"), A, ncol=1, heights = c(1,1.2), labels = c("A","B"),common.legend = TRUE, legend = c("right") )
AB


####PCAS ####
pca_data=read.csv("PCAResults_AppleVBloom.csv")
pca_data_noout=read.csv("PCAResults_AppleVBloom.T5B4_Rem.csv")
head(pca_data)

p1 = ggplot(pca_data, aes(x = PC1, y = PC2, color = Timepoint, shape = treatment)) +
  geom_point(size = 3) +
  scale_colour_gradient2(low = "red", mid= "lightblue",high = "blue", midpoint = 3 , limits = c(1, 5), ) +
  scale_shape_manual(name="",breaks =c("A", "B"), values=c(17,19),labels=c("Low Quality","High Quality"))+
  theme_minimal(base_size=8)


p2 = ggplot(pca_data, aes(x = PC1, y = PC3, color = treatment, shape = treatment )) +
  geom_point(size = 3) +
  scale_shape_manual(name="",breaks =c("A", "B"), values=c(17,19),labels=c("Low Quality","High Quality"))+
  #scale_colour_gradient2(low = "red", mid= "white",high = "blue", midpoint = 3 , limits = c(1, 5), ) +
  scale_color_manual(name="",breaks =c("A", "B"),values = c("grey","black"), labels=c("Low Quality","High Quality"))+
  theme_minimal(base_size=8)


p<-ggarrange(p1, p2, ncol=2, widths = c(1,1), labels = c("A","B"), common.legend = TRUE, legend=c('right') )
p

q1= ggplot(pca_data_noout, aes(x = PC1, y = PC2, color = Timepoint, shape = treatment)) +
  geom_point(size = 3) +
  scale_colour_gradient2(low = "red", mid= "lightblue",high = "blue", midpoint = 3 , limits = c(1, 5), ) +
  scale_shape_manual(breaks =c("A", "B"), values=c(15,19),labels=c("Low","High"))+
  theme_minimal()

q2 = ggplot(pca_data_noout, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 3) +
  #scale_colour_gradient2(low = "red", mid= "white",high = "blue", midpoint = 3 , limits = c(1, 5), ) +
  scale_color_manual(breaks =c("A", "B"),values = c("red","black"), labels=c("Low","High"))+
  theme_minimal()

q<-ggarrange(q1, q2, ncol=2, widths = c(1,1), labels = c("C","D"), common.legend = FALSE, legend=c('bottom') )
q
####HEATMAP####

#Apple cages t1 -> 5 Manhattan Plot
df.sigA = read.csv('./df.sig.A.1_5.csv')
df.sigA <- df.sigA %>% mutate(sigLabel=as.factor(sigLevel),compLabel = as.factor(comparison))

df.sigA = df.sigA %>% filter(abs(afShift) > 0.005)
p.A = df.sigA %>% ggplot(aes(x=pos/1000000))+ theme_minimal() + lims(y=c(0,4))  +
  facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
  theme(legend.position = 'none',strip.text = element_text(size=8),strip.text.y = element_text(size=8),
        axis.text.y=element_text(size=8), axis.text.x=element_text(size=8),
        axis.title = element_text(size = 8)) +
  labs(x="Position (Mb)",y="-log10(FDR-corr. p-value)") +
  suppressWarnings(geom_point(aes(y=-log10(FDR),color=sigLabel),alpha=.8,size=.6)) +
  scale_colour_manual(values = c('grey', 'orange', 'red')) +
  theme(strip.text.y = element_blank())
p.A 
ggsave('../Figures/manhattan.A.1_5.pdf',p.A, height = 10, width = 10)


#Bloomington cages t1 -> 5 Manhattan Plot
df.sigB = read.csv('./df.sig.B.1_5.csv')
df.sigB <- df.sigB %>% mutate(sigLabel=as.factor(sigLevel),compLabel = as.factor(comparison))

df.sigB = df.sigB %>% filter(abs(afShift) > 0.005)
p.B = df.sigB %>% ggplot(aes(x=pos/1000000))+ theme_minimal() + lims(y=c(0,4))  +
  facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
  theme(legend.position = 'none',strip.text = element_text(size=8),strip.text.y = element_text(size=8),
        axis.text.y=element_text(size=8), axis.text.x=element_text(size=8),
        axis.title = element_text(size = 8)) +
  labs(x="Position (Mb)",y="-log10(FDR-corr. p-value)") +
  suppressWarnings(geom_point(aes(y=-log10(FDR),color=sigLabel),alpha=.8,size=.6)) +
  scale_colour_manual(values = c('grey', 'orange', 'red')) +
  theme(strip.text.y = element_blank())
p.B 
ggsave('../Figures/manhattan.B.1_5.pdf',p.B, height = 10, width = 10)

####combining plots####

Combined <-ggarrange(p,p.B + rremove("xlab"),p.A,AB, ncol=1, nrow = 4, heights = c(.8,1,1,1), labels = c("", "C", "D",""))
Combined

Combined_no_outlier<-ggarrange(q,AB, ncol=1, nrow = 2, heights = c(1, .8), labels = c("","C"))
Combined_no_outlier

