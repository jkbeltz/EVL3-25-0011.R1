#### orchard 2020 16s analysis #########################
## started 9/14/21 by jk beltz
## post qiime2 analysis ########
install.packages("ggpubr")
library("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

BiocManager::install("dada2", version = "3.17")

install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.17")

#!/usr/bin/env Rscript
library("tidyverse")
library(qiime2R)
library(tidyverse)
library(Biostrings)
library(dplyr)
install.packages("dada2")
library(dada2)
library(qiime2R)
library(ggplot2)
library(ggpubr)
library(plyr)
library(cowplot)
library(nlme)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
library(ggbreak)
library(patchwork)
library(ggpubr)

### online tutorial https://www.google.com/search?q=qiimer+vs.+qiime2r&oq=qiimer+vs.+qiime2r&aqs=chrome..69i57j0i22i30.11627j0j7&sourceid=chrome&ie=UTF-8
##Reading QZA files#######
##SET WD to HOME/DOCUMENTS/16S_QIIME2

###removed wolbachia in qiime, see "useful Qiime Commands"###
##only removed 1 sample from analysis and it was a control!

observed_features_nw=read_qza("raw/O20_Qiime_Metrics/observed_features_vector.qza")
observed_features_nw<-observed_features_nw$data %>% rownames_to_column("SampleID")
write.csv(observed_features_nw, "raw/O20_Qiime_Metrics/observed_features_nw.csv")

##reload all data with wolbachia samples removed

evenvess_nw=scan("raw/O20_Qiime_Metrics/evenness_vector.qza") ##error 
shannon_nw=read_qza("raw/O20_Qiime_Metrics/shannon_vector.qza")
shannon_nw<-shannon_nw$data %>% rownames_to_column("SampleID")
#View(shannon_nw)
observed_features_nw=read_qza("raw/O20_Qiime_Metrics/observed_features_vector.qza")
observed_features_nw<-observed_features_nw$data %>% rownames_to_column("SampleID")
faith_pd_nw=read_qza("raw/O20_Qiime_Metrics/faith_pd_vector.qza")
faith_pd_nw<-faith_pd_nw$data %>% rownames_to_column("SampleID")
rarefied_nw=read_qza("raw/O20_Qiime_Metrics/rarefied_table.qza")

## cant rename column b/c sample names are variables 
uu_distance_matrix_nw=read_qza("raw/O20_Qiime_Metrics/unweighted_unifrac_distance_matrix.qza")
uu_distance_matrix_nw=uu_distance_matrix_nw$data
bray_distance_matrix_nw=read_qza("raw/O20_Qiime_Metrics/bray_curtis_distance_matrix.qza")
bray_distance_matrix_nw=bray_distance_matrix_nw$data
bray_curtis_pcoa_nw=read_qza("raw/O20_Qiime_Metrics/bray_curtis_pcoa_results.qza")
bray_curtis_pcoa_nw=bray_curtis_pcoa_nw$data
jaccard_matrix_nw=read_qza("raw/O20_Qiime_Metrics/jaccard_distance_matrix.qza")
jaccard_matrix_nw=jaccard_matrix_nw$data
jaccard_pcoa_nw=read_qza("raw/O20_Qiime_Metrics/jaccard_pcoa_results.qza")
jaccard_pcoa_nw=jaccard_pcoa_nw$data
uu_pcoa_nw=read_qza("raw/O20_Qiime_Metrics/unweighted_unifrac_pcoa_results.qza")
uu_pcoa_nw=uu_pcoa_nw$data
wu_distance_matrix_nw=read_qza("raw/O20_Qiime_Metrics/weighted_unifrac_distance_matrix.qza")
wu_distance_matrix_nw=wu_distance_matrix_nw$data
wu_pcoa_nw=read_qza("raw/O20_Qiime_Metrics/weighted_unifrac_pcoa_results.qza")
#View(wu_pcoa_nw)

####reread otutable 

SVs_nw=read_qza("raw/O20_Qiime_Metrics/table-no-wolbachia-exact.qza")
OTUTABLE_nw=SVs_nw$data
typeof(OTUTABLE_nw)
View(taxon)

#remerge old tacxon table with new otu values
otu.taxon_nw<- merge(taxon, OTUTABLE_nw, by ='row.names', all=TRUE) ## successful merge, names double checked

#View(otu.taxon_nw)
######metadata work#############
metadata=read.csv(file = "raw/O20_Qiime_Metrics/orchard2016smetadata.csv", header = TRUE)
colnames(metadata)[1]= "SampleID"
write.csv(metadata, "raw/O20_Qiime_Metrics/orchard2016smetadata.csv")
View(metadata)
metadata[,20]=tolower(metadata[,20])
metadata[,22]=tolower(metadata[,22])

### add alpha diveristy metrics to meta###################
a=metadata %>% 
  merge(observed_features_nw, all.y=TRUE) 
b=a %>% 
  merge(faith_pd_nw, all.y=TRUE) 
allmeta=b %>% 
  merge(shannon_nw, all.y=TRUE) 
View (allmeta)
alphadata = allmeta
alphadata <- allmeta[ -c(2:7, 9:13, 23:36, 39:49) ]
alphadata[101,5]="Axenic_LB"
View(alphadata)

alphadatano5<-alphadata[(alphadata$time_point == "T2") | (alphadata$time_point == "T4"),]

###subsetting####
APBloomalpha<-alphadatano5[(alphadatano5$cage_treatment == "Apple") | (alphadatano5$cage_treatment == "Bloom"),]
APATLBalpha<-alphadatano5[(alphadatano5$cage_treatment == "Apple") | (alphadatano5$cage_treatment == "Apple_AT") | (alphadatano5$cage_treatment == "Apple_LB"),]

#cagefood
APBloomfoodalpha=APBloomalpha[(APBloomalpha$sample_type_detail== "cage food"),]
APATLBfoodalpha=APATLBalpha[(APATLBalpha$sample_type_detail== "cage food"),]

#wild flies 
APBloomwildalpha=APBloomalpha[(APBloomalpha$sample_type_detail== "flies"),]
APATLBwildalpha=APATLBalpha[(APATLBalpha$sample_type_detail== "flies"),]

#F1s
APBloomF1salpha=APBloomalpha[(APBloomalpha$sample_type_detail== "dry stored f1s"),]
APATLBF1salpha=APATLBalpha[(APATLBalpha$sample_type_detail== "dry stored f1s"),]


alphacageenv= subset(alphadata, alphadata$SampleType=="Environmental control")
alphaallflies=subset(alphadata, alphadata$SampleType=="Whole fly")
alphawildflies=subset(alphadata, alphadata$sample_type_detail== "flies" )
#alphawildF1s=subset(alphadata, alphadata$sample_type_detail= c("flies", "dry stored f1s"))
alphawildfood<-alphadata[(alphadata$sample_type_detail== "flies") | (alphadata$sample_type_detail == "cage food"),]
alphawildF1s <-alphadata[(alphadata$sample_type_detail == "flies") | (alphadata$sample_type_detail == "dry stored f1s"),]
alphawildfoodf1s <-alphadata[(alphadata$sample_type_detail == "flies") | (alphadata$sample_type_detail =="dry stored f1s") | (alphadata$sample_type_detail== "cage food"),]
alphaf2s=subset(alphadata, alphadata$sample_type_detail == "dry stored f2s")

alphaallflies$cage_treatment= as.factor(alphaallflies$cage_treatment)
alphaallflies$time_point=as.factor(alphaallflies$time_point)
alphaallflies$sample_type_detail=as.factor(alphaallflies$sample_type_detail)
###line plotting#####

alphawildflies$cage_treatment= as.factor(alphawildflies$cage_treatment)
alphawildflies$time_point=as.factor(alphawildflies$time_point)
alphawildflies$sample_type_detail=as.factor(alphawildflies$sample_type_detail)
alphawildflies$study_group=as.factor(alphawildflies$study_group)

#View(APATLBF1salpha)

#### best plot, just change dataframe #################
View(APBloomF1salpha)

OF_F1_AB=APBloomF1salpha%>%
  filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=observed_features, group=interaction(as.factor(cage_treatment)), color=as.factor(cage_treatment))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=6, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,25)+
  xlab("Time Point") +
  ylab("OTU Richness")+
  ggtitle("F1 Common Garden Flies")+
  #scale_x_discrete(labels= Xlabs)+
  scale_colour_manual(name = "Cage Treatment", 
    breaks = c("Apple","Bloom"),
    labels = c("LQ","HQ"),
    values = c("grey","black"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
  axis.title=element_text(size=14))+
  theme(legend.position="none")
OF_F1_AB

OF_CF_AB=APBloomfoodalpha%>%
  filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=observed_features, group=interaction(as.factor(cage_treatment)), color=as.factor(cage_treatment))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=6, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,25)+
  xlab("Time Point") +
  ylab("OTU Richness")+
  ggtitle("      Cage Food Supply")+
  #scale_x_discrete(labels= Xlabs)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(legend.position="none")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
OF_CF_AB

OF_W_AB=APBloomwildalpha%>%
  #filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=observed_features, group=interaction(as.factor(cage_treatment)), color=as.factor(cage_treatment))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=6, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,25)+
  xlab("Time Point") +
  ylab("OTU Richness")+
  ggtitle("       Wild Collected Flies")+
  scale_x_discrete(labels= Xlabs)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(legend.position="none")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

SH_F1_AB=APBloomF1salpha%>%
  #filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=shannon_entropy, group=interaction(as.factor(cage_treatment)), color=as.factor(cage_treatment))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=6, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,5)+
  xlab("Time Point") +
  ylab("Shannon Diversity")+
  scale_x_discrete(labels= Xlabs)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(legend.position="none") 

SH_CF_AB=APBloomfoodalpha%>%
  filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=shannon_entropy, group=interaction(as.factor(cage_treatment)), color=as.factor(cage_treatment))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=6, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,5)+
  xlab("Time Point") +
  ylab("Shannon Diversity")+
  scale_x_discrete(labels= Xlabs)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(legend.position="none") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

SH_W_AB=APBloomwildalpha%>%
  #filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=shannon_entropy, group=interaction(as.factor(cage_treatment)), color=as.factor(cage_treatment))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=6, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,5)+
  xlab("Time Point") +
  ylab("Shannon Diversity")+
  scale_x_discrete(labels= Xlabs)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(legend.position="none")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

FP_F1_AB=APBloomF1salpha%>%
  #filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=faith_pd, group=interaction(as.factor(cage_treatment)), color=as.factor(cage_treatment))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=6, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,6)+
  xlab("Time Point") +
  ylab("Phylogenetic Distance")+
  scale_x_discrete(labels= Xlabs)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(legend.position="none")

FP_CF_AB=APBloomfoodalpha%>%
  #filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=faith_pd, group=interaction(as.factor(cage_treatment)), color=as.factor(cage_treatment))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=6, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,6)+
  xlab("Time Point") +
  ylab("Phylogenetic Distance")+
  scale_x_discrete(labels= Xlabs)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(legend.position="none")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

FP_W_AB=APBloomwildalpha%>%
  #filter (!sample_type_detail == 'cage env') %>%
  #filter (!sample_type_detail == 'wild flies') %>%
  ggplot(aes(x=as.factor(time_point), y=faith_pd, group=interaction(as.factor(cage_treatment)), color=as.factor(cage_treatment))) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se, size=6, shape=19) +
  geom_jitter(shape=19, width=0.15, height=0)+
  ylim(0,6)+
  xlab("Time Point") +
  ylab("Phylogenetic Distance")+
  scale_x_discrete(labels= Xlabs)+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(legend.position="none")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

SUPPbac.div <- ggarrange(OF_F1_AB + rremove("xlab"), OF_CF_AB + rremove("xlab") , OF_W_AB + rremove("xlab") ,
                         SH_F1_AB + rremove("xlab"), SH_CF_AB + rremove("xlab") , SH_W_AB + rremove("xlab") ,
                         FP_F1_AB , FP_CF_AB , FP_W_AB ,
                  labels = c("A", " B", " C", "D", " E", " F", "G", " H", " I"),
                  ncol = 3, nrow = 3, common.legend = FALSE, font.label=list(color="black",size=16, face = "bold"))
SUPPbac.div


#########################PCAs#############################################################################
metadata<-read_q2metadata("metric/sample-metadata.tsv")
uwunifrac<-read_qza("metrics/unweighted_unifrac_pcoa_results.qza")
shannon<-read_qza("metrics/shannon_vector.qza")$data %>% rownames_to_column("SampleID") 

uualphameta=uu_pcoa_nw$Vectors %>%
  select(SampleID, PC1, PC2, PC3, PC4) %>%
  left_join(alphadata)

wualphameta=wu_pcoa_nw$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3, PC4) %>%
  left_join(alphadata)

View(uualphameta)
View(alphadata)
##subsets #########


###VIZ PCAS

#### THIS PCA AB ####
uualphawildcagef1s=subset(uualphameta, uualphameta$sample_type_detail %in% c("flies","cage food","dry stored f1s"))
uualphawildcagef1sAB=subset(uualphawildcagef1s, uualphawildcagef1s$cage_treatment %in% c("Apple","Bloom"))

wualphawildcagef1s=subset(wualphameta, wualphameta$sample_type_detail %in% c("flies","cage food","dry stored f1s"))
wualphawildcagef1sAB=subset(wualphawildcagef1s, wualphawildcagef1s$cage_treatment %in% c("Apple","Bloom"))


WU_AB=wualphawildcagef1sAB%>%
  filter (!cage_treatment == 'Apple_AT') %>% 
  filter (!cage_treatment == 'Apple_LB') %>% 
  filter (!cage_treatment == 'base Bloom') %>% 
  filter (!sample_type_detail == 'cage env') %>% 
  filter (!time_point == 'T5') %>% 
  filter (!time_point == 'T0') %>% 
  mutate(sample_type_detail=recode(sample_type_detail,'cage food'= "Cage Food Supply",'dry stored f1s'="F1 Common Garden Flies ", 'flies'= "Wild Collected Flies"))

WU_AB_plot=WU_AB%>%
  ggplot(aes(x=PC1, y=PC2, color= cage_treatment, shape=time_point, group=cage_treatment)) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + 
  stat_ellipse(type='t',size =1)+#alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_wrap(. ~sample_type_detail)+
  scale_shape_manual(values=c(19,8),name="Timepoint", labels = c("TP2","TP4")) + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes 8,0,15,5,1,16
  scale_size_continuous(name="Shannon Diversity") +
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_bw()  


WU_AB_plot


UU_AB=uualphawildcagef1sAB%>%
  filter (!cage_treatment == 'Apple_AT') %>% 
  filter (!cage_treatment == 'Apple_LB') %>% 
  filter (!cage_treatment == 'base Bloom') %>% 
  filter (!sample_type_detail == 'cage env') %>% 
  filter (!time_point == 'T5') %>% 
  filter (!time_point == 'T0') %>% 
  mutate(sample_type_detail=recode(sample_type_detail,'cage food'= "Cage Food Supply",'dry stored f1s'="F1 Common Garden Flies ", 'flies'= "Wild Collected Flies"))

UU_AB_plot=UU_AB%>%
  ggplot(aes(x=PC1, y=PC2, color= cage_treatment, shape=time_point, group=cage_treatment)) +     #, size=shannon_entropy
  geom_point(alpha=1, size=4) + 
  stat_ellipse(type='t',size =1)+#alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_wrap(. ~sample_type_detail)+
  scale_shape_manual(values=c(19,8),name="Timepoint", labels = c("TP2","TP4")) + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes 8,0,15,5,1,16
  scale_size_continuous(name="Shannon Diversity") +
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("Apple","Bloom"),
                      labels = c("LQ","HQ"),
                      values = c("grey","black"))+
  theme_bw()  


UU_AB_plot


ggarrange(UU_AB_plot, WU_AB_plot, ncol = 1, labels =c('A','B'), common.legend = TRUE, legend="bottom")

#########################BACTERIAL ABUDNACE ON FOOD##############
##  EXPORT!!!!
alphafood24=subset(metadata, metadata$SampleType== "Environmental control" )

a=as.data.frame(colSums(OTUTABLE_nw[,-1]))
colnames(a)[1] = "readcounts"
library(tibble)
readcounts <- tibble::rownames_to_column(a, "SampleID")

alphadata=alphadata %>% 
  merge(readcounts, all.y=TRUE) 

#alphadata = metadata
alphadata <- transform(alphadata, readcountsadjusted = alphadata$readcounts / alphadata$final_library_conc_ng_ul)
alphadata = subset(alphadata, !(study_day %in% c("None", NA)))

alphadatanoout <- subset(alphadata, alphadata$readcountsadjusted < 2000)

alphadatanoout %>%
  as.data.frame() %>%
  filter(`sample_type_detail`== "cage food")  %>%
  ggplot() +
  #geom_boxplot(aes(x=time_point, y=observed_features, fill=cage_treatment)) +
  geom_point(aes(x=time_point , y=observed_features, group=study_group, fill=cage_treatment))

### THIS PLOT 
View(alphadatanoout)
alphadatanoout1=alphadatanoout %>%
  as.data.frame() %>%
  filter (!sample_type_detail == 'dry stored f2s')%>%
  filter (!cage_treatment == 'Apple_LB') %>% 
  filter (!cage_treatment == 'Apple_AT') %>% 
  filter (!cage_treatment == 'base Bloom') %>% 
  filter (!sample_type_detail == 'cage env') %>% 
  filter (!time_point == 'T5') %>%   
  filter (!time_point == 'T0') %>%   
  mutate(sample_type_detail=recode(sample_type_detail,'cage food'= "Cage Food Supply",'dry stored f1s'="F1 Common Garden Flies ", 'flies'= "Wild Collected Flies"))%>% 
  mutate(time_point = recode(time_point,'T2'= "TP2",'T4'="TP4 "))

alphadatanoout_plot=alphadatanoout1 %>%
  ggplot(aes(x=time_point, y=readcountsadjusted, fill=cage_treatment, group=interaction(time_point, cage_treatment))) +
  geom_boxplot()  + 
  ylim(0,1000) +
  facet_wrap(. ~sample_type_detail)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  ylab("Total Bacterial Abundance")+
  xlab(" ")+
  scale_fill_manual(name = "Cage Treatment", 
                   breaks = c("Apple","Bloom"),
                   labels = c("Low Quality","High Quality"),
                   values = c("grey", "black"))+
  theme_bw()  

alphadatanoout_plot

######ALTERNATIVE APPROACH top 15 read abundance######
library(tidyverse)
library(dplyr)
library(qiime2R)

metadata<-alphadata
alphadataSVs<-read_qza("table-no-wolbachia-exact.qza")$data
taxonomy<-read_qza("correcttaxonomy.qza")$data
library(tibble)
taxonomy <- tibble::rownames_to_column(taxonomy, "Feature.ID")

SVs<-apply(SVs, 2, function(x) x/sum(x)*100) #convert to percent
View(SVs)

SVsToPlot<-  
  data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(15, MeanAbundance) %>%
  pull(Feature.ID) #extract only the names from the table

BIGSVs=SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  #summarize(Abundance=sum(Abundance)) %>%
  left_join(alphadata, by= "SampleID") %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  left_join(taxonomy, by= "Feature.ID") %>%
  mutate(Feature=paste(Feature.ID, Taxon)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature))  # trim out leading text from taxonomy string
View(BIGSVs)

SmallSVs=BIGSVs[!(BIGSVs$Feature.ID == "Remainder"),]
test=SmallSVs[sample(nrow(SmallSVs), 10), ]
test$Feature=gsub(".*Bacteria;","",as.character(test$Feature)) ####fixes trimming issue of feature variable

SmallSVs=BIGSVs[!(BIGSVs$Feature.ID == "Remainder"),]##remove remainder
SmallSVs$Feature=gsub(".*Bacteria;","",as.character(SmallSVs$Feature))

BIGSVs$Feature=gsub(".*Bacteria;","",as.character(BIGSVs$Feature))

##wild flies treatmentxtimepoint## 
wildSVs=subset(BIGSVs, BIGSVs$sample_type_detail == "flies", na.rm=TRUE)

wildSVs=wildSVs%>%
  mutate(cage_treatment=recode(cage_treatment,'Apple'= "LQ",'Bloom'="HQ"))%>%
  mutate(time_point=recode(time_point,'T2'= "TP2",'T4'="TP4"))

a=ggplot(subset(wildSVs,cage_treatment %in% c("LQ" , "HQ")), aes(x=study_group, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`cage_treatment` + `time_point`, scales="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
  theme(legend.position = "bottom") +
  ggtitle("Wild Collected Flies") +
  xlab("Cage Treatment x Season")+ 
  ylab(" ")+
  scale_fill_viridis_c(option = "plasma" , name="Abundance",  )+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

a 
#ggsave("wildflyheatmap.pdf  ", height=4, width=10, device="pdf") # save a PDF ? inches by ?  inches



##cage food treatmentxtimepoint## 
foodSVs=subset(BIGSVs, BIGSVs$sample_type_detail == "cage food")
foodSVs=foodSVs%>%
  mutate(cage_treatment=recode(cage_treatment,'Apple'= "LQ",'Bloom'="HQ"))%>%
  mutate(time_point=recode(time_point,'T2'= "TP2",'T4'="TP4"))


b= ggplot(subset(foodSVs,cage_treatment %in% c("LQ" , "HQ")), aes(x=study_group, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`cage_treatment` + `time_point`, scales="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
  theme(legend.position = "bottom") +
  ggtitle("Cage Food Supply") +
  xlab("Cage Treatment x Time Point")+ 
  ylab(" ")+
  scale_fill_viridis_c(option = "plasma" , name="Abundance",  )+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


b

##f1s treatment###
f1sSVs=subset(BIGSVs, BIGSVs$sample_type_detail == "dry stored f1s")
f1sSVs=subset(f1sSVs, f1sSVs$time_point %in% c("T2","T4"),)

f1sSVs=f1sSVs%>%
  mutate(cage_treatment=recode(cage_treatment,'Apple'= "LQ",'Bloom'="HQ"))%>%
  mutate(time_point=recode(time_point,'T2'= "TP2",'T4'="TP4"))


c= ggplot(subset(f1sSVs,cage_treatment %in% c("LQ" , "HQ")), aes(x=study_group, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`cage_treatment` + `time_point`, scales="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1, face="bold")) +
  theme(legend.position = "bottom") +
  ggtitle("F1 Common Garden Flies") +
  ylab(" ")+
  xlab("Cage Treatment x Season") + 
  scale_fill_viridis_c(option = "plasma" , name="Abundance",  )


c

ggarrange(a + rremove("xlab"), b + rremove("xlab"),c + rremove("xlab"), ncol=1, common.legend = TRUE, legend = "right")

