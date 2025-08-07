####Figure 2####
install.packages("patchwork")

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

####Import####
LD = read.csv("raw/O20AB_LD_raw.csv")
StarvM = read.csv("raw/O20AB_SRm_raw.csv")
StarvF = read.csv("raw/O20AB_SRf_raw.csv")
Fecund = read.csv("raw/O20AB_F_raw.csv")
Viab = read.csv("raw/O20AB_V_raw.csv")
Thorax = read.csv("raw/O20AB_TL_raw.csv")

####Make Summary dfs####

######LD####
##by treatment / timepoint##
#View(LD)
LD.sum=LD %>%
  group_by(Group) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   LD=mean(Mean.LD,na.rm = TRUE),
                   LDsd=sd(Mean.LD,na.rm = TRUE),
                   LDse=LDsd/sqrt(N),
                   LDm=mean(Mean.M.LD,na.rm = TRUE),
                   LDmsd=sd(Mean.M.LD,na.rm = TRUE),
                   LDmse=LDmsd/sqrt(N),
                   LDf=mean(Mean.F.LD,na.rm = TRUE),
                   LDfsd=sd(Mean.F.LD,na.rm = TRUE),
                   LDfse=LDfsd/sqrt(N),
  )
LD.sum=unique(LD.sum)

##by cage / timepoint / treatment##
LD.cage.sum=LD %>%
  group_by(Group2) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage=as.factor(Cage),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   LD=mean(Mean.LD,na.rm = TRUE),
                   LDsd=sd(Mean.LD,na.rm = TRUE),
                   LDse=LDsd/sqrt(N),
                   LDm=mean(Mean.M.LD,na.rm = TRUE),
                   LDmsd=sd(Mean.M.LD,na.rm = TRUE),
                   LDmse=LDmsd/sqrt(N),
                   LDf=mean(Mean.F.LD,na.rm = TRUE),
                   LDfsd=sd(Mean.F.LD,na.rm = TRUE),
                   LDfse=LDfsd/sqrt(N),
  )
LD.cage.sum=unique(LD.cage.sum)
#View(LD.cage.sum)

######SRM####
##by treatment / time point##
#View(StarvM)
SRM.sum=StarvM %>%
  group_by(Group) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   SRM=mean(Mean.SRm,na.rm = TRUE),
                   SRMsd=sd(Mean.SRm,na.rm = TRUE),
                   SRMse=SRMsd/sqrt(N),
  )
SRM.sum=unique(SRM.sum)
#View(SRM.sum)

##by cage / timepoint / treatment##
SRM.cage.sum=StarvM %>%
  group_by(Group2) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage=as.factor(Cage),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   SRM=mean(Mean.SRm,na.rm = TRUE),
                   SRMsd=sd(Mean.SRm,na.rm = TRUE),
                   SRMse=SRMsd/sqrt(N),
  )
SRM.cage.sum=unique(SRM.cage.sum)
#View(SRM.cage.sum)

######SRF####
##by treatment / time point##
#View(StarvF)
SRF.sum=StarvF %>%
  group_by(Group) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   SRF=mean(Mean.SRf,na.rm = TRUE),
                   SRFsd=sd(Mean.SRf,na.rm = TRUE),
                   SRFse=SRFsd/sqrt(N),
                   SRFCV=( SRFsd/SRF)*100,
  )
SRF.sum=unique(SRF.sum)
View(SRF.sum)

##by cage / timepoint / treatment##
SRF.cage.sum=StarvF %>%
  group_by(Group2) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage=as.factor(Cage),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   SRF=mean(Mean.SRf,na.rm = TRUE),
                   SRFsd=sd(Mean.SRf,na.rm = TRUE),
                   SRFse=SRFsd/sqrt(N),
  )
SRF.cage.sum=unique(SRF.cage.sum)

######Fecundity####
#View(Fecund)

Fec.sum=Fecund %>%
  group_by(Group) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   Fec=mean(Mean.eggperday,na.rm = TRUE),
                   Fecsd=sd(Mean.eggperday,na.rm = TRUE),
                   Fecse= Fecsd/sqrt(N),
                   FECCV=(Fecsd/Fec)*100,
  )
Fec.sum=unique(Fec.sum)
View(Fec.sum)

##by cage/ timepoint/ treatment##

Fec.cage.sum=Fecund %>%
  group_by(Group2) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage=as.factor(Cage),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   Fec=mean(Mean.eggperday,na.rm = TRUE),
                   Fecsd=sd(Mean.eggperday,na.rm = TRUE),
                   Fecse= Fecsd/sqrt(N),
                   TOTF=mean(totalfecund,na.rm = TRUE),
                   TOTFsd=sd(totalfecund,na.rm = TRUE),
                   TOTFse=TOTFsd/sqrt(N),
  )
Fec.cage.sum=unique(Fec.cage.sum)
#View(Fec.cage.sum)

######Viability####
#View(Viab)
Via.sum=Viab %>%
  group_by(Group) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   Viab=mean(Viability,na.rm = TRUE),
                   Viabsd=sd(Viability,na.rm = TRUE),
                   Viabse= Viabsd/sqrt(N),
                   ViabCV=( Viabsd/Viab)*100,
  )
Via.sum=unique(Via.sum)
View(Via.sum)

##by cage/ timepoint/ treatment##

Via.cage.sum=Viab %>%
  group_by(Group2) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   Cage=as.factor(Cage),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   Viab=mean(Viability,na.rm = TRUE),
                   Viabsd=sd(Viability,na.rm = TRUE),
                   Viabse= Viabsd/sqrt(N),
             
  )
Via.cage.sum=unique(Via.cage.sum)
#View(Via.cage.sum)

######Body Size######

#View(Thorax)
BS.sum=Thorax %>%
  group_by(Group) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   BS=mean(Mean.TL,na.rm = TRUE),
                   BSsd=sd(Mean.TL,na.rm = TRUE),
                   BSse= BSsd/sqrt(N),
                   BSCV=(BSsd/BS)*100,
  )
BS.sum=unique(BS.sum)
View(BS.sum)

##by cage/ timepoint/ treatment##

BS.cage.sum=Thorax %>%
  group_by(Group2) %>%
  dplyr::summarise(Group=as.factor(Group),
                   N=n(),
                   TP=as.factor(Timepoint),
                   Cage.Treat=as.factor(Cage.Treat),
                   Pheno.Treat=as.factor(Pheno.Treat),
                   Cage=as.factor(Cage),
                   Cage.Pheno=as.factor(Cage.Pheno),
                   BS=mean(Mean.TL,na.rm = TRUE),
                   BSse=mean(SE.TL,na.rm = TRUE),
                   
  )
BS.cage.sum=unique(BS.cage.sum)
#View(BS.cage.sum)






####Merging Dfs####
##by tp/treatment##

Phenos.list <- list(LD.sum, SRM.sum, SRF.sum, Fec.sum, Via.sum, BS.sum)      
#View(Phenos.list)
#merge all data frames together
Phenos=Phenos.list %>% reduce(full_join, by=c('Group','TP','Cage.Treat','Pheno.Treat','Cage.Pheno'))
View(Phenos)
##Remove N columns
Phenos <- Phenos[, -c(2,16,20,25,30,35)]
Phenos <- Phenos[-17,]
View(Phenos)

write.csv (Phenos, "Phenos.csv", row.names=FALSE)
Phenos=read.csv("Phenos.csv")
#View(Phenos) #### All phenotypic data in one table by treatment and time point

##by tp/cage/treatment##

Phenos.cage.list <- list(LD.cage.sum, SRM.cage.sum, SRF.cage.sum, Fec.cage.sum, Via.cage.sum, BS.cage.sum)      

#merge all data frames together
Phenos.cage=Phenos.cage.list %>% reduce(full_join, by=c('Group2','Group','TP','Pheno.Treat','Cage','Cage.Pheno'))


##Remove N columns
Phenos.cage <- Phenos.cage[, -c(3,18,19,23,24,28,29,36,37,41,42)]
Phenos.cage <- Phenos.cage[-c(87:88),]
#View(Phenos.cage) #### All phenotypic data in one table by cage, treatment and time point
write.csv (Phenos.cage, "Phenos.cage.csv", row.names=FALSE)
Phenos.cage=read.csv("Phenos.cage.csv")
##Subsetting##
founder= subset(Phenos, Phenos$TP == 0)
founder.cage=subset(Phenos.cage, Phenos.cage$TP == 0)
#View(founder)
Pheno=subset(Phenos, Phenos$TP %in% c(1,2,3,4,5))
Pheno.cage=subset(Phenos.cage, Phenos.cage$TP %in% c(1,2,3,4,5))
#View(Pheno.cage)

Phenos24=subset(Phenos, Phenos$TP %in% c(2,4))
Phenos.cage24=subset(Phenos.cage, Phenos.cage$TP%in% c(2,4))

pheno = subset(Pheno, Pheno$Cage.Pheno %in% c("aa","bb"))
pheno.cage = subset(Phenos.cage, Phenos.cage$Cage.Pheno %in% c("aa","bb"))
#View(pheno.cage)



pheno.cage.pca<-prcomp(pheno.cage[c(11,14,17,20,26,29,32)],scale=TRUE)
pheno.cage.pca.df<-data.frame(Cage=pheno.cage$Cage,Time=pheno.cage$TP,Treat=pheno.cage$Cage.Treat.x,pheno.cage.pca$x)
#View(pheno.cage.pca.df)

####PCA PHENOTYPE####

##no shapes by time
Pheno.pca=ggplot(pheno.cage.pca.df, aes(x=PC1, y=PC2, color = Treat, shape = Treat, group=Treat)) +
  geom_point(size=5)+
  stat_ellipse()+
  scale_colour_manual(name = "Nutritional Environement", 
                      breaks = c("a", "b"),
                      labels = c("Low Quality","High Quality"),
                      values = c("grey","black"))+
  scale_shape_manual(name = "Nutritional Environement", 
                      breaks = c("a", "b"),
                      labels = c("Low Quality","High Quality"),
                      values = c(19,17))+
  theme_bw()+
  theme(
                       legend.position = c(.99, .99),
                       legend.justification = c("right", "top"),
                       legend.box.just = "right",
                       legend.margin = margin(.5, .5, .5,.5), 
                       legend.text=element_text(size=8),
                       legend.title = element_blank()
                     )
pheno.cage.pca.df.manov<-manova(cbind(PC1,PC2)~Treat,pheno.cage.pca.df)
summary(pheno.cage.pca.df.manov)
Pheno.pca

Pheno.pca_time=ggplot(pheno.cage.pca.df, aes(x=PC1, y=PC2, shape=as.factor(Time), color = Treat, group=Treat)) +
  geom_point(size=5)+
  stat_ellipse()+
  scale_colour_manual(name = "Nutritional Environement", 
                      breaks = c("a", "b"),
                      labels = c("Low Quality","High Quality"),
                      values = c("grey","black"))+
  theme_bw()+
  scale_shape_manual(breaks=c("0","1","2","3","4","5"), 
                     label=c("Founder","Timepoint 1","Timepoint 2","Timepoint 3","Timepoint 4","Timepoint 5"), 
                     values =c(8,15,16,17,6,7))+ 
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6), 
    legend.text=element_text(size=8)
  )

Pheno.pca_time
#View(pheno)
aa.pheno.cage = subset(Phenos.cage, Phenos.cage$Cage.Pheno %in% c("aa"))
bb.pheno.cage = subset(Phenos.cage, Phenos.cage$Cage.Pheno %in% c("bb"))

aa.pheno24.cage = subset(Phenos.cage24, Phenos.cage24$Pheno.Treat %in% c("a"))
bb.pheno24.cage = subset(Phenos.cage24, Phenos.cage24$Pheno.Treat %in% c("b"))

#View(aa.pheno24.cage)
####GRAPHING Fig 2 ####

##Development Time Male + Female Fig##

LD1 <- ggplot(pheno, aes(x=as.factor(TP), y=LD, shape=Pheno.Treat, group=Pheno.Treat))+
  geom_point(data = founder, aes(x=as.factor(TP), y=LD,color = Pheno.Treat), size=5)+
  geom_linerange(data=founder, aes(ymin=LD-LDse,ymax=LD+LDse, color=Pheno.Treat), size=1)+
  geom_point(aes(color=Pheno.Treat),size=5)+
  geom_line(aes(color=Pheno.Treat))+
  geom_linerange(aes(ymin=LD-LDse,ymax=LD+LDse, color=Pheno.Treat), size=1)+
  ylab("Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("F","t1", "t2", "t3", "t4", "t5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("a", "b"),
                      labels = c("Low","High"),
                      values = c("grey","black"))+
  geom_line(data=aa.pheno.cage, aes(x=as.factor(TP), y=LD, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4) + 
  geom_line(data=bb.pheno.cage, aes(x=as.factor(TP), y=LD, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
LD1

## Male Development Time Fig##

LDm1 <- ggplot(pheno, aes(x=as.factor(TP), y=LDm,shape=Pheno.Treat, group=Pheno.Treat))+
  geom_point(data = founder, aes(x=as.factor(TP), y=LDm,color = Pheno.Treat), size=5)+
  geom_linerange(data=founder, aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=Pheno.Treat), size=1)+
  geom_point(aes(color=Pheno.Treat),size=5)+
  geom_line(aes(color=Pheno.Treat))+
  geom_linerange(aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=Pheno.Treat), size=1)+
  ylab("Male Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("F","t1", "t2", "t3", "t4", "t5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("a", "b"),
                      labels = c("Low","High"),
                      values = c("grey","black"))+
  geom_line(data=aa.pheno.cage, aes(x=as.factor(TP), y=LDm, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4) + 
  geom_line(data=bb.pheno.cage, aes(x=as.factor(TP), y=LDm, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
LDm1

## Female Development Time Fig##

LDf1 <- ggplot(pheno, aes(x=as.factor(TP), y=LDf, shape=Pheno.Treat, group=Pheno.Treat))+
  geom_point(data = founder, aes(x=as.factor(TP), y=LDf,color = Pheno.Treat), size=5)+
  geom_linerange(data=founder, aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=Pheno.Treat), size=1)+
  geom_point(aes(color=Pheno.Treat),size=5)+
  geom_line(aes(color=Pheno.Treat))+
  geom_linerange(aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=Pheno.Treat), size=1)+
  ylab("Female Development Time (hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("F","t1", "t2", "t3", "t4", "t5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("a", "b"),
                      labels = c("Low","High"),
                      values = c("grey","black"))+
  geom_line(data=aa.pheno.cage, aes(x=as.factor(TP), y=LDf, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4) + 
  geom_line(data=bb.pheno.cage, aes(x=as.factor(TP), y=LDf, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
LDf1

##Starvation Resistance Female##

SRf1 <- ggplot(pheno, aes(x=as.factor(TP), y=SRF, shape=Pheno.Treat, group=Pheno.Treat))+
  geom_point(data = founder, aes(x=as.factor(TP), y=SRF,color = Pheno.Treat), size=5)+
  geom_linerange(data=founder, aes(ymin=SRF-SRFse,ymax=SRF+SRFse, color=Pheno.Treat), size=1)+
  geom_point(aes(color=Pheno.Treat),size=5)+
  geom_line(aes(color=Pheno.Treat))+
  geom_linerange(aes(ymin=SRF-SRFse,ymax=SRF+SRFse, color=Pheno.Treat), size=1)+
  ylab("Female Starvation Restistance (Hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("F","t1", "t2", "t3", "t4", "t5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("a", "b"),
                      labels = c("Low","High"),
                      values = c("grey","black"))+
  geom_line(data=aa.pheno.cage, aes(x=as.factor(TP), y=SRF, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4) + 
  geom_line(data=bb.pheno.cage, aes(x=as.factor(TP), y=SRF, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4)+
 theme_cowplot(10)+ 
  theme(legend.position="none")
SRf1

##Starvation Resistance Male##

SRm1 <- ggplot(pheno, aes(x=as.factor(TP), y=SRM, shape=Pheno.Treat, group=Pheno.Treat))+
  geom_point(data = founder, aes(x=as.factor(TP), y=SRM,color = Pheno.Treat), size=5)+
  geom_linerange(data=founder, aes(ymin=SRM-SRMse,ymax=SRM+SRMse, color=Pheno.Treat), size=1)+
  geom_point(aes(color=Pheno.Treat),size=5)+
  geom_line(aes(color=Pheno.Treat))+
  geom_linerange(aes(ymin=SRM-SRMse,ymax=SRM+SRMse, color=Pheno.Treat), size=1)+
  ylab("Male Starvation Restistance (Hrs)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("F","t1", "t2", "t3", "t4", "t5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("a", "b"),
                      labels = c("Low","High"),
                      values = c("grey","black"))+
  geom_line(data=aa.pheno.cage, aes(x=as.factor(TP), y=SRM, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4) + 
  geom_line(data=bb.pheno.cage, aes(x=as.factor(TP), y=SRM, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+
  theme(legend.position="none")
SRm1

##Fecundity per day##

Fec1 <- ggplot(pheno, aes(x=as.factor(TP), y=Fec, shape=Pheno.Treat, group=Pheno.Treat))+
  geom_point(data = founder, aes(x=as.factor(TP), y=Fec,color = Pheno.Treat), size=5)+
  geom_linerange(data=founder, aes(ymin=Fec-Fecse,ymax=Fec+Fecse, color=Pheno.Treat), size=1)+
  geom_point(aes(color=Pheno.Treat),size=5)+
  geom_line(aes(color=Pheno.Treat))+
  geom_linerange(aes(ymin=Fec-Fecse,ymax=Fec+Fecse, color=Pheno.Treat), size=1)+
  ylab("Fecundity (Eggs / Day)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("F","t1", "t2", "t3", "t4", "t5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("a", "b"),
                      labels = c("Low","High"),
                      values = c("grey","black"))+
  geom_line(data=aa.pheno.cage, aes(x=as.factor(TP), y=Fec, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4) + 
  geom_line(data=bb.pheno.cage, aes(x=as.factor(TP), y=Fec, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
Fec1

##Fecundity Total ##

TOTF1 <- ggplot(pheno, aes(x=as.factor(TP), y=TOTF, shape=Pheno.Treat, group=Pheno.Treat))+
  geom_point(data = founder, aes(x=as.factor(TP), y=TOTF,color = Pheno.Treat), size=5)+
  geom_linerange(data=founder, aes(ymin=TOTF-TOTFse,ymax=TOTF+TOTFse, color=Pheno.Treat), size=1)+
  geom_point(aes(color=Pheno.Treat),size=5)+
  geom_line(aes(color=Pheno.Treat))+
  geom_linerange(aes(ymin=TOTF-TOTFse,ymax=TOTF+TOTFse, color=Pheno.Treat), size=1)+
  ylab("Fecundity (Total)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("F","t1", "t2", "t3", "t4", "t5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("a", "b"),
                      labels = c("Low Nutririon","High"),
                      values = c("grey","black"))+
  geom_line(data=aa.pheno.cage, aes(x=as.factor(TP), y=TOTF, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4) + 
  geom_line(data=bb.pheno.cage, aes(x=as.factor(TP), y=TOTF, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
TOTF1

##Viability ##

Viab1 <- ggplot(pheno, aes(x=as.factor(TP), y=Viab, shape=Pheno.Treat, group=Pheno.Treat))+
  geom_point(data = founder, aes(x=as.factor(TP), y=Viab,color = Pheno.Treat), size=5)+
  geom_linerange(data=founder, aes(ymin=Viab-Viabse,ymax=Viab+Viabse, color=Pheno.Treat), size=1)+
  geom_point(aes(color=Pheno.Treat),size=5)+
  geom_line(aes(color=Pheno.Treat))+
  geom_linerange(aes(ymin=Viab-Viabse,ymax=Viab+Viabse, color=Pheno.Treat), size=1)+
  ylab("Egg to Adult Viability")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("F","t1", "t2", "t3", "t4", "t5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("a", "b"),
                      labels = c("Low Nutririon","High"),
                      values = c("grey","black"))+
  ylim(0.5,1.1)+
  geom_line(data=aa.pheno.cage, aes(x=as.factor(TP), y=Viab, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4) + 
  geom_line(data=bb.pheno.cage, aes(x=as.factor(TP), y=Viab, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+ 
  theme(legend.position="none")
Viab1

##Body Size ##

BS1 <- ggplot(pheno, aes(x=as.factor(TP), y=BS, shape=Pheno.Treat, group=Pheno.Treat))+
  geom_point(data = founder, aes(x=as.factor(TP), y=BS,color = Pheno.Treat), size=5)+
  geom_linerange(data=founder, aes(ymin=BS-BSse,ymax=BS+BSse, color=Pheno.Treat), size=1)+
  geom_point(aes(color=Pheno.Treat),size=5)+
  geom_line(aes(color=Pheno.Treat))+
  geom_linerange(aes(ymin=BS-BSse,ymax=BS+BSse, color=Pheno.Treat), size=1)+
  ylab("Female Thorax Length (uM)")+
  xlab("Sample Timepoints")+
  scale_x_discrete(labels=c("F","t1", "t2", "t3", "t4", "t5"))+
  scale_colour_manual(name = "Cage Treatment", 
                      breaks = c("a", "b"),
                      labels = c("Low Nutririon","High"),
                      values = c("grey","black"))+
  geom_line(data=aa.pheno.cage, aes(x=as.factor(TP), y=BS, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4) + 
  geom_line(data=bb.pheno.cage, aes(x=as.factor(TP), y=BS, group=as.factor(Cage), color=Pheno.Treat), 
            size=.5, alpha= .4)+
  theme_cowplot(10)+
  theme(legend.position="none")
BS1

######Combine####

library(ggpubr)



Fig2 <- ggarrange(LDm1 + rremove("xlab"), LDf1 + rremove("xlab"), SRm1+ rremove("xlab"),SRf1+ rremove("xlab"),Viab1,BS1,TOTF1,Pheno.pca,
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                  ncol = 4, nrow = 2, common.legend = FALSE, font.label=list(color="black",size=16, face = "bold"))
Fig2
save(Fig2, file = "fig2.rdata")


######Linear Mixed Model Analysis F####
###Found in Table 2###


ABLmmdata=subset(Pheno.cage, Pheno.cage$Cage.Pheno %in% c("aa","bb"))
ABLmmdata$TP=as.numeric(ABLmmdata$TP, order=TRUE)
ABLmmdata$Pheno.Treat=as.factor(ABLmmdata$Pheno.Treat)
View(ABLmmdata)

library(tidyverse)
library(ggpubr)
library(rstatix)

#NORMALITY TEST 

ABLmmdata2=ABLmmdata %>% mutate(arcViab=asin(Viab))
shapiro_test(ABLmmdata2$logBS)

ABLmmdata %>% hist(BS)
ggqqplot(ABLmmdata$LDf)

shapiro_test(ABLmmdata$BS)

LDm1 <- lmer(LDm ~ TP*Pheno.Treat + (1|Cage), data = ABLmmdata)  # the syntax stays the same, but now the nesting is taken into account
anova(LDm1)
summary(LDm1)


pairs=lsmeans(LDm1, pairwise ~ TP:Pheno.Treat)
summ=summary(pairs$contrasts)

View(summ)
LDf1 <- lmer(LDf ~ TP*Pheno.Treat + (1|Cage) , data = ABLmmdata)  # the syntax stays the same, but now the nesting is taken into account
anova(LDf1)
summary(LDf1)

V1 <- lmer(arcViab ~ TP*Pheno.Treat + (1|Cage) , data = ABLmmdata2)  # the syntax stays the same, but now the nesting is taken into account
anova(V1)


SRM1 <- lmer(SRM ~ TP*Pheno.Treat + (1|Cage) , data = ABLmmdata)  # the syntax stays the same, but now the nesting is taken into account
anova(SRM1)
summary(SRM1)

SRF1 <- lmer(SRF ~ TP*Pheno.Treat + (1|Cage) , data = ABLmmdata)  # the syntax stays the same, but now the nesting is taken into account
anova(SRF1)
summary(SRF1)

TOTF1 <- lmer(TOTF ~ TP*Pheno.Treat + (1|Cage),  data = ABLmmdata)  # the syntax stays the same, but now the nesting is taken into account
anova(TOTF1)
1summary(TOTF1)

BS1 <- lmer(BS~ TP*Pheno.Treat + (1|Cage) + (1|TP), data = ABLmmdata)  # the syntax stays the same, but now the nesting is taken into account
anova(BS1)
summary(BS1)



















#### Reciprocal Transfer SUPP Figure####

BSCV=ggplot(Phenos24, aes(x=TP, y=BSCV, group= Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  ylab("Body Size Coefficient of Variation (%)")+
  xlab("Sample Timepoints") +
  #scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  #scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  theme_cowplot(font_size = 12) 
BSCV

VIABCV=ggplot(Phenos24, aes(x=TP, y=ViabCV, group= Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  ylab("Viability Coefficient of Variation (%)")+
  xlab("Sample Timepoints") +
  #scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  #scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  theme_cowplot(font_size = 12) 
VIABCV

FECCV=ggplot(Phenos24, aes(x=TP, y=FECCV, group= Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  ylab("Fecundity Coefficient of Variation (%)")+
  xlab("Sample Timepoints") +
  #scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  #scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  theme_cowplot(font_size = 12) 
FECCV


SRFCV=ggplot(Phenos24, aes(x=TP, y=SRFCV, group= Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  ylab("Starvtation Resistance (Female) Coefficient of Variation (%)")+
  xlab("Sample Timepoints") +
  #scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  #scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  theme_cowplot(font_size = 12) 
SRFCV

test = ggarrange(BSCV, VIABCV, SRFCV, FECCV, ncol=2, nrow=2, common.legend = TRUE)
test

view(Phenos)
#view(founder)
#View(Phenos.cage24)


Phenos24.a=subset(Phenos24,Phenos24$Cage.Treat=="a")
Phenos24.b=subset(Phenos24,Phenos24$Cage.Treat=="b")
#View(Phenos24.b)
founder.a=subset(founder,founder$Cage.Treat=="a")
founder.b=subset(founder,founder$Cage.Treat=="b")
aa.pheno24.cage.a=subset(aa.pheno24.cage, aa.pheno24.cage$Cage.Treat.x == "a")
aa.pheno24.cage.b=subset(aa.pheno24.cage, aa.pheno24.cage$Cage.Treat.x == "b")
bb.pheno24.cage.a=subset(bb.pheno24.cage, bb.pheno24.cage$Cage.Treat.x =="a")
bb.pheno24.cage.b=subset(bb.pheno24.cage, bb.pheno24.cage$Cage.Treat.x =="b")



##male development time####

LDm1a <- ggplot(Phenos24.a, aes(x=as.numeric(TP), y=LDm, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=Cage.Treat), size=1)+
  ylab("Male Development Time (hrs)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.a, aes(x=as.numeric(TP), y=LDm,color = Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.a, aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.a, aes(x=as.numeric(TP), y=LDm, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.a, aes(x=as.numeric(TP), y=LDm, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x =element_blank())

LDm1a

LDm1b <- ggplot(Phenos24.b, aes(x=as.numeric(TP), y=LDm, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=Cage.Treat), size=1)+
  ylab("Male Development Time (hrs)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,-1),limits = c(-1, 4.5), breaks = c(0,2,4),labels=c("Founder","t2","t4"))+
  scale_x_break(c(.5, 1.9), scales = 5)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.b, aes(x=as.numeric(TP), y=LDm,color=Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.b, aes(ymin=LDm-LDmse,ymax=LDm+LDmse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.b, aes(x=as.numeric(TP), y=LDm, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.b, aes(x=as.numeric(TP), y=LDm, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(legend.position="none")+
  theme(axis.title.x =element_blank())

LDm1b
LDm1= aplot::plot_list(LDm1b, LDm1a, nrow = 1, ncol = 2, heights = c(5,5,5,3))
LDm1

##female developement time####
LDf1a <- ggplot(Phenos24.a, aes(x=as.numeric(TP), y=LDf, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=Cage.Treat), size=1)+
  ylab("Female Development Time (hrs)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.a, aes(x=as.numeric(TP), y=LDf,color = Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.a, aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.a, aes(x=as.numeric(TP), y=LDf, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.a, aes(x=as.numeric(TP), y=LDf, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x =element_blank())
 

LDf1a

LDf1b <- ggplot(Phenos24.b, aes(x=as.numeric(TP), y=LDf, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=Cage.Treat), size=1)+
  ylab("Female Development Time (hrs)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,-1),limits = c(-1, 4.5), breaks = c(0,2,4),labels=c("Founder","t2","t4"))+
  scale_x_break(c(.5, 1.9), scales = 5)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.b, aes(x=as.numeric(TP), y=LDf,color=Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.b, aes(ymin=LDf-LDfse,ymax=LDf+LDfse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.b, aes(x=as.numeric(TP), y=LDf, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.b, aes(x=as.numeric(TP), y=LDf, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(legend.position="none")+
  theme(axis.title.x =element_blank())
  

LDf1b
LDf1= aplot::plot_list(LDf1b, LDf1a, nrow = 1, ncol = 2, heights = c(5,5,5,3))
LDf1

##Male Starvation Resistance####
SRM1a <- ggplot(Phenos24.a, aes(x=as.numeric(TP), y=SRM, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=SRM-SRMse,ymax=SRM+SRMse, color=Cage.Treat), size=1)+
  ylab("Male Starvation Resistance (hrs)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.a, aes(x=as.numeric(TP), y=SRM,color = Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.a, aes(ymin=SRM-SRMse,ymax=SRM+SRMse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.a, aes(x=as.numeric(TP), y=SRM, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.a, aes(x=as.numeric(TP), y=SRM, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x =element_blank())

SRM1a

SRM1b <- ggplot(Phenos24.b, aes(x=as.numeric(TP), y=SRM, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=SRM-SRMse,ymax=SRM+SRMse, color=Cage.Treat), size=1)+
  ylab("Male Starvation Resistance (hrs)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,-1),limits = c(-1, 4.5), breaks = c(0,2,4),labels=c("Founder","t2","t4"))+
  scale_x_break(c(.5, 1.9), scales = 5)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.b, aes(x=as.numeric(TP), y=SRM,color=Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.b, aes(ymin=SRM-SRMse,ymax=SRM+SRMse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.b, aes(x=as.numeric(TP), y=SRM, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.b, aes(x=as.numeric(TP), y=SRM, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(legend.position="none")+
  theme(axis.title.x =element_blank())

SRM1b
SRM1= aplot::plot_list(SRM1b, SRM1a, nrow = 1, ncol = 2, heights = c(5,5,5,3))
SRM1

##Female Starvation Resistance####

SRF1a <- ggplot(Phenos24.a, aes(x=as.numeric(TP), y=SRF, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=SRF-SRFse,ymax=SRF+SRFse, color=Cage.Treat), size=1)+
  ylab("Female Starvation Resistance (hrs)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.a, aes(x=as.numeric(TP), y=SRF,color = Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.a, aes(ymin=SRF-SRFse,ymax=SRF+SRFse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.a, aes(x=as.numeric(TP), y=SRF, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.a, aes(x=as.numeric(TP), y=SRF, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")+
  theme(axis.title.x =element_blank())

SRF1a

SRF1b <- ggplot(Phenos24.b, aes(x=as.numeric(TP), y=SRF, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=SRF-SRFse,ymax=SRF+SRFse, color=Cage.Treat), size=1)+
  ylab("Female Starvation Resistance (hrs)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,-1),limits = c(-1, 4.5), breaks = c(0,2,4),labels=c("Founder","t2","t4"))+
  scale_x_break(c(.5, 1.9), scales = 5)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.b, aes(x=as.numeric(TP), y=SRF,color=Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.b, aes(ymin=SRF-SRFse,ymax=SRF+SRFse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.b, aes(x=as.numeric(TP), y=SRF, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.b, aes(x=as.numeric(TP), y=SRF, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(legend.position="none")+
  theme(axis.title.x =element_blank())

SRF1b
SRF1= aplot::plot_list(SRF1b, SRF1a, nrow = 1, ncol = 2, heights = c(5,5,5,3))
SRF1

##Fecundity total ####
TOTF1a <- ggplot(Phenos24.a, aes(x=as.numeric(TP), y=TOTF, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=TOTF-TOTFse,ymax=TOTF+TOTFse, color=Cage.Treat), size=1)+
  ylab("Total Fecundity")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.a, aes(x=as.numeric(TP), y=TOTF,color = Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.a, aes(ymin=TOTF-TOTFse,ymax=TOTF+TOTFse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.a, aes(x=as.numeric(TP), y=TOTF, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.a, aes(x=as.numeric(TP), y=TOTF, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")

TOTF1a

TOTF1b <- ggplot(Phenos24.b, aes(x=as.numeric(TP), y=TOTF, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=TOTF-TOTFse,ymax=TOTF+TOTFse, color=Cage.Treat), size=1)+
  ylab("Total Fecundity")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,-1),limits = c(-1, 4.5), breaks = c(0,2,4),labels=c("Founder","t2","t4"))+
  scale_x_break(c(.5, 1.9), scales = 5)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.b, aes(x=as.numeric(TP), y=TOTF,color=Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.b, aes(ymin=TOTF-TOTFse,ymax=TOTF+TOTFse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.b, aes(x=as.numeric(TP), y=TOTF, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.b, aes(x=as.numeric(TP), y=TOTF, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(legend.position="none")

TOTF1b
TOTF1= aplot::plot_list(TOTF1b, TOTF1a, nrow = 1, ncol = 2, heights = c(5,5,5,3))
TOTF1

##Viability#### 
Viab1a <- ggplot(Phenos24.a, aes(x=as.numeric(TP), y=Viab, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=Viab-Viabse,ymax=Viab+Viabse, color=Cage.Treat), size=1)+
  ylab("Viability")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.a, aes(x=as.numeric(TP), y=Viab,color = Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.a, aes(ymin=Viab-Viabse,ymax=Viab+Viabse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.a, aes(x=as.numeric(TP), y=Viab, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.a, aes(x=as.numeric(TP), y=Viab, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")

Viab1a

Viab1b <- ggplot(Phenos24.b, aes(x=as.numeric(TP), y=Viab, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=Viab-Viabse,ymax=Viab+Viabse, color=Cage.Treat), size=1)+
  ylab("Viability")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,-1),limits = c(-1, 4.5), breaks = c(0,2,4),labels=c("Founder","t2","t4"))+
  scale_x_break(c(.5, 1.9), scales = 5)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.b, aes(x=as.numeric(TP), y=Viab,color=Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.b, aes(ymin=Viab-Viabse,ymax=Viab+Viabse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.b, aes(x=as.numeric(TP), y=Viab, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.b, aes(x=as.numeric(TP), y=Viab, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(legend.position="none")

Viab1b
Viab1= aplot::plot_list(Viab1b, Viab1a, nrow = 1, ncol = 2, heights = c(5,5,5,3))
Viab1


##Body Size####

BS1a <- ggplot(Phenos24.a, aes(x=as.numeric(TP), y=BS, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=BS-BSse,ymax=BS+BSse, color=Cage.Treat), size=1)+
  ylab("Female Thorax Length (um)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,0),limits = c(-.5, 4.5), breaks = c(2,4),labels=c("t2","t4"))+
  scale_x_break(c(1.0, 1.9), scales = 10)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.a, aes(x=as.numeric(TP), y=BS,color = Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.a, aes(ymin=BS-BSse,ymax=BS+BSse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.a, aes(x=as.numeric(TP), y=BS, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.a, aes(x=as.numeric(TP), y=BS, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")

BS1a

BS1b <- ggplot(Phenos24.b, aes(x=as.numeric(TP), y=BS, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=BS-BSse,ymax=BS+BSse, color=Cage.Treat), size=1)+
  ylab("Female Thorax Length (um)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,-1),limits = c(-1, 4.5), breaks = c(0,2,4),labels=c("Founder","t2","t4"))+
  scale_x_break(c(.5, 1.9), scales = 5)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.b, aes(x=as.numeric(TP), y=BS,color=Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.b, aes(ymin=BS-BSse,ymax=BS+BSse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.b, aes(x=as.numeric(TP), y=BS, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.b, aes(x=as.numeric(TP), y=BS, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 10) +
  theme(legend.position="none")

BS1b
BS1= aplot::plot_list(BS1b, BS1a, nrow = 1, ncol = 2, heights = c(5,5,5,3))
BS1


######Combine####

###TEST TO MAKE LEGEND ALONE###

TEST <- ggplot(Phenos24.a, aes(x=as.numeric(TP), y=BS, group=Cage.Pheno))+
  geom_point(aes(color=Cage.Treat,shape=Pheno.Treat),size=5)+
  geom_line(aes(color=Cage.Treat))+
  geom_linerange(aes(ymin=BS-BSse,ymax=BS+BSse, color=Cage.Treat), size=1)+
  ylab("Development Time (hrs)")+
  xlab("Sample Timepoints") +
  scale_x_continuous(expand = c(0,0),limits = c(.8, 5.5), breaks = c(3,5),labels=c("TP2","TP4"))+
  scale_x_break(c(1.1, 2.9), scales = 200)+
  scale_colour_manual(name = "Cage Treatment ", 
                      breaks = c("a", "b"),
                      labels = c("Low", "High"),
                      values = c("red","black", "green"), 
                      drop = FALSE)+
  scale_shape_manual(name = "Phenotyping Treatment ", 
                     breaks = c("a", "b"),
                     labels = c("Low","High"),
                     values = c(a=8,b=19),
                     drop = FALSE)+
  scale_linetype_manual(name = "Phenotyping Treatment ", 
                        breaks = c("a", "b"),
                        labels = c("Low","High"),
                        values=c(a="dashed",b="F1"))+
  geom_point(data=founder.a, aes(x=as.numeric(TP), y=BS,color = Cage.Treat, shape=Pheno.Treat), size=5)+
  geom_linerange(data=founder.a, aes(ymin=BS-BSse,ymax=BS+BSse, color=Cage.Treat), size=1)+
  geom_line(data=aa.pheno24.cage.a, aes(x=as.numeric(TP), y=BS, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) + 
  geom_line(data=bb.pheno24.cage.a, aes(x=as.numeric(TP), y=BS, group=as.factor(Cage), color=Cage.Treat.x, linetype=Pheno.Treat), 
            size=.5, alpha= .3) +
  theme_cowplot(font_size = 12) +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x =element_blank())


TEST


##get legend##

leg <- get_legend(TEST)

# Convert to a ggplot and print
as_ggplot(leg)

##FILLER to center LEGEND
s=plot_spacer() + theme_void()

##MULTIPLOT##
aplot::plot_list(LDm1, LDf1, SRM1, SRF1, Viab1, BS1, TOTF1,leg, nrow = 4, ncol = 2, heights = c(5,5,5,5))

Fig2SUPP <- ggarrange(LDm1, LDf1, SRM1, SRF1, Viab1, BS1, TOTF1,leg, 
                   labels = c("A", "B", "C", "D", "E", "F", "G", ""),
                   widths = c(1,1),
                   ncol = 2, nrow = 4)

Fig2SUPP

Fig2plastic <- ggarrange( LDf1, SRF1, Viab1, BS1, 
                      labels = c("I", "J", "K", "L"),
                      widths = c(1,1),
                      ncol = 2, nrow = 2)

Fig2plastic
Fig2
Fig2all<- ggarrange( Fig2plastic, Fig2,
                     labels = c("",""),
                     Heights = c(1,.7),
                     ncol = 1, nrow = 2)
Fig2all

my_plot_list <- list(Fig2,Fig2plastic)

Fig2all= ggarrange(plotlist = my_plot_list, heights = c(1,.8), nrow = 2)

Fig2all



## LINEAR MODELING FOR Reciprvol transfer SUPP ####
##THREE_WAY COMBINED LINEAR
View(Phenos.cage24)

Phenos.cage24$TP=as.numeric(Phenos.cage24$TP)
Phenos.cage24$Cage.Treat.x=as.factor(Phenos.cage24$Cage.Treat.x)
Phenos.cage24$Pheno.Treat=as.factor(Phenos.cage24$Pheno.Treat)
Phenos.cage24$Cage=as.factor(Phenos.cage24$Cage)

testabLDm1<-lmer(LDm ~ TP*Pheno.Treat*Cage.Treat.x + (1|Cage) + (1|TP), data = Phenos.cage24)

abLDm1<-lmer(LDm ~ TP*Pheno.Treat*Cage.Treat.x, + (1|Cage), data=Phenos.cage24) ##
anova(abLDm1) 
anova(testabLDm1)
summary(abLDm1)
summary(testabLDm1)##

abLDf1<-lme(LDf ~ TP*Pheno.Treat*Cage.Treat.x, random=~1|Cage/TP, data=Phenos.cage24) ##
anova(abLDf1) 
summary(abLDf1)

abV1<-lme(Viab ~ TP*Pheno.Treat*Cage.Treat.x, random=~1|Cage, data=Phenos.cage24) ##
anova(abV1) 
summary(abV1)

abSRM1<-lme(SRM ~ TP*Pheno.Treat*Cage.Treat.x, random=~1|Cage, data=Phenos.cage24) ##
anova(abSRM1) 
0.269summary(abSRM1)
39.59
abSRF1<-lme(SRF ~ TP*Pheno.Treat*Cage.Treat.x, random=~1|Cage, data=Phenos.cage24) ##
anova(abSRF1) 
summary(abSRF1)

abTOTF1<-lme(TOTF ~ TP*Pheno.Treat*Cage.Treat.x, random=~1|Cage, data=Phenos.cage24) ##
anova(abTOTF1) 
summary(abTOTF1)

abBS1<-lme(BS ~ TP*Pheno.Treat*Cage.Treat.x, random=~1|Cage/TP, data=Phenos.cage24) ##
anova(abBS1) 
0.26summary(abBS1)
1.

##LQ CAGES alone

a.Phenos.cage24=subset(Phenos.cage24, Phenos.cage24$Cage.Treat.x == "a")

a.Phenos.cage24$TP=as.factor(a.Phenos.cage24$TP)
a.Phenos.cage24$Pheno.Treat=as.factor(a.Phenos.cage24$Pheno.Treat)


aabLDm1<-lme(LDm ~ TP*Pheno.Treat, random=~1|Cage/TP, data=a.Phenos.cage24) ##
anova(aabLDm1) 
summary(aabLDm1)##

aabLDf1<-lme(LDf ~ TP*Pheno.Treat, random=~1|Cage/TP, data=a.Phenos.cage24) ##
anova(aabLDf1) 
summary(aabLDf1)

aabV1<-lme(Viab ~ TP*Pheno.Treat, random=~1|Cage/TP, data=a.Phenos.cage24) ##
anova(aabV1) 
summary(aabV1)

aabSRM1<-lme(SRM ~ TP*Pheno.Treat, random=~1|Cage/TP, data=a.Phenos.cage24) ##
anova(aabSRM1) 
summary(aabSRM1)

aabSRF1<-lme(SRF ~ TP*Pheno.Treat, random=~1|Cage/TP, data=a.Phenos.cage24) ##
anova(aabSRF1) 
summary(aabSRF1)

aabTOTF1<-lme(TOTF ~ TP*Pheno.Treat, random=~1|Cage/TP, data=a.Phenos.cage24) ##
anova(aabTOTF1) 
summary(aabTOTF1)

aabBS1<-lme(BS ~ TP*Pheno.Treat, random=~1|Cage/TP, data=a.Phenos.cage24) ##
anova(aabBS1) 
summary(aabBS1)




##HQ CAGES alone

b.Phenos.cage24=subset(Phenos.cage24, Phenos.cage24$Cage.Treat.x == "b")

b.Phenos.cage24$TP=as.factor(b.Phenos.cage24$TP)
b.Phenos.cage24$Pheno.Treat=as.factor(b.Phenos.cage24$Pheno.Treat)


babLDm1<-lme(LDm ~ TP*Pheno.Treat, random=~1|Cage/TP, data=b.Phenos.cage24) ##
anova(babLDm1) 
summary(babLDm1)##

babLDf1<-lme(LDf ~ TP*Pheno.Treat, random=~1|Cage/TP, data=b.Phenos.cage24) ##
anova(babLDf1) 
summary(babLDf1)

babV1<-lme(Viab ~ TP*Pheno.Treat, random=~1|Cage/TP, data=b.Phenos.cage24) ##
anova(babV1) 
summary(babV1)

babSRM1<-lme(SRM ~ TP*Pheno.Treat, random=~1|Cage/TP, data=b.Phenos.cage24) ##
anova(babSRM1) 
summary(babSRM1)

babSRF1<-lme(SRF ~ TP*Pheno.Treat, random=~1|Cage/TP, data=b.Phenos.cage24) ##
anova(babSRF1) 
summary(babSRF1)

babTOTF1<-lme(TOTF ~ TP*Pheno.Treat, random=~1|Cage/TP, data=b.Phenos.cage24) ##
anova(babTOTF1) 
summary(babTOTF1)

babBS1<-lme(BS ~ TP*Pheno.Treat, random=~1|Cage/TP, data=b.Phenos.cage24) ##
anova(babBS1) 
summary(babBS1)





