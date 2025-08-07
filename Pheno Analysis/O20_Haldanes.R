#### FFIGURE 1 ####

####Library Packages####
#install.packages("lubridate")
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(cowplot)
library(nlme)
library(lme4)
library(lubridate)
library(ggplot2)
library(reshape2)
library(cowplot)

####Import Files####
weather = read_csv("raw/O20_weather_raw.csv")
View(weather)

##separating date and time##

weather <- separate(weather,TimeStamp, into = c("Date", "Time"), sep = "^\\S*\\K\\s+")
weather$Date = format(as.Date(weather$Date, "%m/%d/%Y"), "20%y-%m-%d")

### create summary df####
#View(weather)
weather.sum=weather %>%
  group_by(Date) %>%
  dplyr::summarise(N=n(),
                   MaxAirTC=max(AirTC_Avg,na.rm = TRUE),
                   MinAIRTC=min(AirTC_Avg,na.rm =TRUE),
                   AirTC=mean(AirTC_Avg,na.rm = TRUE),
                   AirTCsd=sd(AirTC_Avg,na.rm = TRUE),
                   AirTCse= AirTCsd/sqrt(N),
                   DegreeDays = (((MaxAirTC + MinAIRTC)/2) - 10),
                   Generations = (DegreeDays / 180)
  )
wdata=unique(weather.sum)
wdata = dplyr::mutate(wdata, ExpDAY = row_number())
View(wdata)
#View(wdata)
#sum(wdata$Generations)
TPF1=subset(wdata, wdata$ExpDAY %in% (1:30))
sum(TPF1$Generations) ##3.150167
TP12=subset(wdata, wdata$ExpDAY %in% (30:55))
sum(TP12$Generations) ##2.284917
TP23=subset(wdata, wdata$ExpDAY %in% (56:76))
sum(TP23$Generations) ##1.267639
TP34=subset(wdata, wdata$ExpDAY %in% (78:118))
sum(TP34$Generations) ##1.300628
TP45=subset(wdata, wdata$ExpDAY %in% (118:133))
sum(TP45$Generations) ##0.1009917
TP15=subset(wdata, wdata$ExpDAY %in% (30:133))
sum(TP15$Generations) ##4.984981
TPF5=subset(wdata, wdata$ExpDAY %in% (1:133))
sum(TPF5$Generations) ##8.044814
###### 25c 7 day generation time = 180 degree days per generation = standardized works out to 8 genrations per season 

Phenos=read.csv("Phenos.csv")
#View(Pheno)
 ## convert all pheno means and sd to natural log
Phenos$logLD=log(Phenos$LD)
Phenos$logLDsd=log(Phenos$LDsd)
Phenos$logLDm=log(Phenos$LDm)
Phenos$logLDmsd=log(Phenos$LDmsd)
Phenos$logLDf=log(Phenos$LDf)
Phenos$logLDfsd=log(Phenos$LDfsd)
Phenos$logSRM=log(Phenos$SRM)
Phenos$logSRMsd=log(Phenos$SRMsd)
Phenos$logSRF=log(Phenos$SRF)
Phenos$logSRFsd=log(Phenos$SRFsd)
Phenos$logFEC=log(Phenos$Fec)
Phenos$logFECsd=log(Phenos$Fecsd)
Phenos$logVIAB=log(Phenos$Viab)
Phenos$logVIABsd=log(Phenos$Viabsd)
Phenos$logBS=log(Phenos$BS)
Phenos$logBSsd=log(Phenos$BSsd)

#View(Phenos)

#### Each of these creates a d term SP term for each phenotype, I term is determiend from degreee day calcultions
####A01####
A01HALD= Phenos %>% 
  filter( TP %in% c(0, 1)) %>% 
  filter( Pheno.Treat %in% c("a")) %>% 
  dplyr::summarise(Treatment= "Low",
                   Interval = "01",
                   LDd = (logLD[which(TP==1)]) - (logLD[which(TP == 0)]),
                   LDMd = (logLDm[which(TP==1)]) - (logLDm[which(TP == 0)]),
                   LDFd = (logLDf[which(TP==1)]) - (logLDf[which(TP == 0)]),
                   SRMd = (logSRM[which(TP==1)]) - (logSRM[which(TP == 0)]),
                   SRFd = (logSRF[which(TP==1)]) - (logSRF[which(TP == 0)]),
                   FECd = (logFEC[which(TP==1)]) - (logFEC[which(TP == 0)]),
                   VIABd = (logVIAB[which(TP==1)]) - (logVIAB[which(TP == 0)]),
                   BSd = (logBS[which(TP==1)]) - (logBS[which(TP == 0)]),
                   LDSp = sqrt((((logLD[which(TP==0)]-1)*(logLDsd[which(TP==0)]^2)) + ((logLD[which(TP==1)]-1)*(logLDsd[which(TP==1)]^2)))/
                                  ((logLD[which(TP==0)] + logLD[which(TP==1)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==0)]-1)*(logLDmsd[which(TP==0)]^2)) + ((logLDm[which(TP==1)]-1)*(logLDmsd[which(TP==1)]^2)))/
                                 ((logLDm[which(TP==0)] + logLDm[which(TP==1)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==0)]-1)*(logLDfsd[which(TP==0)]^2)) + ((logLDf[which(TP==1)]-1)*(logLDfsd[which(TP==1)]^2)))/
                                  ((logLDf[which(TP==0)] + logLDf[which(TP==1)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==0)]-1)*(logSRMsd[which(TP==0)]^2)) + ((logSRM[which(TP==1)]-1)*(logSRMsd[which(TP==1)]^2)))/
                                  ((logSRM[which(TP==0)] + logSRM[which(TP==1)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==0)]-1)*(logSRFsd[which(TP==0)]^2)) + ((logSRF[which(TP==1)]-1)*(logSRFsd[which(TP==1)]^2)))/
                                  ((logSRF[which(TP==0)] + logSRF[which(TP==1)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==0)]-1)*(logFECsd[which(TP==0)]^2)) + ((logFEC[which(TP==1)]-1)*(logFECsd[which(TP==1)]^2)))/
                                  ((logFEC[which(TP==0)] + logFEC[which(TP==1)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==0)]-1)*(logVIABsd[which(TP==0)]^2)) + ((logVIAB[which(TP==1)]-1)*(logVIABsd[which(TP==1)]^2)))/
                                  ((logVIAB[which(TP==0)] + logVIAB[which(TP==1)])-2)),
                   BSSp = sqrt((((logBS[which(TP==0)]-1)*(logBSsd[which(TP==0)]^2)) + ((logBS[which(TP==1)]-1)*(logBSsd[which(TP==1)]^2)))/
                                  ((logBS[which(TP==0)] + logBS[which(TP==1)])-2)),
                   I=3.150167,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
A01HALD <- A01HALD[ -c(3:19) ]                   
View(A01HALD)

####B01####
B01HALD= Phenos %>% 
  filter( TP %in% c(0, 1)) %>% 
  filter( Pheno.Treat %in% c("b")) %>% 
  dplyr::summarise(Treatment= "High",
                   Interval = "01",
                   LDd = (logLD[which(TP==1)]) - (logLD[which(TP == 0)]),
                   LDMd = (logLDm[which(TP==1)]) - (logLDm[which(TP == 0)]),
                   LDFd = (logLDf[which(TP==1)]) - (logLDf[which(TP == 0)]),
                   SRMd = (logSRM[which(TP==1)]) - (logSRM[which(TP == 0)]),
                   SRFd = (logSRF[which(TP==1)]) - (logSRF[which(TP == 0)]),
                   FECd = (logFEC[which(TP==1)]) - (logFEC[which(TP == 0)]),
                   VIABd = (logVIAB[which(TP==1)]) - (logVIAB[which(TP == 0)]),
                   BSd = (logBS[which(TP==1)]) - (logBS[which(TP == 0)]),
                   LDSp = sqrt((((logLD[which(TP==0)]-1)*(logLDsd[which(TP==0)]^2)) + ((logLD[which(TP==1)]-1)*(logLDsd[which(TP==1)]^2)))/
                                 ((logLD[which(TP==0)] + logLD[which(TP==1)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==0)]-1)*(logLDmsd[which(TP==0)]^2)) + ((logLDm[which(TP==1)]-1)*(logLDmsd[which(TP==1)]^2)))/
                                  ((logLDm[which(TP==0)] + logLDm[which(TP==1)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==0)]-1)*(logLDfsd[which(TP==0)]^2)) + ((logLDf[which(TP==1)]-1)*(logLDfsd[which(TP==1)]^2)))/
                                  ((logLDf[which(TP==0)] + logLDf[which(TP==1)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==0)]-1)*(logSRMsd[which(TP==0)]^2)) + ((logSRM[which(TP==1)]-1)*(logSRMsd[which(TP==1)]^2)))/
                                  ((logSRM[which(TP==0)] + logSRM[which(TP==1)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==0)]-1)*(logSRFsd[which(TP==0)]^2)) + ((logSRF[which(TP==1)]-1)*(logSRFsd[which(TP==1)]^2)))/
                                  ((logSRF[which(TP==0)] + logSRF[which(TP==1)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==0)]-1)*(logFECsd[which(TP==0)]^2)) + ((logFEC[which(TP==1)]-1)*(logFECsd[which(TP==1)]^2)))/
                                  ((logFEC[which(TP==0)] + logFEC[which(TP==1)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==0)]-1)*(logVIABsd[which(TP==0)]^2)) + ((logVIAB[which(TP==1)]-1)*(logVIABsd[which(TP==1)]^2)))/
                                   ((logVIAB[which(TP==0)] + logVIAB[which(TP==1)])-2)),
                   BSSp = sqrt((((logBS[which(TP==0)]-1)*(logBSsd[which(TP==0)]^2)) + ((logBS[which(TP==1)]-1)*(logBSsd[which(TP==1)]^2)))/
                                 ((logBS[which(TP==0)] + logBS[which(TP==1)])-2)),
                   I=3.150167,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
B01HALD <- B01HALD[ -c(3:19) ]                   

####A12####
A12HALD= Phenos %>% 
  filter( TP %in% c(1, 2)) %>% 
  filter( Cage.Pheno %in% c("aa")) %>% 
  dplyr::summarise(Treatment= "Low",
                   Interval = "12",
                   LDd = (logLD[which(TP==2)]) - (logLD[which(TP ==1)]),
                   LDMd = (logLDm[which(TP==2)]) - (logLDm[which(TP ==1)]),
                   LDFd = (logLDf[which(TP==2)]) - (logLDf[which(TP ==1)]),
                   SRMd = (logSRM[which(TP==2)]) - (logSRM[which(TP ==1)]),
                   SRFd = (logSRF[which(TP==2)]) - (logSRF[which(TP ==1)]),
                   FECd = (logFEC[which(TP==2)]) - (logFEC[which(TP ==1)]),
                   VIABd = (logVIAB[which(TP==2)]) - (logVIAB[which(TP ==1)]),
                   BSd = (logBS[which(TP==2)]) - (logBS[which(TP ==1)]),
                   LDSp = sqrt((((logLD[which(TP==1)]-1)*(logLDsd[which(TP==1)]^2)) + ((logLD[which(TP==2)]-1)*(logLDsd[which(TP==2)]^2)))/
                                 ((logLD[which(TP==1)] + logLD[which(TP==2)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==1)]-1)*(logLDmsd[which(TP==1)]^2)) + ((logLDm[which(TP==2)]-1)*(logLDmsd[which(TP==2)]^2)))/
                                  ((logLDm[which(TP==1)] + logLDm[which(TP==2)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==1)]-1)*(logLDfsd[which(TP==1)]^2)) + ((logLDf[which(TP==2)]-1)*(logLDfsd[which(TP==2)]^2)))/
                                  ((logLDf[which(TP==1)] + logLDf[which(TP==2)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==1)]-1)*(logSRMsd[which(TP==1)]^2)) + ((logSRM[which(TP==2)]-1)*(logSRMsd[which(TP==2)]^2)))/
                                  ((logSRM[which(TP==1)] + logSRM[which(TP==2)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==1)]-1)*(logSRFsd[which(TP==1)]^2)) + ((logSRF[which(TP==2)]-1)*(logSRFsd[which(TP==2)]^2)))/
                                  ((logSRF[which(TP==1)] + logSRF[which(TP==2)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==1)]-1)*(logFECsd[which(TP==1)]^2)) + ((logFEC[which(TP==2)]-1)*(logFECsd[which(TP==2)]^2)))/
                                  ((logFEC[which(TP==1)] + logFEC[which(TP==2)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==1)]-1)*(logVIABsd[which(TP==1)]^2)) + ((logVIAB[which(TP==2)]-1)*(logVIABsd[which(TP==2)]^2)))/
                                   ((logVIAB[which(TP==1)] + logVIAB[which(TP==2)])-2)),
                   BSSp = sqrt((((logBS[which(TP==1)]-1)*(logBSsd[which(TP==1)]^2)) + ((logBS[which(TP==2)]-1)*(logBSsd[which(TP==2)]^2)))/
                                 ((logBS[which(TP==1)] + logBS[which(TP==2)])-2)),
                   I=2.284917,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
A12HALD <- A12HALD[ -c(3:19) ]                   
View(A12HALD)

####B12####
B12HALD= Phenos %>% 
  filter( TP %in% c(1, 2)) %>% 
  filter( Cage.Pheno %in% c("bb")) %>% 
  dplyr::summarise(Treatment= "High",
                   Interval = "12",
                   LDd = (logLD[which(TP==2)]) - (logLD[which(TP ==1)]),
                   LDMd = (logLDm[which(TP==2)]) - (logLDm[which(TP ==1)]),
                   LDFd = (logLDf[which(TP==2)]) - (logLDf[which(TP ==1)]),
                   SRMd = (logSRM[which(TP==2)]) - (logSRM[which(TP ==1)]),
                   SRFd = (logSRF[which(TP==2)]) - (logSRF[which(TP ==1)]),
                   FECd = (logFEC[which(TP==2)]) - (logFEC[which(TP ==1)]),
                   VIABd = (logVIAB[which(TP==2)]) - (logVIAB[which(TP ==1)]),
                   BSd = (logBS[which(TP==2)]) - (logBS[which(TP ==1)]),
                   LDSp = sqrt((((logLD[which(TP==1)]-1)*(logLDsd[which(TP==1)]^2)) + ((logLD[which(TP==2)]-1)*(logLDsd[which(TP==2)]^2)))/
                                 ((logLD[which(TP==1)] + logLD[which(TP==2)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==1)]-1)*(logLDmsd[which(TP==1)]^2)) + ((logLDm[which(TP==2)]-1)*(logLDmsd[which(TP==2)]^2)))/
                                  ((logLDm[which(TP==1)] + logLDm[which(TP==2)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==1)]-1)*(logLDfsd[which(TP==1)]^2)) + ((logLDf[which(TP==2)]-1)*(logLDfsd[which(TP==2)]^2)))/
                                  ((logLDf[which(TP==1)] + logLDf[which(TP==2)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==1)]-1)*(logSRMsd[which(TP==1)]^2)) + ((logSRM[which(TP==2)]-1)*(logSRMsd[which(TP==2)]^2)))/
                                  ((logSRM[which(TP==1)] + logSRM[which(TP==2)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==1)]-1)*(logSRFsd[which(TP==1)]^2)) + ((logSRF[which(TP==2)]-1)*(logSRFsd[which(TP==2)]^2)))/
                                  ((logSRF[which(TP==1)] + logSRF[which(TP==2)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==1)]-1)*(logFECsd[which(TP==1)]^2)) + ((logFEC[which(TP==2)]-1)*(logFECsd[which(TP==2)]^2)))/
                                  ((logFEC[which(TP==1)] + logFEC[which(TP==2)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==1)]-1)*(logVIABsd[which(TP==1)]^2)) + ((logVIAB[which(TP==2)]-1)*(logVIABsd[which(TP==2)]^2)))/
                                   ((logVIAB[which(TP==1)] + logVIAB[which(TP==2)])-2)),
                   BSSp = sqrt((((logBS[which(TP==1)]-1)*(logBSsd[which(TP==1)]^2)) + ((logBS[which(TP==2)]-1)*(logBSsd[which(TP==2)]^2)))/
                                 ((logBS[which(TP==1)] + logBS[which(TP==2)])-2)),
                   I=2.284917,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
B12HALD <- B12HALD[ -c(3:19) ]                   
#View(B12HALD)

####A23####
A23HALD= Phenos %>% 
  filter( TP %in% c(2, 3)) %>% 
  filter( Cage.Pheno %in% c("aa")) %>% 
  dplyr::summarise(Treatment= "Low",
                   Interval = "23",
                   LDd = (logLD[which(TP==3)]) - (logLD[which(TP==2)]),
                   LDMd = (logLDm[which(TP==3)]) - (logLDm[which(TP==2)]),
                   LDFd = (logLDf[which(TP==3)]) - (logLDf[which(TP==2)]),
                   SRMd = (logSRM[which(TP==3)]) - (logSRM[which(TP==2)]),
                   SRFd = (logSRF[which(TP==3)]) - (logSRF[which(TP==2)]),
                   FECd = (logFEC[which(TP==3)]) - (logFEC[which(TP==2)]),
                   VIABd = (logVIAB[which(TP==3)]) - (logVIAB[which(TP==2)]),
                   BSd = (logBS[which(TP==3)]) - (logBS[which(TP==2)]),
                   LDSp = sqrt((((logLD[which(TP==2)]-1)*(logLDsd[which(TP==2)]^2)) + ((logLD[which(TP==3)]-1)*(logLDsd[which(TP==3)]^2)))/
                                 ((logLD[which(TP==2)] + logLD[which(TP==3)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==2)]-1)*(logLDmsd[which(TP==2)]^2)) + ((logLDm[which(TP==3)]-1)*(logLDmsd[which(TP==3)]^2)))/
                                  ((logLDm[which(TP==2)] + logLDm[which(TP==3)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==2)]-1)*(logLDfsd[which(TP==2)]^2)) + ((logLDf[which(TP==3)]-1)*(logLDfsd[which(TP==3)]^2)))/
                                  ((logLDf[which(TP==2)] + logLDf[which(TP==3)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==2)]-1)*(logSRMsd[which(TP==2)]^2)) + ((logSRM[which(TP==3)]-1)*(logSRMsd[which(TP==3)]^2)))/
                                  ((logSRM[which(TP==2)] + logSRM[which(TP==3)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==2)]-1)*(logSRFsd[which(TP==2)]^2)) + ((logSRF[which(TP==3)]-1)*(logSRFsd[which(TP==3)]^2)))/
                                  ((logSRF[which(TP==2)] + logSRF[which(TP==3)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==2)]-1)*(logFECsd[which(TP==2)]^2)) + ((logFEC[which(TP==3)]-1)*(logFECsd[which(TP==3)]^2)))/
                                  ((logFEC[which(TP==2)] + logFEC[which(TP==3)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==2)]-1)*(logVIABsd[which(TP==2)]^2)) + ((logVIAB[which(TP==3)]-1)*(logVIABsd[which(TP==3)]^2)))/
                                   ((logVIAB[which(TP==2)] + logVIAB[which(TP==3)])-2)),
                   BSSp = sqrt((((logBS[which(TP==2)]-1)*(logBSsd[which(TP==2)]^2)) + ((logBS[which(TP==3)]-1)*(logBSsd[which(TP==3)]^2)))/
                                 ((logBS[which(TP==2)] + logBS[which(TP==3)])-2)),
                   I=1.267639,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
A23HALD <- A23HALD[ -c(3:19) ]                   
#View(A23HALD)

####B23####
B23HALD= Phenos %>% 
  filter( TP %in% c(2, 3)) %>% 
  filter( Cage.Pheno %in% c("bb")) %>% 
  dplyr::summarise(Treatment= "High",
                   Interval = "23",
                   LDd = (logLD[which(TP==3)]) - (logLD[which(TP==2)]),
                   LDMd = (logLDm[which(TP==3)]) - (logLDm[which(TP==2)]),
                   LDFd = (logLDf[which(TP==3)]) - (logLDf[which(TP==2)]),
                   SRMd = (logSRM[which(TP==3)]) - (logSRM[which(TP==2)]),
                   SRFd = (logSRF[which(TP==3)]) - (logSRF[which(TP==2)]),
                   FECd = (logFEC[which(TP==3)]) - (logFEC[which(TP==2)]),
                   VIABd = (logVIAB[which(TP==3)]) - (logVIAB[which(TP==2)]),
                   BSd = (logBS[which(TP==3)]) - (logBS[which(TP==2)]),
                   LDSp = sqrt((((logLD[which(TP==2)]-1)*(logLDsd[which(TP==2)]^2)) + ((logLD[which(TP==3)]-1)*(logLDsd[which(TP==3)]^2)))/
                                 ((logLD[which(TP==2)] + logLD[which(TP==3)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==2)]-1)*(logLDmsd[which(TP==2)]^2)) + ((logLDm[which(TP==3)]-1)*(logLDmsd[which(TP==3)]^2)))/
                                  ((logLDm[which(TP==2)] + logLDm[which(TP==3)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==2)]-1)*(logLDfsd[which(TP==2)]^2)) + ((logLDf[which(TP==3)]-1)*(logLDfsd[which(TP==3)]^2)))/
                                  ((logLDf[which(TP==2)] + logLDf[which(TP==3)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==2)]-1)*(logSRMsd[which(TP==2)]^2)) + ((logSRM[which(TP==3)]-1)*(logSRMsd[which(TP==3)]^2)))/
                                  ((logSRM[which(TP==2)] + logSRM[which(TP==3)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==2)]-1)*(logSRFsd[which(TP==2)]^2)) + ((logSRF[which(TP==3)]-1)*(logSRFsd[which(TP==3)]^2)))/
                                  ((logSRF[which(TP==2)] + logSRF[which(TP==3)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==2)]-1)*(logFECsd[which(TP==2)]^2)) + ((logFEC[which(TP==3)]-1)*(logFECsd[which(TP==3)]^2)))/
                                  ((logFEC[which(TP==2)] + logFEC[which(TP==3)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==2)]-1)*(logVIABsd[which(TP==2)]^2)) + ((logVIAB[which(TP==3)]-1)*(logVIABsd[which(TP==3)]^2)))/
                                   ((logVIAB[which(TP==2)] + logVIAB[which(TP==3)])-2)),
                   BSSp = sqrt((((logBS[which(TP==2)]-1)*(logBSsd[which(TP==2)]^2)) + ((logBS[which(TP==3)]-1)*(logBSsd[which(TP==3)]^2)))/
                                 ((logBS[which(TP==2)] + logBS[which(TP==3)])-2)),
                   I=1.267639,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
B23HALD <- B23HALD[ -c(3:19) ]                   
#View(B23HALD)

####A34####
A34HALD= Phenos %>% 
  filter( TP %in% c(3, 4)) %>% 
  filter( Cage.Pheno %in% c("aa")) %>% 
  dplyr::summarise(Treatment= "Low",
                   Interval = "34",
                   LDd = (logLD[which(TP==4)]) - (logLD[which(TP==3)]),
                   LDMd = (logLDm[which(TP==4)]) - (logLDm[which(TP==3)]),
                   LDFd = (logLDf[which(TP==4)]) - (logLDf[which(TP==3)]),
                   SRMd = (logSRM[which(TP==4)]) - (logSRM[which(TP==3)]),
                   SRFd = (logSRF[which(TP==4)]) - (logSRF[which(TP==3)]),
                   FECd = (logFEC[which(TP==4)]) - (logFEC[which(TP==3)]),
                   VIABd = (logVIAB[which(TP==4)]) - (logVIAB[which(TP==3)]),
                   BSd = (logBS[which(TP==4)]) - (logBS[which(TP==3)]),
                   LDSp = sqrt((((logLD[which(TP==3)]-1)*(logLDsd[which(TP==3)]^2)) + ((logLD[which(TP==4)]-1)*(logLDsd[which(TP==4)]^2)))/
                                 ((logLD[which(TP==3)] + logLD[which(TP==4)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==3)]-1)*(logLDmsd[which(TP==3)]^2)) + ((logLDm[which(TP==4)]-1)*(logLDmsd[which(TP==4)]^2)))/
                                  ((logLDm[which(TP==3)] + logLDm[which(TP==4)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==3)]-1)*(logLDfsd[which(TP==3)]^2)) + ((logLDf[which(TP==4)]-1)*(logLDfsd[which(TP==4)]^2)))/
                                  ((logLDf[which(TP==3)] + logLDf[which(TP==4)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==3)]-1)*(logSRMsd[which(TP==3)]^2)) + ((logSRM[which(TP==4)]-1)*(logSRMsd[which(TP==4)]^2)))/
                                  ((logSRM[which(TP==3)] + logSRM[which(TP==4)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==3)]-1)*(logSRFsd[which(TP==3)]^2)) + ((logSRF[which(TP==4)]-1)*(logSRFsd[which(TP==4)]^2)))/
                                  ((logSRF[which(TP==3)] + logSRF[which(TP==4)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==3)]-1)*(logFECsd[which(TP==3)]^2)) + ((logFEC[which(TP==4)]-1)*(logFECsd[which(TP==4)]^2)))/
                                  ((logFEC[which(TP==3)] + logFEC[which(TP==4)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==3)]-1)*(logVIABsd[which(TP==3)]^2)) + ((logVIAB[which(TP==4)]-1)*(logVIABsd[which(TP==4)]^2)))/
                                   ((logVIAB[which(TP==3)] + logVIAB[which(TP==4)])-2)),
                   BSSp = sqrt((((logBS[which(TP==3)]-1)*(logBSsd[which(TP==3)]^2)) + ((logBS[which(TP==4)]-1)*(logBSsd[which(TP==4)]^2)))/
                                 ((logBS[which(TP==3)] + logBS[which(TP==4)])-2)),
                   I=1.300628,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
A34HALD <- A34HALD[ -c(3:19) ]                   
#View(A34HALD)

####B34####
B34HALD= Phenos %>% 
  filter( TP %in% c(3, 4)) %>% 
  filter( Cage.Pheno %in% c("bb")) %>% 
  dplyr::summarise(Treatment= "High",
                   Interval = "34",
                   LDd = (logLD[which(TP==4)]) - (logLD[which(TP==3)]),
                   LDMd = (logLDm[which(TP==4)]) - (logLDm[which(TP==3)]),
                   LDFd = (logLDf[which(TP==4)]) - (logLDf[which(TP==3)]),
                   SRMd = (logSRM[which(TP==4)]) - (logSRM[which(TP==3)]),
                   SRFd = (logSRF[which(TP==4)]) - (logSRF[which(TP==3)]),
                   FECd = (logFEC[which(TP==4)]) - (logFEC[which(TP==3)]),
                   VIABd = (logVIAB[which(TP==4)]) - (logVIAB[which(TP==3)]),
                   BSd = (logBS[which(TP==4)]) - (logBS[which(TP==3)]),
                   LDSp = sqrt((((logLD[which(TP==3)]-1)*(logLDsd[which(TP==3)]^2)) + ((logLD[which(TP==4)]-1)*(logLDsd[which(TP==4)]^2)))/
                                 ((logLD[which(TP==3)] + logLD[which(TP==4)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==3)]-1)*(logLDmsd[which(TP==3)]^2)) + ((logLDm[which(TP==4)]-1)*(logLDmsd[which(TP==4)]^2)))/
                                  ((logLDm[which(TP==3)] + logLDm[which(TP==4)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==3)]-1)*(logLDfsd[which(TP==3)]^2)) + ((logLDf[which(TP==4)]-1)*(logLDfsd[which(TP==4)]^2)))/
                                  ((logLDf[which(TP==3)] + logLDf[which(TP==4)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==3)]-1)*(logSRMsd[which(TP==3)]^2)) + ((logSRM[which(TP==4)]-1)*(logSRMsd[which(TP==4)]^2)))/
                                  ((logSRM[which(TP==3)] + logSRM[which(TP==4)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==3)]-1)*(logSRFsd[which(TP==3)]^2)) + ((logSRF[which(TP==4)]-1)*(logSRFsd[which(TP==4)]^2)))/
                                  ((logSRF[which(TP==3)] + logSRF[which(TP==4)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==3)]-1)*(logFECsd[which(TP==3)]^2)) + ((logFEC[which(TP==4)]-1)*(logFECsd[which(TP==4)]^2)))/
                                  ((logFEC[which(TP==3)] + logFEC[which(TP==4)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==3)]-1)*(logVIABsd[which(TP==3)]^2)) + ((logVIAB[which(TP==4)]-1)*(logVIABsd[which(TP==4)]^2)))/
                                   ((logVIAB[which(TP==3)] + logVIAB[which(TP==4)])-2)),
                   BSSp = sqrt((((logBS[which(TP==3)]-1)*(logBSsd[which(TP==3)]^2)) + ((logBS[which(TP==4)]-1)*(logBSsd[which(TP==4)]^2)))/
                                 ((logBS[which(TP==3)] + logBS[which(TP==4)])-2)),
                   I=1.300628,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
B34HALD <- B34HALD[ -c(3:19) ]                   
#View(B34HALD)

####A45####
A45HALD= Phenos %>% 
  filter( TP %in% c(4, 5)) %>% 
  filter( Cage.Pheno %in% c("aa")) %>% 
  dplyr::summarise(Treatment= "Low",
                   Interval = "45",
                   LDd = (logLD[which(TP==5)]) - (logLD[which(TP==4)]),
                   LDMd = (logLDm[which(TP==5)]) - (logLDm[which(TP==4)]),
                   LDFd = (logLDf[which(TP==5)]) - (logLDf[which(TP==4)]),
                   SRMd = (logSRM[which(TP==5)]) - (logSRM[which(TP==4)]),
                   SRFd = (logSRF[which(TP==5)]) - (logSRF[which(TP==4)]),
                   FECd = (logFEC[which(TP==5)]) - (logFEC[which(TP==4)]),
                   VIABd = (logVIAB[which(TP==5)]) - (logVIAB[which(TP==4)]),
                   BSd = (logBS[which(TP==5)]) - (logBS[which(TP==4)]),
                   LDSp = sqrt((((logLD[which(TP==4)]-1)*(logLDsd[which(TP==4)]^2)) + ((logLD[which(TP==5)]-1)*(logLDsd[which(TP==5)]^2)))/
                                 ((logLD[which(TP==4)] + logLD[which(TP==5)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==4)]-1)*(logLDmsd[which(TP==4)]^2)) + ((logLDm[which(TP==5)]-1)*(logLDmsd[which(TP==5)]^2)))/
                                  ((logLDm[which(TP==4)] + logLDm[which(TP==5)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==4)]-1)*(logLDfsd[which(TP==4)]^2)) + ((logLDf[which(TP==5)]-1)*(logLDfsd[which(TP==5)]^2)))/
                                  ((logLDf[which(TP==4)] + logLDf[which(TP==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==4)]-1)*(logSRMsd[which(TP==4)]^2)) + ((logSRM[which(TP==5)]-1)*(logSRMsd[which(TP==5)]^2)))/
                                  ((logSRM[which(TP==4)] + logSRM[which(TP==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==4)]-1)*(logSRFsd[which(TP==4)]^2)) + ((logSRF[which(TP==5)]-1)*(logSRFsd[which(TP==5)]^2)))/
                                  ((logSRF[which(TP==4)] + logSRF[which(TP==5)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==4)]-1)*(logFECsd[which(TP==4)]^2)) + ((logFEC[which(TP==5)]-1)*(logFECsd[which(TP==5)]^2)))/
                                  ((logFEC[which(TP==4)] + logFEC[which(TP==5)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==4)]-1)*(logVIABsd[which(TP==4)]^2)) + ((logVIAB[which(TP==5)]-1)*(logVIABsd[which(TP==5)]^2)))/
                                   ((logVIAB[which(TP==4)] + logVIAB[which(TP==5)])-2)),
                   BSSp = sqrt((((logBS[which(TP==4)]-1)*(logBSsd[which(TP==4)]^2)) + ((logBS[which(TP==5)]-1)*(logBSsd[which(TP==5)]^2)))/
                                 ((logBS[which(TP==4)] + logBS[which(TP==5)])-2)),
                   I=0.1009917,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
A45HALD <- A45HALD[ -c(3:19) ]                   


####B45####
B45HALD= Phenos %>% 
  filter( TP %in% c(4, 5)) %>% 
  filter( Cage.Pheno %in% c("bb")) %>% 
  dplyr::summarise(Treatment= "High",
                   Interval = "45",
                   LDd = (logLD[which(TP==5)]) - (logLD[which(TP==4)]),
                   LDMd = (logLDm[which(TP==5)]) - (logLDm[which(TP==4)]),
                   LDFd = (logLDf[which(TP==5)]) - (logLDf[which(TP==4)]),
                   SRMd = (logSRM[which(TP==5)]) - (logSRM[which(TP==4)]),
                   SRFd = (logSRF[which(TP==5)]) - (logSRF[which(TP==4)]),
                   FECd = (logFEC[which(TP==5)]) - (logFEC[which(TP==4)]),
                   VIABd = (logVIAB[which(TP==5)]) - (logVIAB[which(TP==4)]),
                   BSd = (logBS[which(TP==5)]) - (logBS[which(TP==4)]),
                   LDSp = sqrt((((logLD[which(TP==4)]-1)*(logLDsd[which(TP==4)]^2)) + ((logLD[which(TP==5)]-1)*(logLDsd[which(TP==5)]^2)))/
                                 ((logLD[which(TP==4)] + logLD[which(TP==5)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==4)]-1)*(logLDmsd[which(TP==4)]^2)) + ((logLDm[which(TP==5)]-1)*(logLDmsd[which(TP==5)]^2)))/
                                  ((logLDm[which(TP==4)] + logLDm[which(TP==5)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==4)]-1)*(logLDfsd[which(TP==4)]^2)) + ((logLDf[which(TP==5)]-1)*(logLDfsd[which(TP==5)]^2)))/
                                  ((logLDf[which(TP==4)] + logLDf[which(TP==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==4)]-1)*(logSRMsd[which(TP==4)]^2)) + ((logSRM[which(TP==5)]-1)*(logSRMsd[which(TP==5)]^2)))/
                                  ((logSRM[which(TP==4)] + logSRM[which(TP==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==4)]-1)*(logSRFsd[which(TP==4)]^2)) + ((logSRF[which(TP==5)]-1)*(logSRFsd[which(TP==5)]^2)))/
                                  ((logSRF[which(TP==4)] + logSRF[which(TP==5)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==4)]-1)*(logFECsd[which(TP==4)]^2)) + ((logFEC[which(TP==5)]-1)*(logFECsd[which(TP==5)]^2)))/
                                  ((logFEC[which(TP==4)] + logFEC[which(TP==5)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==4)]-1)*(logVIABsd[which(TP==4)]^2)) + ((logVIAB[which(TP==5)]-1)*(logVIABsd[which(TP==5)]^2)))/
                                   ((logVIAB[which(TP==4)] + logVIAB[which(TP==5)])-2)),
                   BSSp = sqrt((((logBS[which(TP==4)]-1)*(logBSsd[which(TP==4)]^2)) + ((logBS[which(TP==5)]-1)*(logBSsd[which(TP==5)]^2)))/
                                 ((logBS[which(TP==4)] + logBS[which(TP==5)])-2)),
                   I=0.1009917,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
B45HALD <- B45HALD[ -c(3:19) ]                   
#View(B45HALD)

####A15####
A15HALD= Phenos %>% 
  filter( TP %in% c(1, 5)) %>% 
  filter( Cage.Pheno %in% c("aa")) %>% 
  dplyr::summarise(Treatment= "Low",
                   Interval = "15",
                   LDd = (logLD[which(TP==5)]) - (logLD[which(TP==1)]),
                   LDMd = (logLDm[which(TP==5)]) - (logLDm[which(TP==1)]),
                   LDFd = (logLDf[which(TP==5)]) - (logLDf[which(TP==1)]),
                   SRMd = (logSRM[which(TP==5)]) - (logSRM[which(TP==1)]),
                   SRFd = (logSRF[which(TP==5)]) - (logSRF[which(TP==1)]),
                   FECd = (logFEC[which(TP==5)]) - (logFEC[which(TP==1)]),
                   VIABd = (logVIAB[which(TP==5)]) - (logVIAB[which(TP==1)]),
                   BSd = (logBS[which(TP==5)]) - (logBS[which(TP==1)]),
                   LDSp = sqrt((((logLD[which(TP==1)]-1)*(logLDsd[which(TP==1)]^2)) + ((logLD[which(TP==5)]-1)*(logLDsd[which(TP==5)]^2)))/
                                 ((logLD[which(TP==1)] + logLD[which(TP==5)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==1)]-1)*(logLDmsd[which(TP==1)]^2)) + ((logLDm[which(TP==5)]-1)*(logLDmsd[which(TP==5)]^2)))/
                                  ((logLDm[which(TP==1)] + logLDm[which(TP==5)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==1)]-1)*(logLDfsd[which(TP==1)]^2)) + ((logLDf[which(TP==5)]-1)*(logLDfsd[which(TP==5)]^2)))/
                                  ((logLDf[which(TP==1)] + logLDf[which(TP==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==1)]-1)*(logSRMsd[which(TP==1)]^2)) + ((logSRM[which(TP==5)]-1)*(logSRMsd[which(TP==5)]^2)))/
                                  ((logSRM[which(TP==1)] + logSRM[which(TP==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==1)]-1)*(logSRFsd[which(TP==1)]^2)) + ((logSRF[which(TP==5)]-1)*(logSRFsd[which(TP==5)]^2)))/
                                  ((logSRF[which(TP==1)] + logSRF[which(TP==5)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==1)]-1)*(logFECsd[which(TP==1)]^2)) + ((logFEC[which(TP==5)]-1)*(logFECsd[which(TP==5)]^2)))/
                                  ((logFEC[which(TP==1)] + logFEC[which(TP==5)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==1)]-1)*(logVIABsd[which(TP==1)]^2)) + ((logVIAB[which(TP==5)]-1)*(logVIABsd[which(TP==5)]^2)))/
                                   ((logVIAB[which(TP==1)] + logVIAB[which(TP==5)])-2)),
                   BSSp = sqrt((((logBS[which(TP==1)]-1)*(logBSsd[which(TP==1)]^2)) + ((logBS[which(TP==5)]-1)*(logBSsd[which(TP==5)]^2)))/
                                 ((logBS[which(TP==1)] + logBS[which(TP==5)])-2)),
                   I=4.984981,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
A15HALD <- A15HALD[ -c(3:19) ]                   
#View(A15HALD)

####B15####
B15HALD= Phenos %>% 
  filter( TP %in% c(1, 5)) %>% 
  filter( Cage.Pheno %in% c("bb")) %>% 
  dplyr::summarise(Treatment= "High",
                   Interval = "15",
                   LDd = (logLD[which(TP==5)]) - (logLD[which(TP==1)]),
                   LDMd = (logLDm[which(TP==5)]) - (logLDm[which(TP==1)]),
                   LDFd = (logLDf[which(TP==5)]) - (logLDf[which(TP==1)]),
                   SRMd = (logSRM[which(TP==5)]) - (logSRM[which(TP==1)]),
                   SRFd = (logSRF[which(TP==5)]) - (logSRF[which(TP==1)]),
                   FECd = (logFEC[which(TP==5)]) - (logFEC[which(TP==1)]),
                   VIABd = (logVIAB[which(TP==5)]) - (logVIAB[which(TP==1)]),
                   BSd = (logBS[which(TP==5)]) - (logBS[which(TP==1)]),
                   LDSp = sqrt((((logLD[which(TP==1)]-1)*(logLDsd[which(TP==1)]^2)) + ((logLD[which(TP==5)]-1)*(logLDsd[which(TP==5)]^2)))/
                                 ((logLD[which(TP==1)] + logLD[which(TP==5)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==1)]-1)*(logLDmsd[which(TP==1)]^2)) + ((logLDm[which(TP==5)]-1)*(logLDmsd[which(TP==5)]^2)))/
                                  ((logLDm[which(TP==1)] + logLDm[which(TP==5)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==1)]-1)*(logLDfsd[which(TP==1)]^2)) + ((logLDf[which(TP==5)]-1)*(logLDfsd[which(TP==5)]^2)))/
                                  ((logLDf[which(TP==1)] + logLDf[which(TP==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==1)]-1)*(logSRMsd[which(TP==1)]^2)) + ((logSRM[which(TP==5)]-1)*(logSRMsd[which(TP==5)]^2)))/
                                  ((logSRM[which(TP==1)] + logSRM[which(TP==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==1)]-1)*(logSRFsd[which(TP==1)]^2)) + ((logSRF[which(TP==5)]-1)*(logSRFsd[which(TP==5)]^2)))/
                                  ((logSRF[which(TP==1)] + logSRF[which(TP==5)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==1)]-1)*(logFECsd[which(TP==1)]^2)) + ((logFEC[which(TP==5)]-1)*(logFECsd[which(TP==5)]^2)))/
                                  ((logFEC[which(TP==1)] + logFEC[which(TP==5)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==1)]-1)*(logVIABsd[which(TP==1)]^2)) + ((logVIAB[which(TP==5)]-1)*(logVIABsd[which(TP==5)]^2)))/
                                   ((logVIAB[which(TP==1)] + logVIAB[which(TP==5)])-2)),
                   BSSp = sqrt((((logBS[which(TP==1)]-1)*(logBSsd[which(TP==1)]^2)) + ((logBS[which(TP==5)]-1)*(logBSsd[which(TP==5)]^2)))/
                                 ((logBS[which(TP==1)] + logBS[which(TP==5)])-2)),
                   I=4.984981,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
B15HALD <- B15HALD[ -c(3:19) ]                   
View(B15HALD)

####A05####
A05HALD= Phenos %>% 
  filter( TP %in% c(0, 5)) %>% 
  filter( Pheno.Treat %in% c("a")) %>% 
  dplyr::summarise(Treatment= "Low",
                   Interval = "05",
                   LDd = (logLD[which(TP==5)]) - (logLD[which(TP==0)]),
                   LDMd = (logLDm[which(TP==5)]) - (logLDm[which(TP==0)]),
                   LDFd = (logLDf[which(TP==5)]) - (logLDf[which(TP==0)]),
                   SRMd = (logSRM[which(TP==5)]) - (logSRM[which(TP==0)]),
                   SRFd = (logSRF[which(TP==5)]) - (logSRF[which(TP==0)]),
                   FECd = (logFEC[which(TP==5)]) - (logFEC[which(TP==0)]),
                   VIABd = (logVIAB[which(TP==5)]) - (logVIAB[which(TP==0)]),
                   BSd = (logBS[which(TP==5)]) - (logBS[which(TP==0)]),
                   LDSp = sqrt((((logLD[which(TP==0)]-1)*(logLDsd[which(TP==0)]^2)) + ((logLD[which(TP==5)]-1)*(logLDsd[which(TP==5)]^2)))/
                                 ((logLD[which(TP==0)] + logLD[which(TP==5)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==0)]-1)*(logLDmsd[which(TP==0)]^2)) + ((logLDm[which(TP==5)]-1)*(logLDmsd[which(TP==5)]^2)))/
                                  ((logLDm[which(TP==0)] + logLDm[which(TP==5)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==0)]-1)*(logLDfsd[which(TP==0)]^2)) + ((logLDf[which(TP==5)]-1)*(logLDfsd[which(TP==5)]^2)))/
                                  ((logLDf[which(TP==0)] + logLDf[which(TP==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==0)]-1)*(logSRMsd[which(TP==0)]^2)) + ((logSRM[which(TP==5)]-1)*(logSRMsd[which(TP==5)]^2)))/
                                  ((logSRM[which(TP==0)] + logSRM[which(TP==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==0)]-1)*(logSRFsd[which(TP==0)]^2)) + ((logSRF[which(TP==5)]-1)*(logSRFsd[which(TP==5)]^2)))/
                                  ((logSRF[which(TP==0)] + logSRF[which(TP==5)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==0)]-1)*(logFECsd[which(TP==0)]^2)) + ((logFEC[which(TP==5)]-1)*(logFECsd[which(TP==5)]^2)))/
                                  ((logFEC[which(TP==0)] + logFEC[which(TP==5)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==0)]-1)*(logVIABsd[which(TP==0)]^2)) + ((logVIAB[which(TP==5)]-1)*(logVIABsd[which(TP==5)]^2)))/
                                   ((logVIAB[which(TP==0)] + logVIAB[which(TP==5)])-2)),
                   BSSp = sqrt((((logBS[which(TP==0)]-1)*(logBSsd[which(TP==0)]^2)) + ((logBS[which(TP==5)]-1)*(logBSsd[which(TP==5)]^2)))/
                                 ((logBS[which(TP==0)] + logBS[which(TP==5)])-2)),
                   I=4.984981,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
A05HALD <- A05HALD[ -c(3:19) ]                   
View(A05HALD)

####B05####
B05HALD= Phenos %>% 
  filter( TP %in% c(0, 5)) %>% 
  filter( Pheno.Treat %in% c("b")) %>% 
  dplyr::summarise(Treatment= "High",
                   Interval = "05",
                   LDd = (logLD[which(TP==5)]) - (logLD[which(TP==0)]),
                   LDMd = (logLDm[which(TP==5)]) - (logLDm[which(TP==0)]),
                   LDFd = (logLDf[which(TP==5)]) - (logLDf[which(TP==0)]),
                   SRMd = (logSRM[which(TP==5)]) - (logSRM[which(TP==0)]),
                   SRFd = (logSRF[which(TP==5)]) - (logSRF[which(TP==0)]),
                   FECd = (logFEC[which(TP==5)]) - (logFEC[which(TP==0)]),
                   VIABd = (logVIAB[which(TP==5)]) - (logVIAB[which(TP==0)]),
                   BSd = (logBS[which(TP==5)]) - (logBS[which(TP==0)]),
                   LDSp = sqrt((((logLD[which(TP==0)]-1)*(logLDsd[which(TP==0)]^2)) + ((logLD[which(TP==5)]-1)*(logLDsd[which(TP==5)]^2)))/
                                 ((logLD[which(TP==0)] + logLD[which(TP==5)])-2)),
                   LDMSp = sqrt((((logLDm[which(TP==0)]-1)*(logLDmsd[which(TP==0)]^2)) + ((logLDm[which(TP==5)]-1)*(logLDmsd[which(TP==5)]^2)))/
                                  ((logLDm[which(TP==0)] + logLDm[which(TP==5)])-2)),
                   LDFSp = sqrt((((logLDf[which(TP==0)]-1)*(logLDfsd[which(TP==0)]^2)) + ((logLDf[which(TP==5)]-1)*(logLDfsd[which(TP==5)]^2)))/
                                  ((logLDf[which(TP==0)] + logLDf[which(TP==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(TP==0)]-1)*(logSRMsd[which(TP==0)]^2)) + ((logSRM[which(TP==5)]-1)*(logSRMsd[which(TP==5)]^2)))/
                                  ((logSRM[which(TP==0)] + logSRM[which(TP==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(TP==0)]-1)*(logSRFsd[which(TP==0)]^2)) + ((logSRF[which(TP==5)]-1)*(logSRFsd[which(TP==5)]^2)))/
                                  ((logSRF[which(TP==0)] + logSRF[which(TP==5)])-2)),
                   FECSp = sqrt((((logFEC[which(TP==0)]-1)*(logFECsd[which(TP==0)]^2)) + ((logFEC[which(TP==5)]-1)*(logFECsd[which(TP==5)]^2)))/
                                  ((logFEC[which(TP==0)] + logFEC[which(TP==5)])-2)),
                   VIABSp = sqrt((((logVIAB[which(TP==0)]-1)*(logVIABsd[which(TP==0)]^2)) + ((logVIAB[which(TP==5)]-1)*(logVIABsd[which(TP==5)]^2)))/
                                   ((logVIAB[which(TP==0)] + logVIAB[which(TP==5)])-2)),
                   BSSp = sqrt((((logBS[which(TP==0)]-1)*(logBSsd[which(TP==0)]^2)) + ((logBS[which(TP==5)]-1)*(logBSsd[which(TP==5)]^2)))/
                                 ((logBS[which(TP==0)] + logBS[which(TP==5)])-2)),
                   I=4.984981,
                   LD_Hal = ((LDd/LDSp)/I),
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   BS_HAL = ((BSd/BSSp)/I))
B05HALD <- B05HALD[ -c(3:19) ]                   
View(B05HALD)

####Merge all dfs####
A_HALD <- rbind(A12HALD,A23HALD,A34HALD,A45HALD,A15HALD)
B_HALD <- rbind(B12HALD,B23HALD,B34HALD,B45HALD,B15HALD)
HALD <- rbind(A_HALD,B_HALD)
View(HALD)

### make them long###
A_HALD_long <- melt(A_HALD ,  id.vars = c('Treatment','Interval'), variable.name = 'Phenotypes')
A_HALD_long$absValue <- abs(A_HALD_long$value)
A_HALD_long$logValue <- log(A_HALD_long$absValue)

B_HALD_long <- melt(B_HALD ,  id.vars = c('Treatment','Interval'), variable.name = 'Phenotypes')
B_HALD_long$absValue <- abs(B_HALD_long$value)
B_HALD_long$logValue <- log(B_HALD_long$absValue)

HALD_long <- melt(HALD ,  id.vars = c('Treatment','Interval'), variable.name = 'Phenotypes')
HALD_long$absValue <- abs(HALD_long$value)
HALD_long$logValue <- log10(HALD_long$absValue)
HALD_long$logg10value <- log10(HALD_long$absValue)

HALD_long$Interval <- factor(HALD$Interval, levels=c( "12", "23", "34", "45", "15"))

View(HALD_long)
#create line plot for each column in data frame
HALD_LONG_SUMM= HALD_long %>% 
  group_by(Treatment, Interval, Phenotypes) %>%
  dplyr::summarise(N=n(),
                   log_mean =mean(logValue),
                   log_sd = sd(logValue)
  )  
              
View(HALD_LONG_SUMM)
View(HALD_long)
ggplot(HALD_long, aes(Interval, logValue)) +
  geom_point(aes(colour = Phenotypes, group = Phenotypes), size = 3) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+
  scale_colour_manual(breaks=c("LD_Hal","LDM_HAL","LDF_HAL","SRM_HAL","SRF_HAL","FEC_HAL","VIAB_HAL","BS_HAL"), 
                      labels =c("Dev Time", "Dev Time (male)", "Dev Time (female)", "Starv Time (male)", "Starv Time (female)", "Fecundity", "Viability", "Body Size"), 
                      values =c("red","coral4","coral", "darkgreen", "green", "gold","blue","grey"))+
  scale_x_discrete(breaks = c("01", "12", "23", "34", "45", "15", "05"), labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "1-5", "0-5") )+
  ylab("log(Haldanes)") +
  xlab("Timepoint Interval") +
  facet_wrap(~ Treatment) + 
  theme_cowplot()

ggplot(HALD_long, aes(Interval, logValue, fill= Phenotypes)) +
  geom_point(aes(colour = Phenotypes), size = 3) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+
  scale_colour_manual(breaks=c("LD_Hal","LDM_HAL","LDF_HAL","SRM_HAL","SRF_HAL","FEC_HAL","VIAB_HAL","BS_HAL"), 
                      labels =c("Dev Time", "Dev Time (male)", "Dev Time (female)", "Starv Time (male)", "Starv Time (female)", "Fecundity", "Viability", "Body Size"), 
                      values =c("red","coral4","coral", "darkgreen", "green", "gold","blue","grey"))+
  scale_x_discrete(breaks = c("01", "12", "23", "34", "45", "15", "05"), labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "1-5", "0-5") )+
  ylab("log(Haldanes)") +
  xlab("Timepoint Interval") +
  facet_wrap(~ Treatment) + 
  theme_cowplot()


breaks <- seq(0, 1.1 , by = 24)
subset(HALD_long,Interval != "15")

Hald_combined=ggplot(subset(HALD_long,Interval != "15"), aes(1,absValue)) +
  geom_jitter(aes(shape= Treatment, color=Interval), size=5, width = .003) + 
  xlim(.995,1.005)+
  labs(color = "                                                                                                                                       Timepoint Interval")+
  labs(shape = "               Cage Treatment")+
  scale_shape_manual(breaks=c("Low", "High"), label=c("Low Quality","High Quality"), 
                     values =c(19,17))+
  scale_color_manual(breaks = c( "12", "23", "34", "45"),
                     labels = c( "1-2", "2-3", "3-4", "4-5"),
                     values = c("#D8D8D8","#909090","#000000","red"))+
  scale_y_continuous(name="Rate of Phenotypic Change (Haldanes)", trans = log_trans(100), breaks = c(0.0001, 0.0100, 1.0000), labels = c("0.0001", "0.0100", "1.0000")) +
  theme_cowplot(14) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
Hald_combined


load("fig2.rdata") 
Fig2

Fig2H <- ggarrange(Fig2, Hald_combined,
                  labels = c("", "I"),
                  widths = c(1,.3),
                  ncol = 2, nrow = 1, common.legend = TRUE)
Fig2H

