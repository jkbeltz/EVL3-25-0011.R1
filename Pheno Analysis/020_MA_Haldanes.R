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

#View(wdata)
#sum(wdata$Generations)

TP24=subset(wdata, wdata$ExpDAY %in% (56:118))
sum(TP24$Generations) ##2.640433
TP45=subset(wdata, wdata$ExpDAY %in% (118:133))
sum(TP45$Generations) ##0.1009917
TP25=subset(wdata, wdata$ExpDAY %in% (56:133))
sum(TP25$Generations) ##2.700064

###### 25c 7 day generation time = 180 degree days per generation = standardized works out to 8 genrations per season 

Phenos=read.csv("MAPhenos.csv")
View(Phenos)
## convert all pheno means and sd to natural log
Phenos$logLDM=log(Phenos$LDM)
Phenos$logLDMsd=log(Phenos$LDMsd)
Phenos$logLDF=log(Phenos$LDF)
Phenos$logLDFsd=log(Phenos$LDFsd)
Phenos$logSRM=log(Phenos$SRM)
Phenos$logSRMsd=log(Phenos$SRMsd)
Phenos$logSRF=log(Phenos$SRF)
Phenos$logSRFsd=log(Phenos$SRFsd)
Phenos$logFec=log(Phenos$Fec)
Phenos$logFecsd=log(Phenos$Fecsd)
Phenos$logViab=log(Phenos$Viab)
Phenos$logViabsd=log(Phenos$Viabsd)
Phenos$logDW=log(Phenos$DW)
Phenos$logDWsd=log(Phenos$DWsd)
Phenos$logLW=log(Phenos$LW)
Phenos$logLWsd=log(Phenos$LWsd)



#View(Phenos)

#### Each of these creates a d term SP term for each phenotype, I term is determiend from degreee day calcultions
####C24####
View(Phenos)
C24HALD= Phenos %>% 
  filter( timepoint %in% c(2, 4)) %>% 
  filter( cage.treat2 %in% c("C")) %>% 
  dplyr::summarise(Treatment= "Control",
                   Interval = "24",
                   LDMd = (logLDM[which(timepoint==4)]) - (logLDM[which(timepoint==2)]),
                   LDFd = (logLDF[which(timepoint==4)]) - (logLDF[which(timepoint==2)]),
                   SRMd = (logSRM[which(timepoint==4)]) - (logSRM[which(timepoint==2)]),
                   SRFd = (logSRF[which(timepoint==4)]) - (logSRF[which(timepoint==2)]),
                   FECd = (logFec[which(timepoint==4)]) - (logFec[which(timepoint==2)]),
                   VIABd = (logViab[which(timepoint==4)]) - (logViab[which(timepoint==2)]),
                   DWd = (logDW[which(timepoint==4)]) - (logDW[which(timepoint==2)]),
                   LWd = (logLW[which(timepoint==4)]) - (logLW[which(timepoint==2)]),
                   LDMSp = sqrt((((logLDM[which(timepoint==2)]-1)*(logLDMsd[which(timepoint==2)]^2)) + ((logLDM[which(timepoint==4)]-1)*(logLDMsd[which(timepoint==4)]^2)))/
                                  ((logLDM[which(timepoint==2)] + logLDM[which(timepoint==4)])-2)),
                   LDFSp = sqrt((((logLDF[which(timepoint==2)]-1)*(logLDFsd[which(timepoint==2)]^2)) + ((logLDF[which(timepoint==4)]-1)*(logLDFsd[which(timepoint==4)]^2)))/
                                  ((logLDF[which(timepoint==2)] + logLDF[which(timepoint==4)])-2)),
                   SRMSp = sqrt((((logSRM[which(timepoint==2)]-1)*(logSRMsd[which(timepoint==2)]^2)) + ((logSRM[which(timepoint==4)]-1)*(logSRMsd[which(timepoint==4)]^2)))/
                                  ((logSRM[which(timepoint==2)] + logSRM[which(timepoint==4)])-2)),
                   SRFSp = sqrt((((logSRF[which(timepoint==2)]-1)*(logSRFsd[which(timepoint==2)]^2)) + ((logSRF[which(timepoint==4)]-1)*(logSRFsd[which(timepoint==4)]^2)))/
                                  ((logSRF[which(timepoint==2)] + logSRF[which(timepoint==4)])-2)),
                   FECSp = sqrt((((logFec[which(timepoint==2)]-1)*(logFecsd[which(timepoint==2)]^2)) + ((logFec[which(timepoint==4)]-1)*(logFecsd[which(timepoint==4)]^2)))/
                                  ((logFec[which(timepoint==2)] + logFec[which(timepoint==4)])-2)),
                   VIABSp = sqrt((((logViab[which(timepoint==2)]-1)*(logViabsd[which(timepoint==2)]^2)) + ((logViab[which(timepoint==4)]-1)*(logViabsd[which(timepoint==4)]^2)))/
                                   ((logViab[which(timepoint==2)] + logViab[which(timepoint==4)])-2)),
                   DWSp = sqrt((((logDW[which(timepoint==2)]-1)*(logDWsd[which(timepoint==2)]^2)) + ((logDW[which(timepoint==4)]-1)*(logDWsd[which(timepoint==4)]^2)))/
                                 ((logDW[which(timepoint==2)] + logDW[which(timepoint==4)])-2)),
                   LWSp = sqrt((((logLW[which(timepoint==2)]-1)*(logLWsd[which(timepoint==2)]^2)) + ((logLW[which(timepoint==4)]-1)*(logLWsd[which(timepoint==4)]^2)))/
                                 ((logLW[which(timepoint==2)] + logLW[which(timepoint==4)])-2)),
                   I=2.640433,
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   DW_HAL = ((DWd/DWSp)/I),
                   LW_HAL = ((LWd/LWSp)/I))

C24HALD<- C24HALD[ -c(3:19) ] 
View(C24HALD) 
                 
#View(A23HALD)


####AT24####
AT24HALD= Phenos %>% 
  filter( timepoint %in% c(2, 4)) %>% 
  filter( cage.treat2 %in% c("AT")) %>% 
  dplyr::summarise(Treatment= "+AT",
                   Interval = "24",
                   LDMd = (logLDM[which(timepoint==4)]) - (logLDM[which(timepoint==2)]),
                   LDFd = (logLDF[which(timepoint==4)]) - (logLDF[which(timepoint==2)]),
                   SRMd = (logSRM[which(timepoint==4)]) - (logSRM[which(timepoint==2)]),
                   SRFd = (logSRF[which(timepoint==4)]) - (logSRF[which(timepoint==2)]),
                   FECd = (logFec[which(timepoint==4)]) - (logFec[which(timepoint==2)]),
                   VIABd = (logViab[which(timepoint==4)]) - (logViab[which(timepoint==2)]),
                   DWd = (logDW[which(timepoint==4)]) - (logDW[which(timepoint==2)]),
                   LWd = (logLW[which(timepoint==4)]) - (logLW[which(timepoint==2)]),
                   LDMSp = sqrt((((logLDM[which(timepoint==2)]-1)*(logLDMsd[which(timepoint==2)]^2)) + ((logLDM[which(timepoint==4)]-1)*(logLDMsd[which(timepoint==4)]^2)))/
                                  ((logLDM[which(timepoint==2)] + logLDM[which(timepoint==4)])-2)),
                   LDFSp = sqrt((((logLDF[which(timepoint==2)]-1)*(logLDFsd[which(timepoint==2)]^2)) + ((logLDF[which(timepoint==4)]-1)*(logLDFsd[which(timepoint==4)]^2)))/
                                  ((logLDF[which(timepoint==2)] + logLDF[which(timepoint==4)])-2)),
                   SRMSp = sqrt((((logSRM[which(timepoint==2)]-1)*(logSRMsd[which(timepoint==2)]^2)) + ((logSRM[which(timepoint==4)]-1)*(logSRMsd[which(timepoint==4)]^2)))/
                                  ((logSRM[which(timepoint==2)] + logSRM[which(timepoint==4)])-2)),
                   SRFSp = sqrt((((logSRF[which(timepoint==2)]-1)*(logSRFsd[which(timepoint==2)]^2)) + ((logSRF[which(timepoint==4)]-1)*(logSRFsd[which(timepoint==4)]^2)))/
                                  ((logSRF[which(timepoint==2)] + logSRF[which(timepoint==4)])-2)),
                   FECSp = sqrt((((logFec[which(timepoint==2)]-1)*(logFecsd[which(timepoint==2)]^2)) + ((logFec[which(timepoint==4)]-1)*(logFecsd[which(timepoint==4)]^2)))/
                                  ((logFec[which(timepoint==2)] + logFec[which(timepoint==4)])-2)),
                   VIABSp = sqrt((((logViab[which(timepoint==2)]-1)*(logViabsd[which(timepoint==2)]^2)) + ((logViab[which(timepoint==4)]-1)*(logViabsd[which(timepoint==4)]^2)))/
                                   ((logViab[which(timepoint==2)] + logViab[which(timepoint==4)])-2)),
                   DWSp = sqrt((((logDW[which(timepoint==2)]-1)*(logDWsd[which(timepoint==2)]^2)) + ((logDW[which(timepoint==4)]-1)*(logDWsd[which(timepoint==4)]^2)))/
                                 ((logDW[which(timepoint==2)] + logDW[which(timepoint==4)])-2)),
                   LWSp = sqrt((((logLW[which(timepoint==2)]-1)*(logLWsd[which(timepoint==2)]^2)) + ((logLW[which(timepoint==4)]-1)*(logLWsd[which(timepoint==4)]^2)))/
                                 ((logLW[which(timepoint==2)] + logLW[which(timepoint==4)])-2)),
                   I=2.640433,
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   DW_HAL = ((DWd/DWSp)/I),
                   LW_HAL = ((LWd/LWSp)/I))

AT24HALD<- AT24HALD[ -c(3:19) ]                   
View(AT24HALD) 
                  


####LB24####
LB24HALD= Phenos %>% 
  filter( timepoint %in% c(2, 4)) %>% 
  filter( cage.treat2 %in% c("LB")) %>% 
  dplyr::summarise(Treatment= "+LB",
                   Interval = "24",
                   LDMd = (logLDM[which(timepoint==4)]) - (logLDM[which(timepoint==2)]),
                   LDFd = (logLDF[which(timepoint==4)]) - (logLDF[which(timepoint==2)]),
                   SRMd = (logSRM[which(timepoint==4)]) - (logSRM[which(timepoint==2)]),
                   SRFd = (logSRF[which(timepoint==4)]) - (logSRF[which(timepoint==2)]),
                   FECd = (logFec[which(timepoint==4)]) - (logFec[which(timepoint==2)]),
                   VIABd = (logViab[which(timepoint==4)]) - (logViab[which(timepoint==2)]),
                   DWd = (logDW[which(timepoint==4)]) - (logDW[which(timepoint==2)]),
                   LWd = (logLW[which(timepoint==4)]) - (logLW[which(timepoint==2)]),
                   LDMSp = sqrt((((logLDM[which(timepoint==2)]-1)*(logLDMsd[which(timepoint==2)]^2)) + ((logLDM[which(timepoint==4)]-1)*(logLDMsd[which(timepoint==4)]^2)))/
                                  ((logLDM[which(timepoint==2)] + logLDM[which(timepoint==4)])-2)),
                   LDFSp = sqrt((((logLDF[which(timepoint==2)]-1)*(logLDFsd[which(timepoint==2)]^2)) + ((logLDF[which(timepoint==4)]-1)*(logLDFsd[which(timepoint==4)]^2)))/
                                  ((logLDF[which(timepoint==2)] + logLDF[which(timepoint==4)])-2)),
                   SRMSp = sqrt((((logSRM[which(timepoint==2)]-1)*(logSRMsd[which(timepoint==2)]^2)) + ((logSRM[which(timepoint==4)]-1)*(logSRMsd[which(timepoint==4)]^2)))/
                                  ((logSRM[which(timepoint==2)] + logSRM[which(timepoint==4)])-2)),
                   SRFSp = sqrt((((logSRF[which(timepoint==2)]-1)*(logSRFsd[which(timepoint==2)]^2)) + ((logSRF[which(timepoint==4)]-1)*(logSRFsd[which(timepoint==4)]^2)))/
                                  ((logSRF[which(timepoint==2)] + logSRF[which(timepoint==4)])-2)),
                   FECSp = sqrt((((logFec[which(timepoint==2)]-1)*(logFecsd[which(timepoint==2)]^2)) + ((logFec[which(timepoint==4)]-1)*(logFecsd[which(timepoint==4)]^2)))/
                                  ((logFec[which(timepoint==2)] + logFec[which(timepoint==4)])-2)),
                   VIABSp = sqrt((((logViab[which(timepoint==2)]-1)*(logViabsd[which(timepoint==2)]^2)) + ((logViab[which(timepoint==4)]-1)*(logViabsd[which(timepoint==4)]^2)))/
                                   ((logViab[which(timepoint==2)] + logViab[which(timepoint==4)])-2)),
                   DWSp = sqrt((((logDW[which(timepoint==2)]-1)*(logDWsd[which(timepoint==2)]^2)) + ((logDW[which(timepoint==4)]-1)*(logDWsd[which(timepoint==4)]^2)))/
                                 ((logDW[which(timepoint==2)] + logDW[which(timepoint==4)])-2)),
                   LWSp = sqrt((((logLW[which(timepoint==2)]-1)*(logLWsd[which(timepoint==2)]^2)) + ((logLW[which(timepoint==4)]-1)*(logLWsd[which(timepoint==4)]^2)))/
                                 ((logLW[which(timepoint==2)] + logLW[which(timepoint==4)])-2)),
                   I=2.640433,
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I),
                   DW_HAL = ((DWd/DWSp)/I),
                   LW_HAL = ((LWd/LWSp)/I))
LB24HALD <- LB24HALD[ -c(3:19) ]
View(LB24HALD)                   


####C25####
C25HALD = Phenos %>% 
  filter( timepoint %in% c(2, 5)) %>% 
  filter( cage.treat2 %in% c("C")) %>% 
  dplyr::summarise(Treatment= "Control",
                   Interval = "25",
                   LDMd = (logLDM[which(timepoint==5)]) - (logLDM[which(timepoint==2)]),
                   LDFd = (logLDF[which(timepoint==5)]) - (logLDF[which(timepoint==2)]),
                   SRMd = (logSRM[which(timepoint==5)]) - (logSRM[which(timepoint==2)]),
                   SRFd = (logSRF[which(timepoint==5)]) - (logSRF[which(timepoint==2)]),
                   FECd = (logFec[which(timepoint==5)]) - (logFec[which(timepoint==2)]),
                   VIABd = (logViab[which(timepoint==5)]) - (logViab[which(timepoint==2)]),
                   LDMSp = sqrt((((logLDM[which(timepoint==2)]-1)*(logLDMsd[which(timepoint==2)]^2)) + ((logLDM[which(timepoint==5)]-1)*(logLDMsd[which(timepoint==5)]^2)))/
                                  ((logLDM[which(timepoint==2)] + logLDM[which(timepoint==5)])-2)),
                   LDFSp = sqrt((((logLDF[which(timepoint==2)]-1)*(logLDFsd[which(timepoint==2)]^2)) + ((logLDF[which(timepoint==5)]-1)*(logLDFsd[which(timepoint==5)]^2)))/
                                  ((logLDF[which(timepoint==2)] + logLDF[which(timepoint==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(timepoint==2)]-1)*(logSRMsd[which(timepoint==2)]^2)) + ((logSRM[which(timepoint==5)]-1)*(logSRMsd[which(timepoint==5)]^2)))/
                                  ((logSRM[which(timepoint==2)] + logSRM[which(timepoint==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(timepoint==2)]-1)*(logSRFsd[which(timepoint==2)]^2)) + ((logSRF[which(timepoint==5)]-1)*(logSRFsd[which(timepoint==5)]^2)))/
                                  ((logSRF[which(timepoint==2)] + logSRF[which(timepoint==5)])-2)),
                   FECSp = sqrt((((logFec[which(timepoint==2)]-1)*(logFecsd[which(timepoint==2)]^2)) + ((logFec[which(timepoint==5)]-1)*(logFecsd[which(timepoint==5)]^2)))/
                                  ((logFec[which(timepoint==2)] + logFec[which(timepoint==5)])-2)),
                   VIABSp = sqrt((((logViab[which(timepoint==2)]-1)*(logViabsd[which(timepoint==2)]^2)) + ((logViab[which(timepoint==5)]-1)*(logViabsd[which(timepoint==5)]^2)))/
                                   ((logViab[which(timepoint==2)] + logViab[which(timepoint==5)])-2)),
                   I=2.700064,
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I))

C25HALD<- C25HALD[ -c(3:15) ] 
View(C25HALD) 

#View(A23HALD)


####AT25####
AT25HALD = Phenos %>% 
  filter( timepoint %in% c(2, 5)) %>% 
  filter( cage.treat2 %in% c("AT")) %>% 
  dplyr::summarise(Treatment= "+AT",
                   Interval = "25",
                   LDMd = (logLDM[which(timepoint==5)]) - (logLDM[which(timepoint==2)]),
                   LDFd = (logLDF[which(timepoint==5)]) - (logLDF[which(timepoint==2)]),
                   SRMd = (logSRM[which(timepoint==5)]) - (logSRM[which(timepoint==2)]),
                   SRFd = (logSRF[which(timepoint==5)]) - (logSRF[which(timepoint==2)]),
                   FECd = (logFec[which(timepoint==5)]) - (logFec[which(timepoint==2)]),
                   VIABd = (logViab[which(timepoint==5)]) - (logViab[which(timepoint==2)]),
                   LDMSp = sqrt((((logLDM[which(timepoint==2)]-1)*(logLDMsd[which(timepoint==2)]^2)) + ((logLDM[which(timepoint==5)]-1)*(logLDMsd[which(timepoint==5)]^2)))/
                                  ((logLDM[which(timepoint==2)] + logLDM[which(timepoint==5)])-2)),
                   LDFSp = sqrt((((logLDF[which(timepoint==2)]-1)*(logLDFsd[which(timepoint==2)]^2)) + ((logLDF[which(timepoint==5)]-1)*(logLDFsd[which(timepoint==5)]^2)))/
                                  ((logLDF[which(timepoint==2)] + logLDF[which(timepoint==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(timepoint==2)]-1)*(logSRMsd[which(timepoint==2)]^2)) + ((logSRM[which(timepoint==5)]-1)*(logSRMsd[which(timepoint==5)]^2)))/
                                  ((logSRM[which(timepoint==2)] + logSRM[which(timepoint==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(timepoint==2)]-1)*(logSRFsd[which(timepoint==2)]^2)) + ((logSRF[which(timepoint==5)]-1)*(logSRFsd[which(timepoint==5)]^2)))/
                                  ((logSRF[which(timepoint==2)] + logSRF[which(timepoint==5)])-2)),
                   FECSp = sqrt((((logFec[which(timepoint==2)]-1)*(logFecsd[which(timepoint==2)]^2)) + ((logFec[which(timepoint==5)]-1)*(logFecsd[which(timepoint==5)]^2)))/
                                  ((logFec[which(timepoint==2)] + logFec[which(timepoint==5)])-2)),
                   VIABSp = sqrt((((logViab[which(timepoint==2)]-1)*(logViabsd[which(timepoint==2)]^2)) + ((logViab[which(timepoint==5)]-1)*(logViabsd[which(timepoint==5)]^2)))/
                                   ((logViab[which(timepoint==2)] + logViab[which(timepoint==5)])-2)),
                   I=2.700064,
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I))

AT25HALD<- AT25HALD[ -c(3:15) ] 
View(AT25HALD) 



####LB25####
LB25HALD = Phenos %>% 
  filter( timepoint %in% c(2, 5)) %>% 
  filter( cage.treat2 %in% c("LB")) %>% 
  dplyr::summarise(Treatment= "+LB",
                   Interval = "25",
                   LDMd = (logLDM[which(timepoint==5)]) - (logLDM[which(timepoint==2)]),
                   LDFd = (logLDF[which(timepoint==5)]) - (logLDF[which(timepoint==2)]),
                   SRMd = (logSRM[which(timepoint==5)]) - (logSRM[which(timepoint==2)]),
                   SRFd = (logSRF[which(timepoint==5)]) - (logSRF[which(timepoint==2)]),
                   FECd = (logFec[which(timepoint==5)]) - (logFec[which(timepoint==2)]),
                   VIABd = (logViab[which(timepoint==5)]) - (logViab[which(timepoint==2)]),
                   LDMSp = sqrt((((logLDM[which(timepoint==2)]-1)*(logLDMsd[which(timepoint==2)]^2)) + ((logLDM[which(timepoint==5)]-1)*(logLDMsd[which(timepoint==5)]^2)))/
                                  ((logLDM[which(timepoint==2)] + logLDM[which(timepoint==5)])-2)),
                   LDFSp = sqrt((((logLDF[which(timepoint==2)]-1)*(logLDFsd[which(timepoint==2)]^2)) + ((logLDF[which(timepoint==5)]-1)*(logLDFsd[which(timepoint==5)]^2)))/
                                  ((logLDF[which(timepoint==2)] + logLDF[which(timepoint==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(timepoint==2)]-1)*(logSRMsd[which(timepoint==2)]^2)) + ((logSRM[which(timepoint==5)]-1)*(logSRMsd[which(timepoint==5)]^2)))/
                                  ((logSRM[which(timepoint==2)] + logSRM[which(timepoint==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(timepoint==2)]-1)*(logSRFsd[which(timepoint==2)]^2)) + ((logSRF[which(timepoint==5)]-1)*(logSRFsd[which(timepoint==5)]^2)))/
                                  ((logSRF[which(timepoint==2)] + logSRF[which(timepoint==5)])-2)),
                   FECSp = sqrt((((logFec[which(timepoint==2)]-1)*(logFecsd[which(timepoint==2)]^2)) + ((logFec[which(timepoint==5)]-1)*(logFecsd[which(timepoint==5)]^2)))/
                                  ((logFec[which(timepoint==2)] + logFec[which(timepoint==5)])-2)),
                   VIABSp = sqrt((((logViab[which(timepoint==2)]-1)*(logViabsd[which(timepoint==2)]^2)) + ((logViab[which(timepoint==5)]-1)*(logViabsd[which(timepoint==5)]^2)))/
                                   ((logViab[which(timepoint==2)] + logViab[which(timepoint==5)])-2)),
                   I=2.700064,
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I))

LB25HALD<- LB25HALD[ -c(3:15) ] 
View(LB25HALD) 



####C45####
C45HALD = Phenos %>% 
  filter( timepoint %in% c(4, 5)) %>% 
  filter( cage.treat2 %in% c("C")) %>% 
  dplyr::summarise(Treatment= "Control",
                   Interval = "45",
                   LDMd = (logLDM[which(timepoint==5)]) - (logLDM[which(timepoint==4)]),
                   LDFd = (logLDF[which(timepoint==5)]) - (logLDF[which(timepoint==4)]),
                   SRMd = (logSRM[which(timepoint==5)]) - (logSRM[which(timepoint==4)]),
                   SRFd = (logSRF[which(timepoint==5)]) - (logSRF[which(timepoint==4)]),
                   FECd = (logFec[which(timepoint==5)]) - (logFec[which(timepoint==4)]),
                   VIABd = (logViab[which(timepoint==5)]) - (logViab[which(timepoint==4)]),
                   LDMSp = sqrt((((logLDM[which(timepoint==4)]-1)*(logLDMsd[which(timepoint==4)]^2)) + ((logLDM[which(timepoint==5)]-1)*(logLDMsd[which(timepoint==5)]^2)))/
                                  ((logLDM[which(timepoint==4)] + logLDM[which(timepoint==5)])-2)),
                   LDFSp = sqrt((((logLDF[which(timepoint==4)]-1)*(logLDFsd[which(timepoint==4)]^2)) + ((logLDF[which(timepoint==5)]-1)*(logLDFsd[which(timepoint==5)]^2)))/
                                  ((logLDF[which(timepoint==4)] + logLDF[which(timepoint==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(timepoint==4)]-1)*(logSRMsd[which(timepoint==4)]^2)) + ((logSRM[which(timepoint==5)]-1)*(logSRMsd[which(timepoint==5)]^2)))/
                                  ((logSRM[which(timepoint==4)] + logSRM[which(timepoint==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(timepoint==4)]-1)*(logSRFsd[which(timepoint==4)]^2)) + ((logSRF[which(timepoint==5)]-1)*(logSRFsd[which(timepoint==5)]^2)))/
                                  ((logSRF[which(timepoint==4)] + logSRF[which(timepoint==5)])-2)),
                   FECSp = sqrt((((logFec[which(timepoint==4)]-1)*(logFecsd[which(timepoint==4)]^2)) + ((logFec[which(timepoint==5)]-1)*(logFecsd[which(timepoint==5)]^2)))/
                                  ((logFec[which(timepoint==4)] + logFec[which(timepoint==5)])-2)),
                   VIABSp = sqrt((((logViab[which(timepoint==4)]-1)*(logViabsd[which(timepoint==4)]^2)) + ((logViab[which(timepoint==5)]-1)*(logViabsd[which(timepoint==5)]^2)))/
                                   ((logViab[which(timepoint==4)] + logViab[which(timepoint==5)])-2)),
                   I=0.1009917,
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I))

C45HALD<- C45HALD[ -c(3:15) ] 
View(C45HALD) 

#View(A23HALD)



####AT45####
AT45HALD = Phenos %>% 
  filter( timepoint %in% c(4, 5)) %>% 
  filter( cage.treat2 %in% c("AT")) %>% 
  dplyr::summarise(Treatment= "+AT",
                   Interval = "45",
                   LDMd = (logLDM[which(timepoint==5)]) - (logLDM[which(timepoint==4)]),
                   LDFd = (logLDF[which(timepoint==5)]) - (logLDF[which(timepoint==4)]),
                   SRMd = (logSRM[which(timepoint==5)]) - (logSRM[which(timepoint==4)]),
                   SRFd = (logSRF[which(timepoint==5)]) - (logSRF[which(timepoint==4)]),
                   FECd = (logFec[which(timepoint==5)]) - (logFec[which(timepoint==4)]),
                   VIABd = (logViab[which(timepoint==5)]) - (logViab[which(timepoint==4)]),
                   LDMSp = sqrt((((logLDM[which(timepoint==4)]-1)*(logLDMsd[which(timepoint==4)]^2)) + ((logLDM[which(timepoint==5)]-1)*(logLDMsd[which(timepoint==5)]^2)))/
                                  ((logLDM[which(timepoint==4)] + logLDM[which(timepoint==5)])-2)),
                   LDFSp = sqrt((((logLDF[which(timepoint==4)]-1)*(logLDFsd[which(timepoint==4)]^2)) + ((logLDF[which(timepoint==5)]-1)*(logLDFsd[which(timepoint==5)]^2)))/
                                  ((logLDF[which(timepoint==4)] + logLDF[which(timepoint==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(timepoint==4)]-1)*(logSRMsd[which(timepoint==4)]^2)) + ((logSRM[which(timepoint==5)]-1)*(logSRMsd[which(timepoint==5)]^2)))/
                                  ((logSRM[which(timepoint==4)] + logSRM[which(timepoint==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(timepoint==4)]-1)*(logSRFsd[which(timepoint==4)]^2)) + ((logSRF[which(timepoint==5)]-1)*(logSRFsd[which(timepoint==5)]^2)))/
                                  ((logSRF[which(timepoint==4)] + logSRF[which(timepoint==5)])-2)),
                   FECSp = sqrt((((logFec[which(timepoint==4)]-1)*(logFecsd[which(timepoint==4)]^2)) + ((logFec[which(timepoint==5)]-1)*(logFecsd[which(timepoint==5)]^2)))/
                                  ((logFec[which(timepoint==4)] + logFec[which(timepoint==5)])-2)),
                   VIABSp = sqrt((((logViab[which(timepoint==4)]-1)*(logViabsd[which(timepoint==4)]^2)) + ((logViab[which(timepoint==5)]-1)*(logViabsd[which(timepoint==5)]^2)))/
                                   ((logViab[which(timepoint==4)] + logViab[which(timepoint==5)])-2)),
                   I=0.1009917,
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I))

AT45HALD<- AT45HALD[ -c(3:15) ] 
View(AT45HALD) 

#View(A23HALD)



####LB45####
LB45HALD = Phenos %>% 
  filter( timepoint %in% c(4, 5)) %>% 
  filter( cage.treat2 %in% c("LB")) %>% 
  dplyr::summarise(Treatment= "+LB",
                   Interval = "45",
                   LDMd = (logLDM[which(timepoint==5)]) - (logLDM[which(timepoint==4)]),
                   LDFd = (logLDF[which(timepoint==5)]) - (logLDF[which(timepoint==4)]),
                   SRMd = (logSRM[which(timepoint==5)]) - (logSRM[which(timepoint==4)]),
                   SRFd = (logSRF[which(timepoint==5)]) - (logSRF[which(timepoint==4)]),
                   FECd = (logFec[which(timepoint==5)]) - (logFec[which(timepoint==4)]),
                   VIABd = (logViab[which(timepoint==5)]) - (logViab[which(timepoint==4)]),
                   LDMSp = sqrt((((logLDM[which(timepoint==4)]-1)*(logLDMsd[which(timepoint==4)]^2)) + ((logLDM[which(timepoint==5)]-1)*(logLDMsd[which(timepoint==5)]^2)))/
                                  ((logLDM[which(timepoint==4)] + logLDM[which(timepoint==5)])-2)),
                   LDFSp = sqrt((((logLDF[which(timepoint==4)]-1)*(logLDFsd[which(timepoint==4)]^2)) + ((logLDF[which(timepoint==5)]-1)*(logLDFsd[which(timepoint==5)]^2)))/
                                  ((logLDF[which(timepoint==4)] + logLDF[which(timepoint==5)])-2)),
                   SRMSp = sqrt((((logSRM[which(timepoint==4)]-1)*(logSRMsd[which(timepoint==4)]^2)) + ((logSRM[which(timepoint==5)]-1)*(logSRMsd[which(timepoint==5)]^2)))/
                                  ((logSRM[which(timepoint==4)] + logSRM[which(timepoint==5)])-2)),
                   SRFSp = sqrt((((logSRF[which(timepoint==4)]-1)*(logSRFsd[which(timepoint==4)]^2)) + ((logSRF[which(timepoint==5)]-1)*(logSRFsd[which(timepoint==5)]^2)))/
                                  ((logSRF[which(timepoint==4)] + logSRF[which(timepoint==5)])-2)),
                   FECSp = sqrt((((logFec[which(timepoint==4)]-1)*(logFecsd[which(timepoint==4)]^2)) + ((logFec[which(timepoint==5)]-1)*(logFecsd[which(timepoint==5)]^2)))/
                                  ((logFec[which(timepoint==4)] + logFec[which(timepoint==5)])-2)),
                   VIABSp = sqrt((((logViab[which(timepoint==4)]-1)*(logViabsd[which(timepoint==4)]^2)) + ((logViab[which(timepoint==5)]-1)*(logViabsd[which(timepoint==5)]^2)))/
                                   ((logViab[which(timepoint==4)] + logViab[which(timepoint==5)])-2)),
                   I=0.1009917,
                   LDM_HAL = ((LDMd/LDMSp)/I),
                   LDF_HAL = ((LDFd/LDFSp)/I),
                   SRM_HAL = ((SRMd/SRMSp)/I),
                   SRF_HAL = ((SRFd/SRFSp)/I),
                   FEC_HAL = ((FECd/FECSp)/I),
                   VIAB_HAL = ((VIABd/VIABSp)/I))

LB45HALD<- LB45HALD[ -c(3:15) ] 
View(LB45HALD) 


####Merge all dfs####
HALD_24 <- rbind(C24HALD,AT24HALD,LB24HALD)
HALD_25 <- rbind(C25HALD,AT25HALD,LB25HALD)
HALD_45 <- rbind(C45HALD,AT45HALD,LB45HALD)

View(HALD_24_long)

### make them long###
HALD_24_long <- melt(HALD_24 ,  id.vars = c('Treatment','Interval'), variable.name = 'Phenotypes')
HALD_24_long$absValue <- abs(HALD_24_long$value)
HALD_24_long$logValue <- log(HALD_24_long$absValue)

HALD_25_long <- melt(HALD_25 ,  id.vars = c('Treatment','Interval'), variable.name = 'Phenotypes')
HALD_25_long$absValue <- abs(HALD_25_long$value)
HALD_25_long$logValue <- log(HALD_25_long$absValue)

HALD_45_long <- melt(HALD_45 ,  id.vars = c('Treatment','Interval'), variable.name = 'Phenotypes')
HALD_45_long$absValue <- abs(HALD_45_long$value)
HALD_45_long$logValue <- log(HALD_45_long$absValue)

HALD_LONG_combined=bind_rows(HALD_24_long,HALD_25_long,HALD_45_long)
View(HALD_LONG_combined)
#create line plot for each column in data frame
ggplot(HALD_24_long, aes(Interval, logValue)) +
  geom_point(aes(colour = Phenotypes), size = 3) + 
  scale_colour_manual(breaks=c("LDM_HAL","LDF_HAL","SRM_HAL","SRF_HAL","FEC_HAL","VIAB_HAL","DW_HAL","LW_HAL","PIG_HAL"), 
                      labels =c("Dev Time (male)", "Dev Time (female)", "Starv Time (male)", "Starv Time (female)", "Fecundity", "Viability", "Dry Weight","Lipid Weight","Pigmentation"), 
                      values =c("red","coral4", "darkgreen", "green", "gold","blue","grey","black","coral"))+
  scale_x_discrete(breaks = c("01", "12", "23", "34", "45", "15", "05"), labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "1-5", "0-5") )+
  ylab("log(Haldanes)") +
  xlab("Timepoint Interval") +
  facet_wrap(~ Treatment) + 
  theme_cowplot()

ggplot(HALD_25_long, aes(Interval, logValue)) +
  geom_point(aes(colour = Phenotypes), size = 3) + 
  scale_colour_manual(breaks=c("LDM_HAL","LDF_HAL","SRM_HAL","SRF_HAL","FEC_HAL","VIAB_HAL","PIG_HAL"), 
                      labels =c("Dev Time (male)", "Dev Time (female)", "Starv Time (male)", "Starv Time (female)", "Fecundity", "Viability","Pigmentation"), 
                      values =c("red","coral4", "darkgreen", "green", "gold","blue","coral"))+
  scale_x_discrete(breaks = c("01", "12", "23", "34", "45", "15", "05"), labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "1-5", "0-5") )+
  ylab("log(Haldanes)") +
  xlab("Timepoint Interval") +
  facet_wrap(~ Treatment) + 
  theme_cowplot()

ggplot(HALD_45_long, aes(Interval, logValue)) +
  geom_point(aes(colour = Phenotypes), size = 3) + 
  scale_colour_manual(breaks=c("LDM_HAL","LDF_HAL","SRM_HAL","SRF_HAL","FEC_HAL","VIAB_HAL","PIG_HAL"), 
                      labels =c("Dev Time (male)", "Dev Time (female)", "Starv Time (male)", "Starv Time (female)", "Fecundity", "Viability","Pigmentation"), 
                      values =c("red","coral4", "darkgreen", "green", "gold","blue","coral"))+
  scale_x_discrete(breaks = c("01", "12", "23", "34", "45", "15", "05"), labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "1-5", "0-5") )+
  ylab("log(Haldanes)") +
  xlab("Timepoint Interval") +
  facet_wrap(~ Treatment) + 
  theme_cowplot()

ggplot(HALD_LONG_combined, aes(Interval, logValue)) +
  geom_point(aes(colour = Phenotypes), size = 3) + 
  scale_colour_manual(breaks=c("LDM_HAL","LDF_HAL","SRM_HAL","SRF_HAL","FEC_HAL","VIAB_HAL","DW_HAL","LW_HAL","PIG_HAL"), 
                      labels =c("Dev Time (male)", "Dev Time (female)", "Starv Time (male)", "Starv Time (female)", "Fecundity", "Viability", "Dry Weight","Lipid Weight","Pigmentation"), 
                      values =c("red","coral4", "darkgreen", "green", "gold","blue","grey","black","coral"))+
  scale_x_discrete(breaks = c("01", "12", "23", "34", "45", "15", "05"), labels = c("0-1", "1-2", "2-3", "3-4", "4-5", "1-5", "0-5") )+
  ylab("log(Haldanes)") +
  xlab("Timepoint Interval") +
  facet_wrap(~ Treatment) + 
  theme_cowplot()


View(HALD_LONG_combined)
HALD_LONG_combined$absValue=as.numeric(HALD_LONG_combined$logValue)

Hald_combinedlog=ggplot(HALD_LONG_combined, aes(Treatment,logValue, color=Interval)) +
  geom_jitter(size=5, width = .3) + 
  scale_color_manual(breaks = c( "24", "25", "45"),
                     labels = c( "1-2", "1-3","2-3"),
                     values = c("#D8D8D8","#909090","red"))+
  scale_y_continuous(name="log(Haldanes)") +
  theme_cowplot(10) 
Hald_combinedlog

Hald_combinedval=ggplot(HALD_LONG_combined, aes(Treatment,value, color=Interval)) +
  geom_jitter(size=5, width = .3) + 
  scale_color_manual(breaks = c( "24", "25", "45"),
                     labels = c( "2-4", "2-5","4-5"),
                     values = c("#D8D8D8","#909090","red"))+
  scale_y_continuous(name="Rate of Phenotypic Change (Haldanes)") +
  theme_cowplot(10) 
Hald_combinedval


load("fig2.rdata") 
Fig2

MAFig2H <- ggarrange(Hald_combinedval, Hald_combinedlog,
                   labels = c("A", "B"),
                   widths = c(1,1),
                   ncol = 2, nrow = 1, common.legend = TRUE)
MAFig2H
save(Hald_combinedlog, file = "Hald_combinedlog.rdata")


###LIinear analysis of haldanes ####

HALDLMM <- lm(logValue ~ Interval*Treatment, data = HALD_LONG_combined)  
anova(HALDLMM)
lsmeans(HALDLMM, pairwise ~ Treatment)
lsmeans(HALDLMM, pairwise ~ Interval)
