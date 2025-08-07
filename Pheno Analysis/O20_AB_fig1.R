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
library(lsmeans)

####Import Files####
weather = read_csv("raw/O20_weather_raw.csv")
census = read_csv("raw/O20AB_CS_raw.csv")
#View(weather)

##separating date and time##

weather <- separate(weather, TimeStamp, into = c("Date", "Time"), sep = "^\\S*\\K\\s+")
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
                   RH=mean(`RH (percent)`,na.rm = TRUE),
                   RHsd=sd(`RH (percent)`,na.rm = TRUE),
                   RHse=RHsd/sqrt(N),
                   Temp_C=mean(Temp_C,na.rm = TRUE),
                   Temp_Csd=sd(Temp_C,na.rm = TRUE),
                   Temp_Cse= Temp_Csd/sqrt(N),
                   Temp_Cavg=mean(Temp_C_Avg,na.rm = TRUE),
                   Temp_Cavgsd=sd(Temp_C_Avg,na.rm = TRUE),
                   Temp_Cavgse= Temp_Cavgsd/sqrt(N),
                   PARden=mean(PAR_Den_Avg,na.rm = TRUE),
                   PARdensd=sd(PAR_Den_Avg,na.rm = TRUE),
                   PARdense= PARdensd/sqrt(N),
                   PARtot=mean(PAR_Tot_Tot,na.rm = TRUE),
                   PARtotsd=sd(PAR_Tot_Tot,na.rm = TRUE),
                   PARtotse= PARtotsd/sqrt(N),
                   WS=mean(WS_ms_Avg,na.rm =TRUE),
                   WSsd=sd(WS_ms_Avg,na.rm = TRUE),
                   WSse= WSsd/sqrt(N),
                   Rain=mean(Rain_mm_Tot,na.rm = TRUE),
                   Rainsd=sd(Rain_mm_Tot,na.rm = TRUE),
                   Rainse=Rainsd/sqrt(N)
  )
wdata=unique(weather.sum)
#View(wdata)

###graph min/max temperature####

AirTCMin <-ggplot(wdata, aes(x=as.factor(Date), y=MaxAirTC, group=1))+
  geom_line(color="Red")+
  #geom_linerange(aes(ymin=AirTC-AirTCse,ymax=AirTC+AirTCse))+
  ylab("Max/Min Daily Temperature C")+
  scale_x_discrete(name=" Experiment Date", 
                   breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
                   labels=c("Start","August", "September", "October", "November", "End")) +
  #xlab("Sample Timepoints")+
  #scale_x_discrete(labels=c("Founder", "Summer", "Fall"))+
  #scale_color_discrete(name = "Cage Treatment", labels = c("Apple Food", "Apple + Ac. thailandicus", "Apple + Lac. brevis"))+
  theme_cowplot()+
  theme(legend.position = "none")

AirTCMinMax<-AirTCMin+geom_line(linetype= "solid" , data=wdata, aes(x=as.factor(Date), y=MinAIRTC, group=1), color="Blue")
AirTCMinMax


####census data####
#View(census)
census= census %>% 
  rename(
    pic.avg = `pic avg`)

census$Date=as.Date(census$Date, format = "%m/%d/%y")

Census.sum=census %>%
  group_by(Date, Cage_treatment) %>%
  dplyr::summarise(N=n(),
                   countavg=mean(Picture_total,na.rm = TRUE),
                   countsd=sd(Picture_total,na.rm = TRUE),
                   countse= countsd/sqrt(N),
                   combestavg=mean(est_total,na.rm = TRUE),
                   combestsd=sd(est_total,na.rm = TRUE),
                   combestse=combestsd/sqrt(N), 
                   picavgest=mean(est_avg,na.rm = TRUE),
                   picavgsd=sd(est_avg,na.rm = TRUE),
                   picavgse=picavgsd/sqrt(N),
                   pictotavg=mean(est_tot,na.rm = TRUE),
                   pictotsd=sd(est_tot,na.rm = TRUE),
                   pictotse= pictotsd/sqrt(N),
  )
#View(Census.sum)

cagecensus.sum=census %>%
  group_by(Date, Cage_treatment, Cage_number) %>%
  dplyr::summarise(N=n(),
                   countavg=mean(Picture_total,na.rm = TRUE),
                   countsd=sd(Picture_total,na.rm = TRUE),
                   countse= countsd/sqrt(N),
                   combestavg=mean(est_total,na.rm = TRUE),
                   combestsd=sd(est_total,na.rm = TRUE),
                   combestse=combestsd/sqrt(N), 
                   picavgest=mean(est_avg,na.rm = TRUE),
                   picavgsd=sd(est_avg,na.rm = TRUE),
                   picavgse=picavgsd/sqrt(N),
                   pictotavg=mean(est_tot,na.rm = TRUE),
                   pictotsd=sd(est_tot,na.rm = TRUE),
                   pictotse= pictotsd/sqrt(N),
                   
                   
  )


cagecensusdata=unique(cagecensus.sum)

##combine dfs####
Census.sum$Date=as.Date(Census.sum$Date)
cagecensusdata$Date=as.Date(cagecensusdata$Date)
wdata$Date=as.Date(wdata$Date)

dplrcombinedtreatment <- left_join(wdata, cagecensusdata, by="Date")
dplrtreatcombinedtreatment <- left_join(wdata, Census.sum, by="Date")
#view(dplrcombinedtreatment)
###make plot###
coeff=8000

test3=ggplot(dplrcombinedtreatment, aes(x=as.factor(Date))) +
  geom_line(aes(y=combestavg, group=dplrcombinedtreatment$Cage_number , color=dplrcombinedtreatment$Cage_treatment), size=1, alpha=.4)+
  #geom_errorbar(aes(ymin=combestavg-combestse,ymax=combestavg+combestse, color=dplrcombinedtreatment$Cage_treatment), size=.5)+
  scale_colour_manual(name = " ", 
                      breaks = c("a", "b"),
                      labels = c("Low Quality", "High Quality"),
                      values = c("grey","black"))+
  #geom_linerange(aes(ymin=combestavg-combestse,ymax=combestavg+combestse), size=.5)+
  geom_line( aes(y=MinAIRTC * coeff), linetype="dotted", size=1, group=1, color="blue2", alpha=.3) + 
  geom_line( aes(y=MaxAirTC * coeff), linetype="dashed",size=1, group=1, color="red1", alpha=.3) +
  scale_x_discrete(
    name=" Experiment Date", 
    breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
    labels=c("Start","August", "September", "October", "November", "End"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
    # Features of the first axis
    name = "Population Size",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Min/ Max Daily Temperature (C°)"))+
  theme_cowplot()

test3+ theme(legend.position = c(.9, .99),
             legend.justification = c("right", "top"),
             legend.box.just = "right",
             legend.margin = margin(.5, .5, .5,.5), 
             legend.text=element_text(size=14),
             legend.title = element_blank())
########
dplrcombinedtreatment2=dplrcombinedtreatment %>% mutate(LOGcount = log10(combestavg) ) 
ggqqplot(dplrcombinedtreatment2$LOGcount)
coeff=6000

test4=ggplot(dplrtreatcombinedtreatment, aes(x=as.factor(Date))) +
  geom_line(aes(y=combestavg, group=dplrtreatcombinedtreatment$Cage_treatment , color=dplrtreatcombinedtreatment$Cage_treatment), size=2)+
  #geom_errorbar(aes(ymin=combestavg-combestse,ymax=combestavg+combestse, color=dplrtreatcombinedtreatment$Cage_treatment), size=1)+
  scale_colour_manual(name = " ", 
                      breaks = c("a", "b"),
                      labels = c("Low Quality", "High Quality"),
                      values = c("grey","black"))+
  
  #geom_line(aes(y=combestavg, group=dplrcombinedtreatment$Cage_number , color=dplrcombinedtreatment$Cage_treatment), size=1)+
  #geom_linerange(aes(ymin=combestavg-combestse,ymax=combestavg+combestse), size=.5)+
  geom_line( aes(y=MinAIRTC * coeff), linetype="dotted", size=1, group=1, color="blue2", alpha=.6) + 
  geom_line( aes(y=MaxAirTC * coeff), linetype="dashed",size=1, group=1, color="red1", alpha=.6) +
  scale_x_discrete(
    name=" Experiment Date", 
    breaks=c("2020-07-21","2020-08-14", "2020-09-18", "2020-10-18", "2020-11-10", "2020-11-21"),
    labels=c("Start","August", "September", "October", "November", "End"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     # Features of the first axis
                     name = "Population Size",
                     # Add a second axis and specify its features
                     sec.axis = sec_axis(~./coeff, name="Min/ Max Daily Temperature (C°)"))+
  theme_cowplot()

test4+ theme(legend.position = c(.9, .99),
             legend.justification = c("right", "top"),
             legend.box.just = "right",
             legend.margin = margin(.5, .5, .5,.5), 
             legend.text=element_text(size=14),
             legend.title = element_blank())

####linear regression analysis####
#
##All Timepoints##
CensusALL = census
names(census)[9] ="picavg"
CensusALL = census
CensusALL$Date=as.Date(CensusALL$Date)

CensusALL=CensusALL %>% 
  mutate(TP = rank(Date)) 

CensusALL$TP=as.numeric(CensusALL$TP)
CensusALL$Date=as.factor(CensusALL$Date)
CensusALL$DateCage <- paste(CensusALL$Date,CensusALL$Cage_number)
View(CensusALL)

CensusALLanv<-lme(picavg ~ Date*Cage_treatment, random=~1|Cage_number, data=CensusALL)
anova(CensusALLanv) ## sig Date, insig cage_treatment and interaction 
summary(CensusALLanv)


pairs=lsmeans(CensusALLanv, pairwise ~ Cage_treatment|Date)
summ=summary(pairs$contrasts)
View(summ)


write.csv(summ, "lsmeans-pairwise-fig1.csv")
### Second Mixed-Effect model
install.packages("lmerTest")
library(lmerTest)

int_model1 <- lmer(picavg ~ Date*Cage_treatment + (1|Cage_number), data=CensusALL)
summary(int_model1)


.##Expanding timepoints##
CensusExpand=subset(census, census$Date == c("2020-08-14") | census$Date == c("2020-08-31"))
View(CensusExpand)
CensusExpand$Date=as.numeric(CensusExpand$Date)
CensusExpandanv<-lme(picavg ~ Date*Cage_treatment, random=~1|Cage_number/Date, data=CensusExpand)
anova(CensusExpandanv) ## Sig Date, insig interaction and treatment 
summary(CensusExpandanv)

##Peak timepoints##
CensusPeak=subset(census, census$Date == c("2020-09-04") | census$Date == c("2020-09-11") | census$Date == c("2020-09-18") | census$Date == c("2020-09-25"))
CensusPeak$Date=as.numeric(CensusPeak$Date)
CensusPeakanv<-lme(picavg ~ Date*Cage_treatment, random=~1|Cage_number, data=CensusPeak)
anova(CensusPeakanv) ## Sig Date  insig interaction and treatment 
summary(CensusPeakanv)

##Collapse timepoints##
CensusCollapse=subset(census, census$Date == c("2020-10-09") | census$Date == c("2020-10-18") | census$Date == c("2020-11-01") | census$Date == c("2020-11-10") | census$Date == c("2020-11-21"))
CensusCollapse$Date=as.numeric(CensusCollapse$Date)
CensusCollapseanv<-lme(picavg ~ Date*Cage_treatment, random=~1|Cage_number/Date, data=CensusCollapse)
anova(CensusCollapseanv) ## Sig Date and Treatment, insig interaction
summary(CensusCollapseanv)









