library(reshape2)
library(plyr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(latex2exp)
library(ggrepel)
library(scales)
library(grid)
library(ggtext)
library(ggpubr)

# PREPARATION
## Setting Working directory

ABpca2 <- Phenos.cage ### this df is produced when generating figure 2, must complete that first. 
ABpca2 <- ABpca2[, -c(8,9,10,12,13,15,16,18,19,21,22,23,24,25,27,28,30,31,33)]## remove SE and SD
ABpca2 <- subset(ABpca2, ABpca2$TP %in% c(0,2,4))
#View(ABpca2)


#### LISTS AND COLOR ASSIGNMENTS####

  
  
  AB2ColorListTime <- c("0"="grey",
                     "2"="yellow",
                     "4"= "blue")
  
  
  
  AB2ColorListPopulation<-c("aa"="#f40006", "ab" = "#8B0000","ba" ="#707070",
                            "bb"="black")
  
  #View(AB2ColorListPopulation)
  
  AB2ColorListMDS=c("aa"="#f40006", "ab" = "#8B0000","ba" ="#707070",
                    "bb"="black")
  
  AB2ColorListTrait<-c(
    "LDf"="black",
    "LDf"="black",
    "Viab"="black",
    "SRM"="black",
    "SRF"="black",
    "TOTF"="black",
    "BS"="black")
  
  
  #display.brewer.pal(8,"Dark2")
  
  AB2TimeList=c("0","2","4")
  AB2PopulationList=c("aa","ab", "ba","bb")
  AB2TraitsList<-data.frame(ABBRV=c("LDm","LDf","Viab","SRM","SRF","TOTF","BS"),LONG=c("Development Time (Males)","Development Time (Females)", "Viability","Starvation Resistance (Males)","Starvation Resistance (Females)","Fecundity","Body Size"))
  AB2SigTrait3List<-data.frame(ABBRV=c("LDf","SRF","Viab","TOTF","BS"),LONG=c("Development Time (Females)","Starvation Resistance (Females)","Viability","Fecundity","Body Size"))
  
  AB2SelectionLabels=c("Founder","TP2","TP4" )
  AB2PopulationLabels=c("LQ","LQ (Phenotyed in HQ)","HQ (Phenotyped in LQ)","HQ" )
  AB2TimeShape=c("0"=8,"2"=17,"4"=19)
  


#SETTING A CONSISTENT IN ALL FIGURES THEME

{
  {
    BASESIZE=8
    TITLESIZE=7
    FONTFAMILY="sans"
    TEXTCOLOR="black"
    AXISTEXTSIZE=7
    AXISTEXTCOLOR="grey40"
    AXISTITLESIZE=7
    LEGENDTEXTSIZE=6
    LABELSIZE=11
  }
  
  theme_evo <- theme_set(theme_bw(base_size=BASESIZE,
                                  base_family=FONTFAMILY))
  theme_evo <- theme_update(panel.background = element_blank(),
                            title = element_text(size=TITLESIZE, face= "bold"),
                            axis.text = element_text(size= AXISTEXTSIZE, color = AXISTEXTCOLOR, family = FONTFAMILY, face= "bold"),
                            axis.title = element_text(size= TITLESIZE, color = TEXTCOLOR, family = FONTFAMILY, face= "bold"),
                            legend.background = element_blank(),
                            legend.box.background = element_blank(),
                            legend.title=element_text(face= "bold",size = LEGENDTEXTSIZE),
                            legend.text=element_text(size = LEGENDTEXTSIZE),
                            legend.key =element_blank()
  )
  
}

# EMPTY PLOT
pBL<-ggplot()+theme_nothing()


# PHENOTYPIC PCA - FIGURE 3
{
  
  #### PREP####
  {
    # FUNTION TO FIND HULLS
    find_hull<-function(df) df[chull(df$PC1,df$PC2),]
    
    # READING THE DATA
   # View(ABpca2)


    # SORT THE DATA
    ABpca2$Time<-factor(ABpca2$TP,AB2TimeList)
    ABpca2$Population=factor(ABpca2$Cage.Pheno,AB2PopulationList)
    ABpca2<-ABpca2 %>% arrange(Cage,Time,Population)
    
    #View(ABpca2)
    #### CALCULATING THE PRINCIPAL COMPONENTS####
    AB2_AA<- ABpca2 %>% subset(Population=="aa")
    AB2_BB<- ABpca2 %>% subset(Population=="bb")
    AB2_AB<- ABpca2 %>% subset(Population=="ab")
    AB2_BA<- ABpca2 %>% subset(Population=="ba")
    

    AB2_0<- ABpca2 %>% subset(Time=="0")
    AB2_2<- ABpca2 %>% subset(Time=="2")
    AB2_4<- ABpca2 %>% subset(Time=="4")
    
    
    #View(AB2_BA)
    ### all variables 

    AB2AA.pca<-prcomp(AB2_AA[8:14],scale=TRUE)
    AB2BB.pca<-prcomp(AB2_BB[8:14],scale=TRUE)
    AB2AB.pca<-prcomp(AB2_AB[8:14],scale=TRUE)
    AB2BA.pca<-prcomp(AB2_BA[8:14],scale=TRUE)
    
    AB2.pca<-prcomp(ABpca2[8:14],scale=TRUE)
    
    
    ### sig variables
    #View(ABpca2)
    
    AB2AA_S.pca<-prcomp(AB2_AA[c(9,11:14)],scale=TRUE)
    AB2BB_S.pca<-prcomp(AB2_BB[c(9,11:14)],scale=TRUE)
    AB2AB_S.pca<-prcomp(AB2_AB[c(9,11:14)],scale=TRUE)
    AB2BA_S.pca<-prcomp(AB2_BA[c(9,11:14)],scale=TRUE)
    
    AB2_S.pca<-prcomp(ABpca2[c(9,11:14)],scale=TRUE)
   
    # GETTING PCA COORDINATES TO A DATA FRAME
    ##All Vars
    AB2AA.pca.df<-data.frame(Cage=AB2_AA$Cage,Population=AB2_AA$Population,Time=AB2_AA$Time,AB2AA.pca$x)
    AB2BB.pca.df<-data.frame(Cage=AB2_BB$Cage,Population=AB2_BB$Population,Time=AB2_BB$Time,AB2BB.pca$x)
    AB2AB.pca.df<-data.frame(Cage=AB2_AB$Cage,Population=AB2_AB$Population,Time=AB2_AB$Time,AB2AB.pca$x)
    AB2BA.pca.df<-data.frame(Cage=AB2_BA$Cage,Population=AB2_BA$Population,Time=AB2_BA$Time,AB2BA.pca$x)
    
    AB2.pca.df<-data.frame(Cage=ABpca2$Cage,Population=ABpca2$Population,Treat=ABpca2$Cage.Treat.x,Pheno=ABpca2$Pheno.Treat,Time=ABpca2$Time,AB2.pca$x)
    
    View(AB2.pca.df)
    
    
    
    ####PCA PHENOTYPES COMBINED####
   Pheno24.pca1=ggplot(AB2.pca.df, aes(x=PC1, y=PC2, shape = Time, color = Pheno, group=Pheno)) +
      stat_ellipse()+
      geom_point(size=3)+
      scale_colour_manual(name = "Assay Diet", 
                          breaks = c("a", "b"),
                          labels = c("Low Quality","High Quality"),
                          values = c("grey","black"))+
      scale_shape_manual(breaks=c("0", "2", "4"), 
                         label=c("Founder","Timepoint 2 (Summer)", "Timepoint 4 (Fall)"), 
                         values =c(21,19,17))+
      theme_cowplot()+
      theme(legend.position = "bottom", legend.box = "horizontal")

    
    Pheno24.pca2=ggplot(AB2.pca.df, aes(x=PC1, y=PC2, shape = Time, color = Treat, group=Treat)) +
      stat_ellipse()+
      geom_point(size=3)+
      scale_colour_manual(name = "Treatment", 
                          breaks = c("a", "b"),
                          labels = c("Low Quality","High Quality"),
                          values = c("grey","black"))+
      scale_shape_manual(breaks=c("0", "2", "4"), 
                         label=c("Founder","Timepoint 2 (Summer)", "Timepoint 4 (Fall)"), 
                         values =c(21,19,17))+
      theme_cowplot()+
      theme(legend.position = "bottom", legend.box = "horizontal")
    
    Pheno24.PCA <- ggarrange(Pheno24.pca1, Pheno24.pca2, ncol=1, labels = c("A","B"), common.legend = FALSE, font.label = list(size = 20))
    
    Pheno24.PCA
    
    ab2.pca.df.manov<-manova(cbind(PC1, PC2)~Treat, AB2.pca.df)
    summary(ab2.pca.df.manov)
    
    ## manova for phenoed on apple / bloom all variables
    
    AB2.pca.AABA.df<-rbind(AB2AA.pca.df, AB2BA.pca.df)
    AB2.AABA.PCA.aov<-manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7)~Time*Population,AB2.pca.AABA.df)
    summary(AB2.AABA.PCA.aov) ### SIG TIME AND POPULATION NOT INTERCEPT
    
    AB2.pca.BBAB.df<-rbind(AB2BB.pca.df,AB2AB.pca.df)
    AB2.BBAB.PCA.aov<-manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7)~Time*Population, AB2.pca.BBAB.df)
    summary(AB2.BBAB.PCA.aov) ###SIG TIME NOT POPULATION OR INTERCEPT
    
    ##Sig Vars 
    AB2AA_S.pca.df<-data.frame(Cage=AB2_AA$Cage,Population=AB2_AA$Population,Time=AB2_AA$Time,AB2AA_S.pca$x)
    AB2BB_S.pca.df<-data.frame(Cage=AB2_BB$Cage,Population=AB2_BB$Population,Time=AB2_BB$Time,AB2BB_S.pca$x)
    AB2AB_S.pca.df<-data.frame(Cage=AB2_AB$Cage,Population=AB2_AB$Population,Time=AB2_AB$Time,AB2AB_S.pca$x)
    AB2BA_S.pca.df<-data.frame(Cage=AB2_BA$Cage,Population=AB2_BA$Population,Time=AB2_BA$Time,AB2BA_S.pca$x)
    
    AB2_S.pca.df<-data.frame(Cage=ABpca2$Cage,Population=ABpca2$Population,Time=ABpca2$Time,AB2_S.pca$x)
    
    ## manova for phenoed on apple / bloom sig variables
    
    AB2_S.pca.AABA.df<-rbind(AB2AA_S.pca.df, AB2BA_S.pca.df)
    AB2_S.AABA.PCA.aov<-manova(cbind(PC1,PC2,PC3,PC4,PC5)~Time*Population,AB2_S.pca.AABA.df)
    summary(AB2_S.AABA.PCA.aov) ### SIG INTERCEPT
    
    AB2_S.pca.BBAB.df<-rbind(AB2BB_S.pca.df,AB2AB_S.pca.df)
    AB2_S.BBAB.PCA.aov<-manova(cbind(PC1,PC2,PC3,PC4,PC5)~Time*Population, AB2_S.pca.BBAB.df)
    summary(AB2_S.BBAB.PCA.aov) ###SIG INTERCEPT
    
   
    
    #### GETTING VECTORS TO A DATA FRAME####
    ##All Vars 
    
    AB2AA.pca.vectors<-as.data.frame(AB2AA.pca$rotation)
    AB2AA.pca.vectors$ABBRV<-rownames(AB2AA.pca.vectors)
    total <- merge(AB2AA.pca.vectors,AB2SigTrait3List,by="ABBRV")
   

    AB2BB.pca.vectors<-as.data.frame(AB2BB.pca$rotation)
    AB2BB.pca.vectors$ABBRV<-rownames(AB2BB.pca.vectors)

    AB2AB.pca.vectors<-as.data.frame(AB2AB.pca$rotation)
    AB2AB.pca.vectors$ABBRV<-rownames(AB2AB.pca.vectors)

    
    AB2BA.pca.vectors<-as.data.frame(AB2BA.pca$rotation)
    AB2BA.pca.vectors$ABBRV<-rownames(AB2BA.pca.vectors)
    
    

    ##SIG VARS ##
    
    AB2AA_S.pca.vectors<-as.data.frame(AB2AA_S.pca$rotation)
    AB2AA_S.pca.vectors$ABBRV<-rownames(AB2AA_S.pca.vectors)
    AB2AA_S.pca.vectors <- merge(AB2AA_S.pca.vectors,AB2SigTrait3List,by="ABBRV")

    
    AB2BB_S.pca.vectors<-as.data.frame(AB2BB_S.pca$rotation)
    AB2BB_S.pca.vectors$ABBRV<-rownames(AB2BB_S.pca.vectors)
    AB2BB_S.pca.vectors <- merge(AB2BB_S.pca.vectors,AB2SigTrait3List,by="ABBRV")

    AB2AB_S.pca.vectors<-as.data.frame(AB2AB_S.pca$rotation)
    AB2AB_S.pca.vectors$ABBRV<-rownames(AB2AB_S.pca.vectors)
    AB2AB_S.pca.vectors <- merge(AB2AB_S.pca.vectors,AB2SigTrait3List,by="ABBRV")

    AB2BA_S.pca.vectors<-as.data.frame(AB2BA_S.pca$rotation)
    AB2BA_S.pca.vectors$ABBRV<-rownames(AB2BA_S.pca.vectors)
    AB2BA_S.pca.vectors <- merge(AB2BA_S.pca.vectors,AB2SigTrait3List,by="ABBRV")

    # CALCULATING VAR EXPLAINED
    ##All Vars 
    AB2AA.pca.var.explained<-AB2AA.pca$sdev^2/sum(AB2AA.pca$sdev^2)
    AB2BB.pca.var.explained<-AB2BB.pca$sdev^2/sum(AB2BB.pca$sdev^2)
    AB2AB.pca.var.explained<-AB2AB.pca$sdev^2/sum(AB2AB.pca$sdev^2)
    AB2BA.pca.var.explained<-AB2BA.pca$sdev^2/sum(AB2BA.pca$sdev^2)
    
    AB2.pca.var.explained<-AB2.pca$sdev^2/sum(AB2.pca$sdev^2)
    
    
    
    ##Sig Vars 
    AB2AA_S.pca.var.explained<-AB2AA_S.pca$sdev^2/sum(AB2AA_S.pca$sdev^2)
    AB2BB_S.pca.var.explained<-AB2BB_S.pca$sdev^2/sum(AB2BB_S.pca$sdev^2)
    AB2AB_S.pca.var.explained<-AB2AB_S.pca$sdev^2/sum(AB2AB_S.pca$sdev^2)
    AB2BA_S.pca.var.explained<-AB2BA_S.pca$sdev^2/sum(AB2BA_S.pca$sdev^2)
    
    AB2_S.pca.var.explained<-AB2_S.pca$sdev^2/sum(AB2_S.pca$sdev^2)
    
    
    ####calculating hulls####
    # SUBSETTING AND CALCULATING HULLS FOR EACH Population FOR PLOTTING
    AB2AA.0.pca.df <-subset(AB2AA.pca.df,Time=="0")
    AB2AA.0.pca.hull<- AB2AA.0.pca.df %>% ddply("Time",find_hull)
    View( AB2BA.4.pca.hull)
    
    AB2AA.2.pca.df <-subset(AB2AA.pca.df,Time=="2")
    AB2AA.2.pca.hull<- AB2AA.2.pca.df %>% ddply("Time",find_hull)
    
    AB2AA.4.pca.df <-subset(AB2AA.pca.df,Time=="4")
    AB2AA.4.pca.hull<- AB2AA.4.pca.df %>% ddply("Time",find_hull)
 
    
    AB2BB.0.pca.df <-subset(AB2BB.pca.df,Time=="0")
    AB2BB.0.pca.hull<- AB2BB.0.pca.df %>% ddply("Time",find_hull)
    
    AB2BB.2.pca.df <-subset(AB2BB.pca.df,Time=="2")
    AB2BB.2.pca.hull<- AB2BB.2.pca.df %>% ddply("Time",find_hull)
    
    AB2BB.4.pca.df <-subset( AB2BB.pca.df,Time=="4")
    AB2BB.4.pca.hull<-  AB2BB.4.pca.df %>% ddply("Time",find_hull)
    
    
    AB2AB.0.pca.df <-subset(AB2AB.pca.df,Time=="0")
    AB2AB.0.pca.hull<- AB2AB.0.pca.df %>% ddply("Time",find_hull)
    
    AB2AB.2.pca.df <-subset(AB2AB.pca.df,Time=="2")
    AB2AB.2.pca.hull<- AB2AB.2.pca.df %>% ddply("Time",find_hull)
    
    AB2AB.4.pca.df <-subset(AB2AB.pca.df,Time=="4")
    AB2AB.4.pca.hull<- AB2AB.4.pca.df %>% ddply("Time",find_hull)
    
    
    AB2BA.0.pca.df <-subset(AB2BA.pca.df,Time=="0")
    AB2BA.0.pca.hull<- AB2BA.0.pca.df %>% ddply("Time",find_hull)
    
    AB2BA.2.pca.df <-subset(AB2BA.pca.df,Time=="2")
    AB2BA.2.pca.hull<- AB2BA.2.pca.df %>% ddply("Time",find_hull)
    
    AB2BA.4.pca.df <-subset(AB2BA.pca.df,Time=="4")
    AB2BA.4.pca.hull<-  AB2BA.4.pca.df %>% ddply("Time",find_hull)
    
    #### HUlls using just sig variables
    
    AB2AA_S.0.pca.df <-subset(AB2AA_S.pca.df,Time=="0")
    AB2AA_S.0.pca.hull<- AB2AA_S.0.pca.df %>% ddply("Time",find_hull)
    
    
    AB2AA_S.2.pca.df <-subset(AB2AA_S.pca.df,Time=="2")
    AB2AA_S.2.pca.hull<- AB2AA_S.2.pca.df %>% ddply("Time",find_hull)
    
    AB2AA_S.4.pca.df <-subset(AB2AA_S.pca.df,Time=="4")
    AB2AA_S.4.pca.hull<- AB2AA_S.4.pca.df %>% ddply("Time",find_hull)
    
    
    AB2BB_S.0.pca.df <-subset(AB2BB_S.pca.df,Time=="0")
    AB2BB_S.0.pca.hull<- AB2BB_S.0.pca.df %>% ddply("Time",find_hull)
    
    AB2BB_S.2.pca.df <-subset(AB2BB_S.pca.df,Time=="2")
    AB2BB_S.2.pca.hull<- AB2BB_S.2.pca.df %>% ddply("Time",find_hull)
    
    AB2BB_S.4.pca.df <-subset(AB2BB_S.pca.df,Time=="4")
    AB2BB_S.4.pca.hull<-  AB2BB_S.4.pca.df %>% ddply("Time",find_hull)
    
    
    AB2AB_S.0.pca.df <-subset(AB2AB_S.pca.df,Time=="0")
    AB2AB_S.0.pca.hull<- AB2AB_S.0.pca.df %>% ddply("Time",find_hull)
    
    AB2AB_S.2.pca.df <-subset(AB2AB_S.pca.df,Time=="2")
    AB2AB_S.2.pca.hull<- AB2AB_S.2.pca.df %>% ddply("Time",find_hull)
    
    AB2AB_S.4.pca.df <-subset(AB2AB_S.pca.df,Time=="4")
    AB2AB_S.4.pca.hull<- AB2AB_S.4.pca.df %>% ddply("Time",find_hull)
    
    
    AB2BA_S.0.pca.df <-subset(AB2BA_S.pca.df,Time=="0")
    AB2BA_S.0.pca.hull<- AB2BA_S.0.pca.df %>% ddply("Time",find_hull)
    
    AB2BA_S.2.pca.df <-subset(AB2BA_S.pca.df,Time=="2")
    AB2BA_S.2.pca.hull<- AB2BA_S.2.pca.df %>% ddply("Time",find_hull)
    
    AB2BA_S.4.pca.df <-subset( AB2BA_S.pca.df,Time=="4")
    AB2BA_S.4.pca.hull<-  AB2BA_S.4.pca.df %>% ddply("Time",find_hull)
    
    
  }
  
  #### PLOT Config####
  {
    # SIZES
    {
      F3_POINTSIZE=2
      F3_LINESIZE=1
      F3_FILLALPHA=0.3
      
      F3_XLLIM=-4.3
      F3_XULIM=3.15
      F3A_YLLIM=-3.0
      F3A_YULIM=4.45
      F3B_YLLIM=-3.45
      F3B_YULIM=4.0
      
      F3_LGDJUST=c(0,1)
      F3_LGDPOS="none" ###chnage if want legends on indivisual plots to c(0,0)
      F3_LGDSPCX=0.02
      F3_LGDSPCY=0.1
      F3_LGDKEYSIZE=0.6
      
      F3B_ARROWSIZE=.5
      F3_ARROWHEAD=0.2
      F3B_TEXTSIZE=2.5
      F3_VECTORSCALEX=1.2
      F3_VECTORSCALEY=1.2
      F3_VECTOR_ALPHA=0.4
      F3_CIRCLE_COLOR="grey90"
      
      F3_VECTORS_POS_L=0.0
      F3_VECTORS_POS_B=0.62
      F3_VECTORS_POS_T=1.02
      F3_VECTORS_POS_R=0.38
      
      F4_VECTORS_POS_L=1.0
      F4_VECTORS_POS_B=0.62
      F4_VECTORS_POS_T=1.02
      F4_VECTORS_POS_R=1.38
      
      
      F3_PDFHW=7.1
      F3_PDFHH=3.5
      F3_PDFVW=3.5
      F3_PDFVH=7.1
      
      # NEEDED FOR THE PCA VECTORS PLOT
      angle <- seq(-pi, pi, length = 50) 
      circle <- data.frame(x = sin(angle), y = cos(angle))
    }
    
    
    #### PLOT with all variables####
    #####On Bloomoington ####
    
    AB2OnBloom<-ggplot()+
      scale_color_manual(values = AB2ColorListPopulation) +
      scale_shape_manual(values = AB2TimeShape) +
      scale_fill_manual(values=AB2ColorListTime)+
      
      geom_polygon(data=AB2BB.0.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BB.0.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2BB.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BB.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2BB.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BB.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      
      geom_polygon(data=AB2AB.0.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AB.0.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2AB.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AB.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2AB.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AB.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      
      geom_point(AB2BB.0.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2BB.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2BB.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      geom_point(AB2AB.0.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2AB.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2AB.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      
    scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3A_YLLIM, F3A_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(AB2.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("PC2 -",percent(AB2.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("Apple and Bloomington, Phenotyped on Bloomington, All Traits") +
      theme(legend.justification=F3_LGDJUST, 
            legend.position = F3_LGDPOS,
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.spacing.x = unit(F3_LGDSPCX,"cm"),
            legend.spacing.y = unit(F3_LGDSPCY,"cm"),
            legend.key.size = unit(F3_LGDKEYSIZE,"cm"))
    
    AB2OnBloom
    
    ##### On Apple ####
    
    AB2OnApple<-ggplot()+
      scale_color_manual(values = AB2ColorListPopulation) +
      scale_shape_manual(values = AB2TimeShape) +
      scale_fill_manual(values=AB2ColorListTime)+
      
      
      geom_polygon(data=AB2AA.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AA.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2AA.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AA.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      
      geom_polygon(data=AB2BA.0.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BA.0.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2BA.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BA.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2BA.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BA.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      
      
      geom_point(AB2AA.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2AA.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      geom_point(AB2BA.0.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2BA.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2BA.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3A_YLLIM, F3A_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(AB2.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("PC2 -",percent(AB2.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("Apple and Bloomington, Phenotyped on Apple, All Traits") +
      theme(legend.justification=F3_LGDJUST, 
            legend.position = F3_LGDPOS,
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.spacing.x = unit(F3_LGDSPCX,"cm"),
            legend.spacing.y = unit(F3_LGDSPCY,"cm"),
            legend.key.size = unit(F3_LGDKEYSIZE,"cm"))
    
    AB2OnApple
    
    #####Vectors####
    
    
    Vec_BB<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data=AB2BB.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2, color=ABBRV), alpha=F3_VECTOR_ALPHA,
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=F3B_ARROWSIZE) + 
      geom_text(data=AB2BB.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*PC2,label=str_wrap(ABBRV, 2),color=ABBRV ),
                size=F3B_TEXTSIZE) + scale_color_manual(values = ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "transparent",colour = "transparent"))
    
    Vec_AB<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data=AB2AB.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2, color=ABBRV), alpha=F3_VECTOR_ALPHA,
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=F3B_ARROWSIZE) + 
      geom_text(data=AB2AB.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*PC2,label=str_wrap(ABBRV, 2),color=ABBRV ),
                size=F3B_TEXTSIZE) + scale_color_manual(values = ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "transparent",colour = "transparent"))
    
    Vec_AA<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data=AB2AA.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2, color=ABBRV), alpha=F3_VECTOR_ALPHA,
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=F3B_ARROWSIZE) + 
      geom_text(data=AB2AA.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*PC2,label=str_wrap(ABBRV, 2),color=ABBRV ),
                size=F3B_TEXTSIZE) + scale_color_manual(values = ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "transparent",colour = "transparent"))
    
    Vec_BA<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data=AB2BA.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2, color=ABBRV), alpha=F3_VECTOR_ALPHA,
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=F3B_ARROWSIZE) + 
      geom_text(data=AB2BA.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*PC2,label=str_wrap(ABBRV, 2),color=ABBRV ),
                size=F3B_TEXTSIZE) + scale_color_manual(values = ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "transparent",colour = "transparent"))
    
    ####PLOT with just sig vars ####
    
    ##### ON Bloomington####
    AB2SOnBloom<-ggplot()+
      scale_color_manual(values = AB2ColorListPopulation) +
      scale_shape_manual(values = AB2TimeShape) +
      scale_fill_manual(values = AB2ColorListTime)+
      
      geom_polygon(data=AB2BB_S.0.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BB_S.0.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2BB_S.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BB_S.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2BB_S.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BB_S.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      
      geom_polygon(data=AB2AB_S.0.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AB_S.0.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2AB_S.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AB_S.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2AB_S.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AB_S.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_point(AB2BB_S.0.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2BB_S.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2BB_S.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      geom_point(AB2AB_S.0.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2AB_S.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2AB_S.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      
      
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3B_YLLIM, F3B_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(AB2_S.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("-PC2 -",percent(AB2_S.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("C") +
      theme(legend.justification=F3_LGDJUST, 
            legend.position = F3_LGDPOS,
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.spacing.y = unit(F3_LGDSPCY,"cm"),
            legend.key.size = unit(F3_LGDKEYSIZE,"cm")
      )
    
    AB2SOnBloom
    
    #####ON Apple####
    
    AB2SOnApple<-ggplot()+
      scale_color_manual(values = AB2ColorListPopulation) +
      scale_shape_manual(values = AB2TimeShape) +
      scale_fill_manual(values = AB2ColorListTime)+
      
      geom_polygon(data=AB2AA_S.0.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AA_S.0.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2AA_S.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AA_S.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2AA_S.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AA_S.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
    
      geom_polygon(data=AB2BA_S.0.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BA_S.0.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2BA_S.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BA_S.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB2BA_S.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BA_S.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_point(AB2AA_S.0.pca.df,
                  mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2AA_S.2.pca.df,
                  mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2AA_S.4.pca.df,
                  mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      geom_point(AB2BA_S.0.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2BA_S.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB2BA_S.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      
      
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3B_YLLIM, F3B_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(AB2_S.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("-PC2 -",percent(AB2_S.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("F") +
      theme(legend.justification=F3_LGDJUST, 
            legend.position = F3_LGDPOS,
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.spacing.y = unit(F3_LGDSPCY,"cm"),
            legend.key.size = unit(F3_LGDKEYSIZE,"cm")
      )
    
    AB2SOnApple
View(AB2BB_S.pca.vectors)
    #####Vectors####
    VecS_BB<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data= AB2BB_S.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2,),
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=.5) + 
      geom_text(data= AB2BB_S.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*PC2,label=str_wrap(LONG),color=ABBRV ),
                size=2) + scale_color_manual(values = AB2ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "transparent",colour = "transparent"))
  VecS_BB
    VecS_AB<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data= AB2AB_S.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2, color="black"),
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=F3B_ARROWSIZE) + 
      geom_text(data= AB2AB_S.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*PC2,label=str_wrap(LONG,1),color=ABBRV ),
                size=2) + scale_color_manual(values = AB2ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "transparent",colour = "transparent"))
    VecS_AB
    
    VecS_AA<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data= AB2AA_S.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2,),
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=F3B_ARROWSIZE) + 
      geom_text(data=AB2AA_S.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*PC2,label=str_wrap(LONG, 2),color=ABBRV ),
                size=2) + scale_color_manual(values = AB2ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "transparent",colour = "transparent"))
    VecS_AA
    
    VecS_BA<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data= AB2BA_S.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2,),
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=F3B_ARROWSIZE) + 
      geom_text(data= AB2BA_S.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*PC2,label=str_wrap(LONG),color=ABBRV ),
                size=2) + scale_color_manual(values = AB2ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "transparent",colour = "transparent"))
    VecS_BA
    

    ####DUMMY PLOT TO MAKE LEGEND####
    
    Legend2_PLOT<-ggplot(AB2.pca.AABA.df,aes(x=PC1,y=PC2,fill=Time,color=Population,shape=Time))+
      geom_point(AB2.pca.AABA.df,
                 mapping = aes(x=PC1,y=PC2,fill=Time,shape=Time,color=Population),size=F3_POINTSIZE)+
      geom_polygon(data=AB2BB.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BB.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      geom_polygon(data=AB2AA.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AA.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      geom_polygon(data=AB2BA.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2BA.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      geom_polygon(data=AB2AB.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=AB2AB.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      scale_color_manual(values =AB2ColorListPopulation, breaks=c("aa","ab","ba","bb"), labels = c("LQ","LQ (Assayed in HQ)","HQ (Assayed in LQ)","LQ" )) +
      scale_shape_manual(guide="legend", values = TimeShape, breaks = c(0,2,4), labels = c("Foudner","TP2","TP4")) +
      scale_fill_manual(guide="legend", values=ColorListTime, breaks = c(0,2,4), labels = c("Founder","TP2","TP4")) +
      theme(legend.position = "bottom")
    Legend2_PLOT
    
    legend2 <- get_legend(Legend2_PLOT)        
    grid.newpage()                               
    grid.draw(legend2)  
    

    ####Compound Plots####
    
    Vec_BBAB <- ggarrange(Vec_BB,Vec_AB, ncol=2, labels = c("Bloomington Pops","Apple Pops"), font.label = list(size = 10))
    Vec_BBAB
    
    AB2OnBloomVEC <-ggarrange(Vec_BBAB, AB2OnBloom, ncol = 1 , heights = c(.3,1))
    
    Vec_AABA <-ggarrange(Vec_AA, Vec_BA, ncol=2, labels =c("Apple Pops", "Bloomington Pops"), font.label = list(size=8))
    
    AB2OnAppleVEC<-ggarrange(Vec_AABA, AB2OnApple, ncol=1, heights = c(.3,1))
    
    AB2ALLTRAITS<-ggarrange (AB2OnBloomVEC, AB2OnAppleVEC, legend2, ncol=3, widths=c(1,1,.2))
    AB2ALLTRAITS 
    
    
    VecS_BBAB <- ggarrange(VecS_BB,VecS_AB, ncol=2, labels = c("A ","B"), font.label = list(size = 8))

    AB2SOnBloomVEC <-ggarrange(VecS_BBAB, AB2SOnBloom, ncol = 1, heights = c(.5,1))
 
    VecS_AABA <-ggarrange(VecS_AA, VecS_BA, ncol=2, labels =c("D", "E"), font.label = list(size=8))
    
    AB2SOnAppleVEC<-ggarrange(VecS_AABA, AB2SOnApple, ncol=1, heights = c(.,1))
    
    AB2SIGTRAIT<-ggarrange (AB2SOnBloomVEC, AB2SOnAppleVEC, ncol=2, widths=c(1,1,.2))
   
    AB2SIGTRAITS <-ggarrange(AB2SIGTRAIT,legend2, ncol=1, heights = c(1,.1))  
    AB2SIGTRAITS
    