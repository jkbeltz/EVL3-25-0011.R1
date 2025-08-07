library(tidyverse)
#library(vegan)

View(Phenos.cage)
####SET UP DATA AND RUN GLM####
# FUNTION TO FIND HULLS
#find_hull<-function(df) df[chull(df$PC1,df$PC2),]
AB1TimeList=c("0","1","2","3","4","5")
AB1PopulationList=c("aa","bb")
AB1TraitsList<-data.frame(ABBRV=c("DT","DT_M","DT_F","Vb","SR_M","SR_F","Fc","BS"),LONG=c("Develeopment Time","Development Time (Males)","Development Time (Females)", "Viability","Starvation Resistance (Males)","Starvation Resistance (Females)","Fecundity","Body Size"))
AB1SigTrait3List<-data.frame(ABBRV=c("DT_F","SR_F","Vb","Fc","BS"),LONG=c("Development Time (Females)","Starvation Resistance (Females)","Viability","Fecundity","Body Size"))

# READING THE DATA
ABpca1<-pheno.cage ### this df is produced when generating figure 2, must complete that first. 
ABpca1 <- ABpca1[, -c(8,9,11,12,14,15,17,18,20,21,23,24,25,26,27,29,30,32)]## remove SE and SD
#View(ABpca1)


# SORT THE DATA
ABpca1$Time<-factor(ABpca1$TP,AB1TimeList)
ABpca1$Population=factor(ABpca1$Cage.Pheno,AB1PopulationList)
ABpca1<-ABpca1 %>% arrange(Cage,Time,Population)
View(ABpca1)

# CALCULATING THE PRINCIPAL COMPONENTS
AB1.pca<-prcomp(ABpca1[8:15],scale=TRUE)

# GETTING PCA COORDINATES TO A DATA FRAME
AB1.pca.df<-data.frame(Cage=ABpca1$Cage,Population=ABpca1$Population,Time=ABpca1$Time,AB1.pca$x)

AB1.var.explained<-AB1.pca$sdev^2/sum(AB1.pca$sdev^2)



# GETTING VECTORS TO A DATA FRAME
AB1.pca.vectors<-as.data.frame(AB1.pca$rotation)
AB1.pca.vectors$ABBRV<-rownames(AB1.pca.vectors)
AB1.pca.vectors<-merge(AB1.pca.vectors,AB1TraitsList,by="ABBRV")


# GETTING SUBSETS FOR REARRANGING THE DATAFRAME
AB1.pca.df.0 <- AB1.pca.df %>% subset( Time=="0") 
AB1.pca.df.1 <- AB1.pca.df %>% subset( Time=="1") 
AB1.pca.df.2 <- AB1.pca.df %>% subset( Time=="2") 
AB1.pca.df.3 <- AB1.pca.df %>% subset( Time=="3") 
AB1.pca.df.4 <- AB1.pca.df %>% subset( Time=="4") 
AB1.pca.df.5 <- AB1.pca.df %>% subset( Time=="5") 

#colnames(InEx.pca.df.B)[colnames(InEx.pca.df.B)=="Time"]<-"Selection"

# REARRANGING THE PCA DATA FRAMES
AB1.pca.conv.012345.df<-rbind(AB1.pca.df.0,AB1.pca.df.1,AB1.pca.df.2,AB1.pca.df.3,AB1.pca.df.4,AB1.pca.df.5)

View(AB1.pca.conv.012345.df)
# MANOVA from PCA 

AB1PCA.aov<-manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8)~Time*Population,AB1.pca.conv.012345.df)
summary(AB1PCA.aov)




####PLOTTING
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
#install.packages("ggtext")
library(ggtext)

# PREPARATION
## Setting Working directory


#### LISTS AND COLOR ASSIGNMENTS####
{
  ColorListTime <- c("0"=brewer.pal(n=11, name="Spectral")[6],
                     "1"=brewer.pal(n=11, name="Spectral")[7],
                     "2"=brewer.pal(n=11, name="Spectral")[8],
                     "3"=brewer.pal(n=11, name="Spectral")[9],
                     "4"=brewer.pal(n=11, name="Spectral")[10],
                     "5"=brewer.pal(n=11, name="Spectral")[11])
  
  ColorListTime <- c("0"="grey",
                     "1"="orange",
                     "2"="yellow",
                     "3"="green",
                     "4"= "blue",
                     "5"= "purple")
  
  
  
  ColorListPopulation<-c("aa"="red",
                         "bb"="black")
  
  ColorListMDS=c("aa"="red",
                 "bb"="black")
  
  ColorListTrait<-c(
    "DT"=brewer.pal(n=12, name="Set3")[1],
    "DT_M"=brewer.pal(n=12, name="Set3")[2],
    "DT_F"=brewer.pal(n=12, name="Set3")[3],
    "Vb"=brewer.pal(n=12, name="Set3")[4],
    "SR_M"=brewer.pal(n=12, name="Set3")[5],
    "SR_F"=brewer.pal(n=12,name="Set3")[6],
    "Fc"=brewer.pal(n=12, name="Set3")[7],
    "BS"=brewer.pal(n=12, name="Set3")[8])
  
  display.brewer.pal(8,"Dark2")
  
  AB1TimeList=c("0","1","2","3","4","5")
  AB1PopulationList=c("aa","bb")
  AB1TraitsList<-data.frame(ABBRV=c("DT","DT_M","DT_F","Vb","SR_M","SR_F","Fc","BS"),LONG=c("Develeopment Time","Development Time (Males)","Development Time (Females)", "Viability","Starvation Resistance (Males)","Starvation Resistance (Females)","Fecundity","Body Size"))
  AB1SigTrait3List<-data.frame(ABBRV=c("DT_F","SR_F","Vb","Fc","BS"),LONG=c("Development Time (Females)","Starvation Resistance (Females)","Viability","Fecundity","Body Size"))
  
  SelectionLabels=c("Founder","TP1","TP2","TP3","TP4","TP5" )

  PopulationLabels=c("Apple","Bloomington")
  TimeShape=c("0"=8,"1"=16,"2"=17,"3"=15,"4"=19,"5"=18)
}


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
    
    #### CALCULATING THE PRINCIPAL COMPONENTS####
    AB1_AP<- ABpca1 %>% subset(Population=="aa")
    AB1_BL<- ABpca1 %>% subset(Population=="bb")
    
    AB1_0<- ABpca1 %>% subset(Time=="0")
    AB1_1<- ABpca1 %>% subset(Time=="1")
    AB1_2<- ABpca1 %>% subset(Time=="2")
    AB1_3<- ABpca1 %>% subset(Time=="3")
    AB1_4<- ABpca1 %>% subset(Time=="4")
    AB1_5<- ABpca1 %>% subset(Time=="5")
    
    ### all variables 
    AB1AP.pca<-prcomp(AB1_AP[8:15],scale=TRUE)
    View(AB1_AP)
    AB1BL.pca<-prcomp(AB1_BL[8:15],scale=TRUE)
    AB1.pca<-prcomp(ABpca1[8:15],scale=TRUE)

    #####
    ### sig variables
    AB1AP_S.pca<-prcomp(AB1_AP[c(9,11:14)],scale=TRUE)
    AB1BL_S.pca<-prcomp(AB1_BL[c(9,11:14)],scale=TRUE)
    AB1_S.pca<-prcomp(ABpca1[c(9,11:14)],scale=TRUE)

  
    # GETTING PCA COORDINATES TO A DATA FRAME
    ##All Vars
    AB1AP.pca.df<-data.frame(Cage=AB1_AP$Cage,Population=AB1_AP$Population,Time=AB1_AP$TP,AB1AP.pca$x)
    AB1BL.pca.df<-data.frame(Cage=AB1_BL$Cage,Population=AB1_BL$Population,Time=AB1_BL$TP,AB1BL.pca$x)
    
    AB1.pca.df<-data.frame(Cage=ABpca1$Cage,Population=ABpca1$Population,Time=ABpca1$TP,AB1.pca$x)
    
    ##Sig Vars 
    AB1AP_S.pca.df<-data.frame(Cage=AB1_AP$Cage,Population=AB1_AP$Population,Time=AB1_AP$TP,AB1AP_S.pca$x)
    AB1BL_S.pca.df<-data.frame(Cage=AB1_BL$Cage,Population=AB1_BL$Population,Time=AB1_BL$TP,AB1BL_S.pca$x)
    
    AB1_S.pca.df<-data.frame(Cage=ABpca1$Cage,Population=ABpca1$Population,Time=ABpca1$TP,AB1_S.pca$x)
    
    
    #### GETTING VECTORS TO A DATA FRAME####
    ##All Vars 
    AB1AP.pca.vectors<-as.data.frame(AB1AP.pca$rotation)
    AB1AP.pca.vectors$ABBRV<-rownames(AB1AP.pca.vectors)
    #AB1AP.pca.vectors<-merge(AB1AP.pca.vectors,TraitsList,by="ABBRV")
    View(AB1AP.pca.vectors)
    
    AB1BL.pca.vectors<-as.data.frame(AB1BL.pca$rotation)
    AB1BL.pca.vectors$ABBRV<-rownames(AB1BL.pca.vectors)
    #AB1BL.pca.vectors<-merge(AB1BL.pca.vectors,TraitsList,by="ABBRV")
    
    AB1.pca.vectors<-as.data.frame(AB1.pca$rotation)
    AB1.pca.vectors$ABBRV<-rownames(AB1.pca.vectors)
    #AB1.pca.vectors<-merge(AB1.pca.vectors,TraitsList,by="ABBRV")
    
    
    ##SIG Vars 
    AB1AP_S.pca.vectors<-as.data.frame(AB1AP_S.pca$rotation)
    AB1AP_S.pca.vectors$ABBRV<-rownames(AB1AP_S.pca.vectors)
    #AB1AP_S.pca.vectors<-merge(AB1AP_S.pca.vectors,SigTraitList ,by="ABBRV")
    
    AB1BL_S.pca.vectors<-as.data.frame(AB1BL_S.pca$rotation)
    AB1BL_S.pca.vectors$ABBRV<-rownames(AB1BL_S.pca.vectors)
    #InExLB_S.pca.vectors<-merge(InExLB_S.pca.vectors,SigTrait3List,by="ABBRV")
    
    AB1_S.pca.vectors<-as.data.frame(AB1_S.pca$rotation)
    AB1_S.pca.vectors$ABBRV<-rownames(AB1_S.pca.vectors)
    #InEx_S.pca.vectors<-merge(InEx_S.pca.vectors,SigTrait3List,by="ABBRV")
    
    # CALCULATING VAR EXPLAINED
    ##All Vars 
    AB1AP.pca.var.explained<-AB1AP.pca$sdev^2/sum(AB1AP.pca$sdev^2)
    AB1BL.pca.var.explained<-AB1BL.pca$sdev^2/sum(AB1BL.pca$sdev^2)
    
    AB1.pca.var.explained<-AB1.pca$sdev^2/sum(AB1.pca$sdev^2)
    
    
    
    ##Sig Vars 
    AB1AP_S.pca.var.explained<-AB1AP_S.pca$sdev^2/sum(AB1AP_S.pca$sdev^2)
    AB1BL_S.pca.var.explained<-AB1BL_S.pca$sdev^2/sum(AB1BL_S.pca$sdev^2)
    
    AB1_S.pca.var.explained<-AB1_S.pca$sdev^2/sum(AB1_S.pca$sdev^2)
    
    
    ####calculating hulls####
    # SUBSETTING AND CALCULATING HULLS FOR EACH Population FOR PLOTTING
    AB1AP.0.pca.df <-subset(AB1AP.pca.df,Time=="0")
    AB1AP.0.pca.hull<- AB1AP.0.pca.df %>% ddply("Time",find_hull)
    
    AB1AP.1.pca.df <-subset(AB1AP.pca.df,Time=="1")
    AB1AP.1.pca.hull<- AB1AP.1.pca.df %>% ddply("Time",find_hull)
    
    AB1AP.2.pca.df <-subset(AB1AP.pca.df,Time=="2")
    AB1AP.2.pca.hull<- AB1AP.2.pca.df %>% ddply("Time",find_hull)
    
    AB1AP.3.pca.df <-subset(AB1AP.pca.df,Time=="3")
    AB1AP.3.pca.hull<- AB1AP.3.pca.df %>% ddply("Time",find_hull)
    
    AB1AP.4.pca.df <-subset(AB1AP.pca.df,Time=="4")
    AB1AP.4.pca.hull<- AB1AP.4.pca.df %>% ddply("Time",find_hull)
    
    AB1AP.5.pca.df <-subset(AB1AP.pca.df,Time=="5")
    AB1AP.5.pca.hull<- AB1AP.5.pca.df %>% ddply("Time",find_hull)
    
    
    
    
    AB1BL.0.pca.df <-subset(AB1BL.pca.df,Time=="0")
    AB1BL.0.pca.hull<- AB1BL.0.pca.df %>% ddply("Time",find_hull)
    
    AB1BL.1.pca.df <-subset(AB1BL.pca.df,Time=="1")
    AB1BL.1.pca.hull<- AB1BL.1.pca.df %>% ddply("Time",find_hull)
    
    AB1BL.2.pca.df <-subset(AB1BL.pca.df,Time=="2")
    AB1BL.2.pca.hull<- AB1BL.2.pca.df %>% ddply("Time",find_hull)
    
    AB1BL.3.pca.df <-subset( AB1BL.pca.df,Time=="3")
    AB1BL.3.pca.hull<-  AB1BL.3.pca.df %>% ddply("Time",find_hull)
    
    AB1BL.4.pca.df <-subset( AB1BL.pca.df,Time=="4")
    AB1BL.4.pca.hull<-  AB1BL.4.pca.df %>% ddply("Time",find_hull)
    
    AB1BL.5.pca.df <-subset( AB1BL.pca.df,Time=="5")
    AB1BL.5.pca.hull<-  AB1BL.5.pca.df %>% ddply("Time",find_hull)
    
    
    #### HUlls using just sig variables
    
    AB1AP_S.0.pca.df <-subset(AB1AP_S.pca.df,Time=="0")
    AB1AP_S.0.pca.hull<- AB1AP_S.0.pca.df %>% ddply("Time",find_hull)
    
    AB1AP_S.1.pca.df <-subset(AB1AP_S.pca.df,Time=="1")
    AB1AP_S.1.pca.hull<- AB1AP_S.1.pca.df %>% ddply("Time",find_hull)
    
    AB1AP_S.2.pca.df <-subset(AB1AP_S.pca.df,Time=="2")
    AB1AP_S.2.pca.hull<- AB1AP_S.2.pca.df %>% ddply("Time",find_hull)
    
    AB1AP_S.3.pca.df <-subset(AB1AP_S.pca.df,Time=="3")
    AB1AP_S.3.pca.hull<- AB1AP_S.3.pca.df %>% ddply("Time",find_hull)
    
    AB1AP_S.4.pca.df <-subset(AB1AP_S.pca.df,Time=="4")
    AB1AP_S.4.pca.hull<- AB1AP_S.4.pca.df %>% ddply("Time",find_hull)
    
    AB1AP_S.5.pca.df <-subset(AB1AP_S.pca.df,Time=="5")
    AB1AP_S.5.pca.hull<- AB1AP_S.5.pca.df %>% ddply("Time",find_hull)
    
    
    
    
    AB1BL_S.0.pca.df <-subset(AB1BL_S.pca.df,Time=="0")
    AB1BL_S.0.pca.hull<- AB1BL_S.0.pca.df %>% ddply("Time",find_hull)
    
    AB1BL_S.1.pca.df <-subset(AB1BL_S.pca.df,Time=="1")
    AB1BL_S.1.pca.hull<- AB1BL_S.1.pca.df %>% ddply("Time",find_hull)
    
    AB1BL_S.2.pca.df <-subset(AB1BL_S.pca.df,Time=="2")
    AB1BL_S.2.pca.hull<- AB1BL_S.2.pca.df %>% ddply("Time",find_hull)
    
    AB1BL_S.3.pca.df <-subset( AB1BL_S.pca.df,Time=="3")
    AB1BL_S.3.pca.hull<-  AB1BL_S.3.pca.df %>% ddply("Time",find_hull)
    
    AB1BL_S.4.pca.df <-subset( AB1BL_S.pca.df,Time=="4")
    AB1BL_S.4.pca.hull<-  AB1BL_S.4.pca.df %>% ddply("Time",find_hull)
    
    AB1BL_S.5.pca.df <-subset( AB1BL_S.pca.df,Time=="5")
    AB1BL_S.5.pca.hull<-  AB1BL_S.5.pca.df %>% ddply("Time",find_hull)
    
    
  }
  
  #### PLOT Config####
  {
    # SIZES
    {
      F3_POINTSIZE=3
      F3_LINESIZE=1.5
      F3_FILLALPHA=0.8
      
      F3_XLLIM=-4.3
      F3_XULIM=3.15
      F3A_YLLIM=-3.0
      F3A_YULIM=4.45
      F3B_YLLIM=-3.45
      F3B_YULIM=4.0
      
      F3_LGDJUST=c(0,1)
      F3_LGDPOS=c(0,0)
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
      
      F3_PDFHW=7.1
      F3_PDFHH=3.5
      F3_PDFVW=3.5
      F3_PDFVH=7.1
      
      # NEEDED FOR THE PCA VECTORS PLOT
      angle <- seq(-pi, pi, length = 50) 
      circle <- data.frame(x = sin(angle), y = cos(angle))
    }
    Time
    SelectionLabels[2:3]
    
    #### PLOT with all variables####
  
    # 
    AB1A1<-ggplot()+
      scale_color_manual(values =ColorListPopulation) +
      scale_shape_manual(values = TimeShape) +
      scale_fill_manual(values=ColorListTime)+
      
    geom_polygon(data=AB1BL.0.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL.0.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB1BL.1.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL.1.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB1BL.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB1BL.3.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL.3.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB1BL.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
    
    geom_point(AB1BL.0.pca.df,
               mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB1BL.1.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB1BL.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB1BL.3.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB1BL.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3A_YLLIM, F3A_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(AB1.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("PC2 -",percent(AB1.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("PCA of Phenotypic Space - Bloomington Populations  - All Traits") +
      theme(legend.justification=F3_LGDJUST, 
            legend.position = F3_LGDPOS,
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.spacing.x = unit(F3_LGDSPCX,"cm"),
            legend.spacing.y = unit(F3_LGDSPCY,"cm"),
            legend.key.size = unit(F3_LGDKEYSIZE,"cm"))
    
    AB1A1
    
    p3A2<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data=AB1BL.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = PC2, color=ABBRV), alpha=F3_VECTOR_ALPHA,
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=F3B_ARROWSIZE) + 
      geom_text(data=AB1BL.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*PC2,label=str_wrap(ABBRV, 2),color=ABBRV ),
                size=F3B_TEXTSIZE) + scale_color_manual(values = ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "transparent",colour = "transparent"))
    
    p3A<-AB1A1+inset_element(p3A2,left = F3_VECTORS_POS_L, bottom = F3_VECTORS_POS_B, right = F3_VECTORS_POS_R, top = F3_VECTORS_POS_T)
    p3A
    
    
    ####PLOT with just sig vars ####
    
    AB1B1<-ggplot()+
      scale_color_manual(values =ColorListPopulation) +
      scale_shape_manual(values = TimeShape) +
      scale_fill_manual(values=ColorListTime)+
      
    geom_polygon(data=AB1BL_S.0.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL_S.0.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB1BL_S.1.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL_S.1.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB1BL_S.2.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL_S.2.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB1BL_S.3.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL_S.3.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
      geom_polygon(data=AB1BL_S.4.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Time),color=Population),size=F3_LINESIZE) +
      geom_polygon(data=AB1BL_S.4.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Time))) +
      
    geom_point(AB1BL_S.0.pca.df,
               mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB1BL_S.1.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB1BL_S.2.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB1BL_S.3.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+
      geom_point(AB1BL_S.4.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Time),shape=as.factor(Time),color=Population),size=F3_POINTSIZE)+ 
      
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3B_YLLIM, F3B_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(AB1_S.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("-PC2 -",percent(AB1_S.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("PCA of Phenotypic Space - Bloomington Populations - Sig Traits") +
      theme(legend.justification=F3_LGDJUST, 
            legend.position = F3_LGDPOS,
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.spacing.y = unit(F3_LGDSPCY,"cm"),
            legend.key.size = unit(F3_LGDKEYSIZE,"cm")
      )
    
    AB1B1
    
    p3B2<-ggplot() +  
      coord_fixed() + theme(legend.position = "none") +
      geom_path(aes(x, y), data = circle, colour=F3_CIRCLE_COLOR)+ 
      geom_segment(data=AB1BL_S.pca.vectors, 
                   aes(x = 0, y = 0, xend = PC1, yend = -PC2, color=ABBRV), alpha=F3_VECTOR_ALPHA,
                   arrow = arrow(length = unit(F3_ARROWHEAD, 'picas'),type="closed"),size=F3B_ARROWSIZE) + 
      geom_text(data=AB1BL_S.pca.vectors, fontface = "bold",
                aes(x=F3_VECTORSCALEX*PC1,y=F3_VECTORSCALEY*-PC2,label=str_wrap(ABBRV, 2),color=ABBRV ),
                size=F3B_TEXTSIZE) + scale_color_manual(values = ColorListTrait) +
      theme(plot.background=element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size=AXISTEXTSIZE),
            panel.background = element_rect(fill = "white",colour = "transparent"))  
    
    p3B<-AB1B1+inset_element(p3B2,left = F3_VECTORS_POS_L, bottom = F3_VECTORS_POS_B, right = F3_VECTORS_POS_R, top = F3_VECTORS_POS_T)
    p3B
    p3A
    
    p3H<-plot_grid(p3A,pBL,p3B,ncol=3,rel_widths = c(1,0.05,1),labels = c("A","","B"),label_size = LABELSIZE)
    p3V<-plot_grid(p3A,p3B,ncol=1,rel_heights = c(1,1),labels = c("A","B"),label_size = LABELSIZE,align='v')
  }
  
  # VISUALIZE
  p3H
  p3V
  # SAVE PLOT
  {
    pdf("Figures/Figure3-Phenotypic-PCAPlot-Horizontal.pdf",
        title="Phenotypic PCA Plot",
        width=F3_PDFHW,height = F3_PDFHH)
    p3H
    dev.off()
    
    pdf("Figures/Figure3-Phenotypic-PCAPlot-Vertical.pdf",
        title="Phenotypic PCA Plot",
        width=F3_PDFVW,height = F3_PDFVH)
    p3V
    dev.off()
  }
}
