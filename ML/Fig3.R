library(pheatmap) 
library(ggplot2)
library(reshape2)  
library(ggpubr) 
  
MBRC<-read.csv("~/METABRIC/Survival Score.csv")
MBRC <- cbind("Survival (MBRC)", MBRC)
colnames(MBRC)<-c("Data","Cell Type","Slope","LI", "HI", "Sign", "Score")

TCGA<-read.csv("~/TCGA/Survival Score.csv")
TCGA <- cbind("Survival (TCGA)", TCGA) 
colnames(TCGA)<-c("Data","Cell Type","Slope","LI", "HI", "Sign", "Score")

pCR <-  read.csv("~/pCR Score.csv")
pCR <- cbind("pCR", pCR) 
colnames(pCR)<-c("Data","Cell Type","Slope","LI", "HI", "Sign", "Score")


All_Scores <- rbind(MBRC, TCGA, pCR)
colnames(All_Scores)<-c("Data","Cell Type","Slope","LI", "HI", "Sign", "Survival Score")

All_ScoresM <- dcast(All_Scores[,-3:-6] , `Cell Type`~ Data )
All_ScoresM$pCR[is.na(All_ScoresM$pCR)] <- 0
All_ScoresM <- na.omit(All_ScoresM)
All_ScoresM <- All_ScoresM[All_ScoresM$`Survival (MBRC)`*All_ScoresM$`Survival (TCGA)`>0,]

All_Scores2 <- All_ScoresM

All_ScoresM$Ave <- All_ScoresM$`Survival (MBRC)`+All_ScoresM$`Survival (TCGA)`
All_ScoresM <- All_ScoresM[order(-All_ScoresM$Ave),]

All_Scores2$`Cell Type` <- as.character(All_Scores2$`Cell Type`)
All_Scores2$`Cell Type` <- factor(All_Scores2$`Cell Type`,levels= All_ScoresM$`Cell Type` )

All_Scores2.m <- melt(All_Scores2)
colnames(All_Scores2.m)<-c("Cell Type","Data", "Score")

All_Scores2.m$Data <-  as.character(All_Scores2.m$Data)
All_Scores2.m$Data <- factor(All_Scores2.m$Data,levels= c("Survival (MBRC)", "Survival (TCGA)" ,  "pCR"))

All_Scores2.m$Sign <- sign(All_Scores2.m$Score) 

All_Scores2.m$Sign[All_Scores2.m$Sign==1] <- "Score > 0"
All_Scores2.m$Sign[All_Scores2.m$Sign==-1] <- "Score < 0"
All_Scores2.m$Score[All_Scores2.m$Score==0] <- NA
All_Scores2.m <- na.omit(All_Scores2.m)

pdf("C:/Users/younesa/Dropbox (Personal)/0Survival/Code/Analysis/Sur vs pCR/All/Fig5a.pdf",width=19,height=7.7)
ggplot(All_Scores2.m,aes(`Cell Type`, forcats::fct_rev(as.factor(Data)),fill=as.factor(Sign),size=abs(`Score`)))+
  geom_point(shape=21,stroke=0) +
  # geom_hline(yintercept = seq(.5, 4.5, 1), size = .9) +
  scale_x_discrete(position = "top") +
  scale_radius(range=c(12,27)) + 
  # scale_fill_viridis(low = "#F8766D", high = "#00BFC4", limits = c(-1, 1)) +
  scale_fill_manual(values=c("#ef8a62","#67a9cf"))+
  theme_bw() +theme(text=element_text(size =30))+
  theme(axis.text.x=element_text(angle=90))+
  theme(legend.position = "top", 
        # panel.grid.major = element_blank(),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18))+
  guides(size=guide_legend(override.aes=list(fill = NA,color="black",stroke=.1), 
                           label.position="bottom",
                           title.position="right",order=1),
         fill = guide_legend(override.aes = list(size = 25),
                             title = NULL))+
  labs(size="Area = |Score|" ,fill="Score",x=NULL,y=NULL)
dev.off()  
 
 
write.csv(All_Scores2.m,file="C:/Users/younesa/Dropbox (Personal)/0Survival/Code/Analysis/Sur vs pCR/All/AllFit.csv",sep="\t",
          row.names=FALSE,quote=FALSE)  

