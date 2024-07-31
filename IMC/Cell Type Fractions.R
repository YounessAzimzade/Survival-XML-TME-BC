setwd("D:/Dropbox (UiO)/0Survival/Code/Analysis/Survival") 
library(survival) 
library(ggfortify)
library(ggplot2)
library(reshape2)
library(survminer)
library(contsurvplot)
library(riskRegression)
library(pammtools)

SingleCells<-read.csv("~/METABRIC/IMC/SingleCells.csv") 

Samples <- as.data.frame(table(SingleCells$metabric_id, SingleCells$cellPhenotype))

colnames(Samples) <- c("Mixture", "Var1", "Var2")
Samples <- dcast(Samples, Mixture~ Var1) 
 
for (i in 1:nrow(Samples)) {Samples[i,2:33] <- Samples[i,2:33]/(sum(Samples[i,2:33]))}


Samples$Mixture <- gsub('-','.',Samples$Mixture) 

METAB<-read.csv("D:/Dropbox (UiO)/Datasets/Metabric/Data.csv" ) 
METAB <- METAB[,c (1,38:51)]
METAB <- METAB[,-14]

#METAB <- METAB[METAB$PAM50=="LumA",]


Merged4 <- merge(Samples, METAB, by = "Mixture", all.x = TRUE)

Merged4 <- na.omit(Merged4)


heatmap(cor(Samples[,2:33]))

Merged2 <- na.omit(Merged)

write.csv(Merged4,file="D:/Dropbox (UiO)/SpTr/METABRIC/DataFCLNA.csv",sep="\t",
          row.names=TRUE,quote=FALSE,fileEncoding = "UTF-8") 

write.table(Samples2,file="D:/Dropbox (UiO)/SpTr/METABRIC/Fractions.txt",sep="\t",
          row.names=TRUE,quote=FALSE) 

Samples<-read.csv("D:/Dropbox (UiO)/SpTr/METABRIC/Fractions.csv",sep = ",", fileEncoding = "UTF-8") 
Samples5<-read.delim("D:/Dropbox (UiO)/SpTr/METABRIC/Fractions.txt") 


heatmap(cor(Samples2))
 