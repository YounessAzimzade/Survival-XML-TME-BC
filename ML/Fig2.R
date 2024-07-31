library(ggplot2)
library(reshape2)  
library(ggpubr)  
library(dplyr) 

Data<-read.csv("~/METABRIC/Data 2.csv")
RSF<-read.csv("~/METABRIC/SHAP RSF.csv")
SSVM<-read.csv("~/METABRIC/SHAP SSVM.csv")
 
DataX<- Data[,3:30]
for(i in 1:ncol(DataX)){ 
  DataX[,i]<-DataX[,i]/(quantile(DataX[,i],probs=0.99))}
DataX[DataX>1] <-1

ShapXS<- SSVM[,2:29]
ShapXR<- RSF[,2:29]

ShapS.m <- melt(ShapXS) 
ShapS.m$value <- scale(ShapS.m$value)

ShapR.m <- melt(ShapXR) 
ShapR.m$value <- scale(ShapR.m$value)

Data.m <- melt(DataX) 

ssvm_Data <- cbind(Data.m,ShapS.m$value)
rsf_Data <- cbind(Data.m,ShapR.m$value)
AllPat  <- cbind(Data.m,-1* (ShapR.m$value+ShapS.m$value)/2)

colnames(ssvm_Data) <- c("Feature", "Fraction", "SHAP value") 
colnames(rsf_Data) <- c("Feature", "Fraction", "SHAP value") 
colnames(AllPat) <- c("Feature", "Fraction", "SHAP value") 


fit_ssvm=ssvm_Data%>%group_by(Feature)%>%do(model=lm(`SHAP value` ~ Fraction, data = .)) 
fit_ssvm$Coef<-0
fit_ssvm$CI<-0 
for(i in 1:nrow(fit_ssvm)){ 
  xx <- as.data.frame(confint(fit_ssvm[[2]][[i]], level = 0.999) )
  fit_ssvm[i,3]<-fit_ssvm[[2]][[i]][["coefficients"]][["Fraction"]]
  fit_ssvm[i,4]<-xx[2,1]*xx[2,2]}


fit_rsf=rsf_Data%>%group_by(Feature)%>%do(model=lm(`SHAP value`~Fraction,data=.)) 
fit_rsf$Coef<-0
fit_rsf$CI<-0 
for(i in 1:nrow(fit_rsf)) {
  xx <- as.data.frame(confint(fit_rsf[[2]][[i]],level=0.999))
  fit_rsf[i,3]<-fit_rsf[[2]][[i]][["coefficients"]][["Fraction"]] 
  fit_rsf[i,4]<-xx[2,1]*xx[2,2]}

ssvm_rsf_combined <- as.data.frame(fit_rsf[,1])

ssvm_rsf_combined$Coef1<-fit_ssvm$Coef [match(ssvm_rsf_combined[,1],fit_ssvm$Feature)]
ssvm_rsf_combined$CI1<-fit_ssvm$CI [match(ssvm_rsf_combined[,1],fit_ssvm$Feature)]

ssvm_rsf_combined$Coef2<-fit_rsf$Coef [match(ssvm_rsf_combined[,1],fit_rsf$Feature)]
ssvm_rsf_combined$CI2<-fit_rsf$CI [match(ssvm_rsf_combined[,1],fit_rsf$Feature)]

ssvm_rsf_combined$RS <-ssvm_rsf_combined$Coef1*ssvm_rsf_combined$Coef2 

colnames(ssvm_rsf_combined)<-c("Feature","Coef1","CI1","Coef2","CI2","RS")

ssvm_rsf_combined<-subset(ssvm_rsf_combined,ssvm_rsf_combined$RS>0 &ssvm_rsf_combined$CI1>0 &ssvm_rsf_combined$CI2>0 )

AllPat <- AllPat[ AllPat$Feature %in%ssvm_rsf_combined$Feature,]

fit =AllPat%>%group_by(Feature)%>%do(model=lm(`SHAP value`~Fraction,data=.)) 
fit $Slope <- 0  
fit $LI <- 0 
fit $HI <- 0 
for (i in 1:nrow(fit)) {
  xx <- as.data.frame(confint(fit [[2]][[i]],level=0.999))
  fit [i,3]<-fit[[2]][[i]][["coefficients"]][["Fraction"]] 
  fit [i,4]<-xx[2,1]
  fit [i,5]<-xx[2,2]}

fit  <- fit [,-2]

fit$Sign<-sign(fit$Slope)

fit<-fit[order(-fit$Slope),]

fit$Feature<-as.character(fit$Feature)
fit$Feature<-factor(fit$Feature,levels=fit$Feature)

colnames(fit)<-c("Cell Type","Slope","LI", "HI", "Sign" ) 
 
fit$`Survival Score`<-fit$Slope/(quantile(fit$Slope,probs=0.95)) 
fit$`Survival Score`[fit$`Survival Score`>1] <- 1
fit$`Survival Score`[fit$`Survival Score`<  -1] <- -1 

fit$LI <- fit$LI*(fit$`Survival Score`/fit$Slope) 
fit$HI <- fit$HI*(fit$`Survival Score`/fit$Slope) 

write.csv(fit,file="~/METABRIC/Survival Score.csv",sep="\t",
          row.names=FALSE,quote=FALSE) 
