
library(survival)  
library(ggplot2)
library(survival) 
library(ggfortify)
library(ggplot2)
library(reshape2)
library(survminer)
library(contsurvplot)
library(riskRegression)
library(pammtools)

METAB<-read.csv("D:/Dropbox/0Survival/Code/ML/Metabric/Data/DataF.csv" ) 
METAB<-METAB[METAB$PAM50=="LumA",]
 
# Subset data for the specified time period
METAB <- METAB[with(METAB, RFSurvival >= time_period_start & RFSurvival <= time_period_end),]

METAB$Score <- METAB$Macrophage 

breaks <- quantile(METAB$Score, probs = c(0,  2/3, 1))

# Create a new column based on quantile breaks
METAB$Score_Label <- cut(METAB$Score, breaks = breaks, labels = c("Low/Int",  "High"), include.lowest = TRUE)
 
METAB2 <- METAB 
METAB2 <- na.omit(METAB2)

# Fit the survival model
model_fit <- survfit(Surv(RFSurvival, RFStatus) ~ Score_Label, data = METAB2)

# Perform a log-rank test to compare survival curves
surv_diff <- survdiff(Surv(RFSurvival, RFStatus) ~ Score_Label, data = METAB2)

# Extract the p-value from the log-rank test
p_value <- 1 - pchisq(surv_diff$chisq, df = 1)

# Print the p-value
cat("Log-rank test p-value:", p_value, "\n")
 
model_fit<-survfit(Surv(RFSurvival,RFStatus)~Score_Label,data=METAB2)
 
pdf("D:/Dropbox (UiO)/0Survival/Code/Analysis/Survival/Macro - LumA/Fig6gMetabTL.pdf",width=8,height=5)
autoplot(model_fit) + 
  labs(x="\n Survival Time (Months)",y="Survival Probabilities\n"  ) +theme_bw() + 
  annotate("text", x = 0, y = 0.1, label = paste(
    "Log-rank test p-value =", round(p_value, 4)), hjust = 0, vjust = 0, size = 5)+
    theme(text = element_text(size =18)) 
dev.off()

 


  