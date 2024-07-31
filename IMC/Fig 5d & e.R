library(ggplot2)
library(reshape2)  
library(readr)
library(stringi)
library(ggpubr)
library(rstatix)
library(dplyr) 

raw_text<-readLines("~/METABRIC/Frequencies.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")
# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")
# Read the data from the text connection
Samples  <- read_csv(combined_text)

Samples <- na.omit(Samples) 


SingleCells<-read.csv("~/METABRIC/SingleCells.csv") 

SingleCells2  <-  SingleCells 
SingleCells   <-  SingleCells2 
SingleCells  <- subset(SingleCells, SingleCells$cellPhenotype=="Macrophages")

SingleCells$ID <- paste(SingleCells$ImageNumber,SingleCells$ObjectNumber,sep="-")
SingleCells$ID <- paste(SingleCells$ID,SingleCells$cellPhenotype,"_",sep="_")

SingleCells$metabric_id <- gsub('-','.',SingleCells$metabric_id) 
SingleCells <- subset(SingleCells, SingleCells$ImageNumber  %in% Samples$Mixture)


 
 
SingleCells$Score <- SingleCells$HLA.ABC



breaks <- quantile(SingleCells$Score, probs = c(0, 1/2, 1))
 
SingleCells$Score_Label <- cut(SingleCells$Score, breaks = breaks, labels = c("Low","High"), include.lowest = TRUE)
 
SubHigh <- subset(SingleCells, SingleCells$Score_Label=="High")
SubHigh2  <- as.data.frame(table(SubHigh$ImageNumber))

SubLow <- subset(SingleCells, SingleCells$Score_Label=="Low")
SubLow2  <- as.data.frame(table(SubLow$ImageNumber))

Samples$`HLA.ABC^{lo}` <- SubLow2$Freq[match(Samples$Mixture, SubLow2$Var1)]
Samples$`HLA.ABC^{hi}` <- SubHigh2$Freq[match(Samples$Mixture, SubHigh2$Var1)]

Samples[is.na(Samples)] <- 0

 

cor_test_result <- cor.test( (Samples$`T_{Reg} & T_{Ex}` ),  (Samples$`HLA.ABC^{hi}` ))

# Extract the p-value and format it in scientific notation
p_value_scientific <- formatC(cor_test_result$p.value, format = "e", digits = 2)
cor_value_formatted <- formatC(cor_test_result$estimate, format = "f", digits = 2)


# Open a PDF device
pdf("~/Fig5e.pdf",width=6,height=5)
# Create the plot and add the p-value annotation
ggplot(Samples, aes(x = `HLA.ABC^{hi}`  , y =  `T_{Reg} & T_{Ex}` )) +
  geom_point(size = 1, alpha = 1) +
  scale_x_log10() + 
  scale_y_log10()+
  geom_smooth(method = "lm", level = 0.95) +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Frequency of HLA-ABC^{hi} macrophages") +
  ylab("Frequency of T_{Reg} & T_{Ex}") +
   annotate("text", x =10, y = 400, label = paste("Pearson r:", cor_value_formatted, "\np-value:", p_value_scientific),   size = 6)
dev.off()


# Perform the correlation test
cor_test_result <- cor.test( Samples$`T_{Reg} & T_{Ex}`   ,  (Samples$`HLA.ABC^{lo}`  ))

# Extract the p-value and format it in scientific notation
p_value_scientific <- formatC(cor_test_result$p.value, format = "e", digits = 2)
cor_value_formatted <- formatC(cor_test_result$estimate, format = "f", digits = 2)


# Open a PDF device
pdf("~/Fig5d.pdf",width=6,height=5)
# Create the plot and add the p-value annotation
ggplot(Samples, aes(x=`HLA.ABC^{lo}`, y=`T_{Reg} & T_{Ex}` )) +
  geom_point(size = 1,alpha=1) +
  scale_x_log10() + 
  scale_y_log10()+
   geom_smooth(method="lm",level=0.95) +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Frequency of HLA-ABC^{lo} macrophages") +
  ylab("Frequency of T_{Reg} & T_{Ex}") +
  annotate("text", x =10, y = 400, label = paste("Pearson r:", cor_value_formatted, "\np-value:", p_value_scientific),   size = 6)
# Close the PDF device
dev.off()


 