# Correlation between Macrophages subsets and T_Ex and T_Reg in scRNA seq
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024

# Load necessary libraries
library(reshape2)
library(ggplot2)

# Load T cells metadata
# Replace the file path with the actual path where your T cells metadata file is located
T_cells_MetaData <- read.delim("D:/Dropbox/scRNA/MetadataTcells.txt")

# Load Macrophages metadata
# Replace the file path with the actual path where your Macrophages metadata file is located
Macs_MetaData <- read.delim("D:/Dropbox/scRNA/MetadataMacs.txt")

# Replace dashes with dots in the Cell_ID column of Macs_MetaData
Macs_MetaData$Cell_ID <- gsub("-", ".", Macs_MetaData$Cell_ID)

# Load scRNA counts for macrophages
# Replace the file path with the actual path where your scRNA counts file is located
Counts_Mac <- read.delim(file="D:/Dropbox/scRNA/scMacsAll.txt")
Counts_Mac <- as.data.frame(Counts_Mac)

# Define the list of markers
marklist <- c("HLA-A", "HLA-B", "HLA-C")

# Extract data for the specified markers and transpose it
Macs_HLA <- data.frame(t(subset(Counts_Mac, rownames(Counts_Mac) %in% marklist)))

# Convert the data to numeric
Macs_HLA_num <- data.frame(lapply(Macs_HLA, as.numeric))

# Apply log2 transformation
Macs_HLA_num <- log2(Macs_HLA_num + 1)

# Calculate the sum of HLA markers for each cell
Macs_HLA_num$Sum <- rowSums(Macs_HLA_num)

# Add cell IDs as a column
Macs_HLA_num$ID <- rownames(Macs_HLA)

# Find the median HLA-ABC level
median <- median(Macs_HLA_num$Sum)

# Subset high and low HLA-ABC expressing macrophages based on the median
Macs_HLA_H <- subset(Macs_HLA_num, Sum > median)
Macs_HLA_L <- subset(Macs_HLA_num, Sum < median)

# Subset macrophages metadata for high and low HLA-ABC expressing cells
MacsHigh <- subset(Macs_MetaData, Cell_ID %in% Macs_HLA_H$ID)
MacsLow <- subset(Macs_MetaData, Cell_ID %in% Macs_HLA_L$ID)

# Subset T cells metadata for a specific T cell subtype
Target_T_cell <- subset(T_cells_MetaData, (celltype_subset == "T_cells_c8_CD8+_LAG3"))
### to include T_Reg, set celltype_subset == "T_cells_c2_CD4+_T-regs_FOXP3" 


# Combine high HLA-ABC macrophages with target T cells
AllH <- rbind(MacsHigh[,1:8], Target_T_cell[,1:8])

# Combine low HLA-ABC macrophages with target T cells
AllL <- rbind(MacsLow[,1:8], Target_T_cell[,1:8])

# Create frequency tables for high and low HLA-ABC groups
Freq_HLA_hi <- as.data.frame(table(AllH$orig.ident, AllH$celltype_minor))
Freq_HLA_hi <- dcast(Freq_HLA_hi, Var1 ~ Var2)

Freq_HLA_lo <- as.data.frame(table(AllL$orig.ident, AllL$celltype_minor))
Freq_HLA_lo <- dcast(Freq_HLA_lo, Var1 ~ Var2)


# Create your plot and add the p-value as an annotation
# Perform the correlation test
cor_test_result <- cor.test( Freq_HLA_lo$`Macrophage` , Freq_HLA_lo$`T cells CD8+`  )

# Extract the p-value and format it in scientific notation
p_value_scientific <- formatC(cor_test_result$p.value, format = "e", digits = 2)
cor_value_formatted <- formatC(cor_test_result$estimate, format = "f", digits = 2)




pdf(" ~/Fig5k.pdf",width=6,height=5 )
ggplot(Freq_HLA_lo,aes(x= `Macrophage` ,y= `T cells CD8+` ))+
  scale_x_log10()+scale_y_log10()+
  geom_point(size=1,alpha=1)+geom_smooth(method="lm", level=0.95)+  
 # facet_wrap(~ Feature)+   
 theme_bw()+ #ylim(-0.16,0.16)+xlim(-0.01,0.95)+
  theme(text=element_text(size=19))+xlab("Frequency of HLA-ABC^{lo} macrophages")+ylab("Frequency of T_{Ex}")+
  annotate("text", x =20, y = 400, label = paste("Pearson r:", cor_value_formatted, "\np-value:", p_value_scientific),   size = 6)
dev.off()




cor_test_result <- cor.test(Freq_HLA_hi$`Macrophage`+1 , Freq_HLA_hi$`T cells CD8+`)

# Extract the p-value and format it in scientific notation
p_value_scientific <- formatC(cor_test_result$p.value, format = "e", digits = 2)
cor_value_formatted <- formatC(cor_test_result$estimate, format = "f", digits = 2)


pdf("~/Fig5j.pdf",width=6,height=5 )
ggplot(Freq_HLA_hi,aes(x=`Macrophage`,y=`T cells CD8+`))+
  scale_x_log10()+scale_y_log10()+
  geom_point(size=1,alpha=1)+geom_smooth(method="lm", level=0.95)+  
  # facet_wrap(~ Feature)+   
  theme_bw()+ #ylim(-0.16,0.16)+xlim(-0.01,0.95)+
  theme(text=element_text(size=19))+xlab("Frequency of HLA-ABC^{hi} macrophages")+ylab("Frequency of T_{Ex}")+
  annotate("text", x =20, y = 400, label = paste("Pearson r:", cor_value_formatted, "\np-value:", p_value_scientific),   size = 6)
dev.off()

 
   