# Analysis of T_Reg & T_Ex  Frequencies at Different Distances from HLA-ABC^{lo} and HLA-ABC^{hi} Macrophages
# in LumA samples
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024

# This script analyzes the frequencies of different cell types at varying distances 
# from macrophages with high MHC I & II expression in the LumA subtype of breast cancer. 
# It visualizes the frequency distribution of T_{Reg} & T_{Ex} cells at different distances.

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(tidyverse)
library(Matrix)
library(proxy)
library(readr)
library(stringi)

# Read and encode the text file for minimum distances
raw_text <- readLines("~/LumA/Min_Distance_All.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
closest_distancesAll <- read_csv(combined_text)
closest_distancesAll <- as.data.frame(closest_distancesAll)

# Read and encode the text file for sample data
raw_text <- readLines("D:/DropBox/SpTr/METABRIC/AllInc.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
Samples <- read_csv(combined_text)
Samples <- subset(Samples, PAM50 == "LumA")

# Extract minimum distance for MHC I & II high cells
distanceSub <- as.data.frame(closest_distancesAll[, "MHC I & II^{hi}"])
rownames(distanceSub) <- closest_distancesAll[, 1]
colnames(distanceSub) <- "Min"

filtered_values <- distanceSub$Min[is.finite(distanceSub$Min)]

# Define breaks using quantiles on the filtered values
breaks <- quantile(filtered_values, probs = c(0, 1/2, 1))

# Create a new column based on quantile breaks
distanceSub$Score_Label <- cut(distanceSub$Min, breaks = breaks, labels = c("close", "dist"), include.lowest = TRUE)

# Label the closest distances
closest_distancesAll$Label <- distanceSub$Score_Label[match(rownames(distanceSub), closest_distancesAll$ID)]

# Subset data based on proximity to MHC I & II high cells
CloseCD <- subset(closest_distancesAll, Label == "close")
DistCD <- subset(closest_distancesAll, Label == "dist")

# Initialize data frames for cell counts
LenDis <- 10
celltx <- colnames(closest_distancesAll[, 2:33])
Cell_CountsC <- data.frame(matrix(NA, ncol = length(celltx), nrow = 200))
colnames(Cell_CountsC) <- celltx
Cell_CountsC <- cbind(Data = "Close to MHC I & II^{hi}", Cell_CountsC)

Cell_CountsD <- data.frame(matrix(NA, ncol = length(celltx), nrow = 200))
colnames(Cell_CountsD) <- celltx
Cell_CountsD <- cbind(Data = "Distant from MHC I & II^{hi}", Cell_CountsD)

# Loop through each cell type and calculate frequencies at different distances
for (i in celltx) {
  print(i)
  cellx <- i
  
  distanceSub <- as.data.frame(closest_distancesAll[[cellx]])
  distanceSubC <- as.data.frame(CloseCD[[cellx]])
  distanceSubD <- as.data.frame(DistCD[[cellx]])
  
  filtered_valuesC <- as.data.frame(distanceSubC[[1]][is.finite(distanceSubC[[1]])])
  filtered_valuesD <- as.data.frame(distanceSubD[[1]][is.finite(distanceSubD[[1]])])
  
  if (nrow(filtered_valuesC) > 0 && nrow(filtered_valuesD) > 0) {
    for (j in 0:100) {
      distan = (-0.5 + j) * LenDis
      Cell_CountsC[j + 1, "Dist"] <- distan
      Cell_CountsD[j + 1, "Dist"] <- distan
      compC <- as.data.frame(filtered_valuesC[filtered_valuesC[, 1] < distan & filtered_valuesC[, 1] > (j - 1) * LenDis, ])
      FreqC <- nrow(compC)
      Cell_CountsC[j + 1, cellx] <- FreqC
      
      compD <- as.data.frame(filtered_valuesD[filtered_valuesD[, 1] < distan & filtered_valuesD[, 1] > (j - 1) * LenDis, ])
      FreqD <- nrow(compD)
      Cell_CountsD[j + 1, cellx] <- FreqD
    }
  }
}

# Combine cell count data
Cell_Counts <- rbind(Cell_CountsC, Cell_CountsD)
Cell_Counts.m <- melt(Cell_Counts[, -2], id.vars = c("Dist", "Data"))

# Save the plot as a PDF file
pdf("~/LumA/Fig5c.pdf", width = 5, height = 5)
ggplot(subset(Cell_Counts.m, variable == "T_{Reg} & T_{Ex}"), aes(x = Dist, y = value, col = Data)) +
  geom_point(size = 6) +
  scale_x_log10() +
  scale_color_manual(values = c("purple", "yellow2")) +
  geom_line(size = 2) +
  theme_bw() +
  xlab(expression("Minimum distance from Macrophages (" * mu * "m)")) +
  ylab("Frequency of T_{Reg} & T_{Ex}") +
  labs(col = "Macrophages") +  # Change the legend title
  theme(text = element_text(size = 16), legend.position = "top") +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1))
dev.off()
