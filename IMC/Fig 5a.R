# Analysis of Marker Levels on Macrophages Close to and Distant from MHC I & II
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024

# This script analyzes the levels of markers on macrophages that are close to and distant from 
# cells with high MHC I & II expression in the LumA subtype of breast cancer. It performs statistical 
# tests and visualizes the differences in marker levels.

# Load necessary libraries
library(ggpubr)
library(rstatix)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(readr)
library(stringi)

# Read and encode the text file for sample data
raw_text <- readLines("D:/Dropbox/SpTr/METABRIC/AllInc.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
Samples <- read_csv(combined_text)

# Subset data for the LumA subtype
Samples <- subset(Samples, Samples$PAM50 == "LumA")

# Read and encode the text file for minimum distances
raw_text <- readLines("~/LumA/Min_Distance_All.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
closest_distancesAll <- read_csv(combined_text)
closest_distancesAll <- as.data.frame(closest_distancesAll)

# Define a new score to classify macrophages
closest_distancesAll$Score2 <- closest_distancesAll$`MHC I & II^{hi}`

filtered_values <- closest_distancesAll$Score2[is.finite(closest_distancesAll$Score2)]

breaks <- quantile(filtered_values, probs = c(0, 1/2, 1))

# Create a new column based on quantile breaks
closest_distancesAll$Score_Label2 <- cut(closest_distancesAll$Score2, breaks = breaks, labels = c("Low", "High"), include.lowest = TRUE)

# Load single cell data from spatial omics
SingleCells <- read.csv("D:/Dropbox/SpTr/METABRIC/SingleCells.csv")
SingleCells <- subset(SingleCells, SingleCells$cellPhenotype == "Macrophages")

# Create unique IDs for single cells
SingleCells$ID <- paste(SingleCells$ImageNumber, SingleCells$ObjectNumber, sep = "-")
SingleCells$ID <- paste(SingleCells$ID, SingleCells$cellPhenotype, "_", sep = "_")

SingleCells$metabric_id <- gsub('-', '.', SingleCells$metabric_id)
SingleCells <- subset(SingleCells, SingleCells$ImageNumber %in% Samples$ImageNumber)

# Keep only marker columns
Marker_Levels <- SingleCells[, c(12:50, 54)]

# Combine DNA1 and DNA2 columns
Marker_Levels$DNA <- Marker_Levels$DNA1 + Marker_Levels$DNA2
Marker_Levels <- Marker_Levels[, -38:-39]

# Match and label the single cells based on macrophage scores
Marker_Levels$Label2 <- closest_distancesAll$Score_Label2[match(closest_distancesAll$ID, Marker_Levels$ID)]

Marker_Levels <- na.omit(Marker_Levels)

# Melt the data for plotting
Marker_Levels.m <- melt(Marker_Levels[, -38])

# Rename columns for clarity
colnames(Marker_Levels.m) <- c("Label", "Feature", "Value")

# Process macrophage labels
Marker_Levels.m$Label <- as.character(Marker_Levels.m$Label)
Marker_Levels.m$Label <- factor(Marker_Levels.m$Label, levels = c("Low", "High"))

# Perform t-tests and adjust p-values for multiple comparisons
stat.test <- Marker_Levels.m %>%
  group_by(Feature) %>%
  t_test(Value ~ Label) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

# Filter significant results
stat.test <- stat.test[stat.test$p.adj < 0.05,]

# Calculate log-transformed p-values with direction
stat.test$statistic2 <- stat.test$statistic / abs(stat.test$statistic)
stat.test$`-Log(p.adj)*Direction` <- -log10(stat.test$p.adj) * stat.test$statistic2
stat.test <- stat.test[order(-stat.test$`-Log(p.adj)*Direction`),]

# Process features for plotting
stat.test$Feature <- as.character(stat.test$Feature)
stat.test$Feature <- factor(stat.test$Feature, levels = stat.test$Feature)

# Assign macrophage proximity labels
stat.test$Macrophage <- 0
stat.test$Macrophage[stat.test$statistic2 == 1] <- "Close to MHC I & II ^{hi}"
stat.test$Macrophage[stat.test$statistic2 == -1] <- "Distant from MHC I & II ^{hi}"

# Save the plot as a PDF file
pdf("~/LumA/Fig5a.pdf", width = 16, height = 8)
ggplot(stat.test, aes(x = Feature, y = `-Log(p.adj)*Direction`, fill = Macrophage)) +
  geom_bar(stat = "identity", width = 0.605) +
  scale_fill_manual(values = c("purple", "yellow2")) +
  scale_color_discrete(name = "") +
  theme(axis.title.y = element_blank()) +
  xlab("Marker") +
  theme_bw() +
  theme(text = element_text(size = 30), axis.text.x = element_text(angle = 90)) 
dev.off()
