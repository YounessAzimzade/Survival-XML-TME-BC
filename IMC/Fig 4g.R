# Analysis of Minimum Distance Between Macrophages and Other Cell Types in LumA Subtype two subgroups defined by 
# fraction of macrophages
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024

 

 
# Load necessary libraries
 library(ggplot2)
library(reshape2) 
library(readr)
library(stringi)
library(survminer)
library(contsurvplot)
library(riskRegression)

# Read and encode the text file for minimum distances
raw_text <- readLines("~/LumA/Min_Distance_All.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
Min_Distance <- read_csv(combined_text)
Min_Distance <- as.data.frame(Min_Distance)

# Extract image number from the data
Min_Distance$Image_Number <- sub("-.*$", "", Min_Distance[, 1])

# Read and encode the text file for LumA samples with macrophage levels
raw_text <- readLines("~/LumA/Samples3L.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
Samples <- read_csv(combined_text)

# Match and label the samples based on macrophage scores
Min_Distance$Label <- Samples$Score_Label[match(Min_Distance$Image_Number, Samples$Image_Number)]

# Melt the data for plotting
Min_Dist.m <- melt(Min_Distance[, -c(1, 2, 34)])

# Rename columns for clarity
colnames(Min_Dist.m) <- c("Macrophage Fraction", "CellType", "Distance")

# Process macrophage fractions
Min_Dist.m$Macrophage <- as.character(Min_Dist.m$Macrophage)
Min_Dist.m$Macrophage <- factor(Min_Dist.m$Macrophage, levels = c("High", "Low/Int"))

# Remove infinite distances and omit NA values
Min_Dist.m$Distance[Min_Dist.m$Distance == "Inf"] <- NA
Min_Dist.m <- na.omit(Min_Dist.m)

# Load additional libraries for statistical analysis
library(ggpubr)
library(rstatix)
library(dplyr)

# Perform t-tests and adjust p-values for multiple comparisons
stat.test <- Min_Dist.m %>%
  group_by(CellType) %>%
  t_test(Distance ~ Macrophage) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") 

# Add xy positions for significant results
stat.test <- stat.test %>%
  add_xy_position(fun = "mean_se", x = "Macrophage") 

# Order and filter significant results
stat.test <- stat.test[order(stat.test$p),]
stat.test <- stat.test[stat.test$p.adj < 0.05,]
stat.test2 <- stat.test

# Calculate log-transformed p-values with direction
stat.test2$statistic2 <- stat.test2$statistic / abs(stat.test2$statistic)
stat.test2$`-Log(p.adj)*Direction` <- -log10(stat.test2$p.adj) * stat.test2$statistic2
stat.test2 <- stat.test2[order(-stat.test2$`-Log(p.adj)*Direction`),]

# Process cell types for plotting
stat.test2$CellType <- as.character(stat.test2$CellType)
stat.test2$CellType <- factor(stat.test2$CellType, levels = stat.test2$CellType)

# Assign macrophage fraction labels
stat.test2$"Macrophage Fraction" <- 0
stat.test2$"Macrophage Fraction"[stat.test2$statistic2 == 1] <- "High"
stat.test2$"Macrophage Fraction"[stat.test2$statistic2 == -1] <- "Low/Int"

# Summarize average distances and standard errors
Min_Dist.m2 <- Min_Dist.m %>%
  group_by(`Macrophage Fraction`, CellType) %>%
  summarize(
    AverageDistance = mean(Distance, na.rm = TRUE),
    StdDevDistance = sd(Distance, na.rm = TRUE),
    SEM = sd(Distance, na.rm = TRUE) / sqrt(n())
  )

# Filter data for significant cell types
Min_Dist.m3 <- subset(Min_Dist.m2, Min_Dist.m2$CellType %in% stat.test2$CellType)

# Process cell types for plotting
Min_Dist.m3$CellType <- as.character(Min_Dist.m3$CellType)
Min_Dist.m3$CellType <- factor(Min_Dist.m3$CellType, levels = stat.test2$CellType)

# Process macrophage fraction labels for plotting
Min_Dist.m3$`Macrophage Fraction` <- as.character(Min_Dist.m3$`Macrophage Fraction`)
Min_Dist.m3$`Macrophage Fraction` <- factor(Min_Dist.m3$`Macrophage Fraction`, levels = c("Low/Int", "High"))

# Save the plot as a PDF file
pdf("~/LumA/Fig4g.pdf", width = 12, height = 12)
ggplot(Min_Dist.m3, aes(x = CellType, y = AverageDistance, color = `Macrophage Fraction`, group = `Macrophage Fraction`)) +
  geom_point(size = 6) +
  geom_line(size = 2) +
  geom_errorbar(aes(ymin = AverageDistance - 1.95 * SEM, ymax = AverageDistance + 1.95 * SEM), width = 0.3) +
  theme_bw() +
  theme(text = element_text(size = 35), axis.text.x = element_text(angle = 90), legend.position = "none") +
  labs(x = "Cell Type", y = "Average Minimum Distance")
dev.off()
