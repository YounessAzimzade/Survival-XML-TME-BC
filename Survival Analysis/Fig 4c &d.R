# Analysis of Macrophage Fractions and pCR Rate and RoR Scores in LumA Subtype
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(EnvStats)
library(ggpubr)
library(rstatix)

# Read the data files
NAC <- read.csv("~NAC/AllSamples.csv")
NAC2 <- read.csv("~NAC/All2000.csv")

# Match and add RoRScore and RoRRisk to NAC data
NAC$RoRScore <- NAC2$RoRScore[match(NAC$Mixture, NAC2$Mixture)]
NAC$RoRRisk  <- NAC2$RoRRisk[match(NAC$Mixture, NAC2$Mixture)]

# Remove NAC2 to free up memory
remove(NAC2)

# Filter data for LumA subtype and ER positive samples
NAC <- NAC[NAC$PAM50 == "LumA", ]
NAC <- NAC[NAC$ER == "Positive", ]

# Set Macrophage score
NAC$Score <- NAC$Macrophage  

# Define breaks using quantile on the Score
breaks <- quantile(NAC$Score, probs = c(0, 2/3, 1))

# Create a new column based on quantile breaks
NAC$Score_Label <- cut(NAC$Score, breaks = breaks, labels = c("Low/Int", "High"), include.lowest = TRUE)

# Create a contingency table for Response and Score_Label
my_table <- as.data.frame(table(NAC$Response, NAC$Score_Label))
my_table$`pCR Rate` <- 3 * my_table$Freq / nrow(NAC)
my_table <- subset(my_table, my_table$Var1 == "pCR")
my_table2 <- my_table[, c(2, 4)]
my_table2[1, 2] <- my_table2[1, 2] * 0.5

# Rename columns for clarity
colnames(my_table2) <- c("Macrophage Fraction", "pCR Rate")

# Convert pCR Rate to numeric
my_table2$`pCR Rate` <- as.numeric(my_table2$`pCR Rate`)

# Plot pCR Rate by Macrophage Fraction
pdf("~NAC/Fig4c.pdf", width = 3.5, height = 6)
ggplot(my_table2, aes(x = `Macrophage Fraction`, y = `pCR Rate`, fill = `Macrophage Fraction`)) +
  geom_bar(stat = "identity", width = 0.6) + 
  xlab("Macrophage Fraction") +
  ylab("pCR Rate") + 
  theme_bw() + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 20))
dev.off()

# Melt the contingency table for further analysis
my_table.m <- melt(my_table)

# Create a contingency table for RoRRisk and Score_Label
my_table2 <- as.data.frame(table(NAC$RoRRisk, NAC$Score_Label))

# Melt the NAC data for plotting
AllPat.m <- melt(NAC[, c(43, 46)])
AllPat.m <- na.omit(AllPat.m)
AllPat.m <- AllPat.m[, -2]
colnames(AllPat.m) <- c("Macrophage_Fraction", "RoR_Score")

# Perform t-test and adjust p-values
stat.test <- AllPat.m %>% 
  t_test(`RoR_Score` ~ `Macrophage_Fraction`) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test <- stat.test %>%
  add_xy_position(fun = "mean_se", x = "Macrophage_Fraction") 

# Set dodge position for violin and box plots
dodge <- position_dodge(width = 0.8)

# Plot RoR Score by Macrophage Fraction
pdf("~NAC/Fig4d.pdf", width = 3.5, height = 6)
ggplot(AllPat.m, aes(y = `RoR_Score`, x = `Macrophage_Fraction`, color = `Macrophage_Fraction`)) +
  geom_violin(position = dodge) + 
  geom_boxplot(width = 0.3, alpha = 0.2, position = dodge) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = TRUE, size = 6) +   
  theme_bw() + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = "none")
dev.off()
