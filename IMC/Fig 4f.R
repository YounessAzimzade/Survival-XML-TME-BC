# Survival Analysis for LumA Subtype Using Spatial Omics Data
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024

# This script performs survival analysis for the LumA subtype of breast cancer
# using spatial omics data. It calculates survival probabilities and performs
# a log-rank test to compare survival curves based on macrophage fractions

# Load necessary libraries
library(survival)
library(ggfortify)
library(ggplot2)
library(reshape2)
library(survminer)
library(contsurvplot)
library(riskRegression)
library(pammtools)
library(stringi)
library(readr)

# Read and encode the text file
raw_text <- readLines("~/Spatial Omics/Samples_All.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
Samples <- read_csv(combined_text)

# Subset data for the LumA subtype
Samples <- subset(Samples, Samples$PAM50 == "LumA")

# Use macrophage scores for analysis
Samples$Score <- Samples$Macrophages

# Calculate quantile breaks for macrophage scores
breaks <- quantile(Samples$Score, probs = c(0, 1/3, 2/3, 1))

# Create a new column based on quantile breaks
Samples$Score_Label <- cut(Samples$Score, breaks = breaks, labels = c("Low/Int", "Low/Int", "High"), include.lowest = TRUE)

# Fit the survival model
model_fit <- survfit(Surv(RFSurvival, RFStatus) ~ Score_Label, data = Samples)

# Perform a log-rank test to compare survival curves
surv_diff <- survdiff(Surv(RFSurvival, RFStatus) ~ Score_Label, data = Samples)

# Extract the p-value from the log-rank test
p_value <- 1 - pchisq(surv_diff$chisq, df = 1)

# Print the p-value
cat("Log-rank test p-value:", p_value, "\n")

# Fit the survival model again for plotting
model_fit <- survfit(Surv(RFSurvival, RFStatus) ~ Score_Label, data = Samples)

# Save the survival plot as a PDF file
pdf("~/LumA/Fig 4f.pdf", width = 8, height = 5)
autoplot(model_fit) + 
  labs(x = "\n Survival Time (Months)", y = "Survival Probabilities\n") + theme_bw() + 
  annotate("text", x = 0, y = 0.1, label = paste("Log-rank test p-value =", round(p_value, 4)), hjust = 0, vjust = 0, size = 6) +
  theme(text = element_text(size = 16))
dev.off()


write.csv(Samples,file="~/LumA/Samples3L.csv"
          ,sep="\t", row.names = F,col.names = TRUE,quote = FALSE)