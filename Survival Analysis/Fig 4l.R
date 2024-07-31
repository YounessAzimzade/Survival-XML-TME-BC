# Analysis of Macrophage Distance and Survival in LumA and LumB Subtypes
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024

# Load necessary libraries
library(survival)
library(ggplot2)
library(reshape2)
library(readr)
library(stringi)
library(survminer)
library(contsurvplot)
library(riskRegression)

# Define the length change variable
lngCh <- 3.0

# Read and process the minimum distances file for LumA
raw_text <- readLines("~LumA/Min_distancesAll.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
closest_distancesAll <- read_csv(combined_text)
closest_distancesAll <- as.data.frame(closest_distancesAll)
closest_distancesAll$ImNum <- sub("-.*$", "", closest_distancesAll[, 1])

# Read and process the LumA samples file
raw_text <- readLines("~/LumA/Samples3L.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
Samples <- read_csv(combined_text)

# Recode Score_Label factor levels
Samples$Score_Label <- recode_factor(Samples$Score_Label, "Low" = "Low/Int", "Int" = "Low/Int")

# Define the cell type pattern to match
celltx <- list("MHC I & II^{hi}")

# Escape special characters for pattern matching
celltx <- gsub("\\^", "\\\\\\^", celltx)
celltx <- gsub("\\{", "\\\\\\{", celltx)
celltx <- gsub("\\}", "\\\\\\}", celltx)
celltx <- gsub("\\+", "\\\\\\+", celltx)
cellx <- celltx

# Initialize a dataframe to store coefficients
Coefs <- data.frame(matrix(ncol = 6))
colnames(Coefs) <- c("radi", "close", "nclose", "dist", "ndist", "data")

# Loop through different radii to calculate coefficients for LumA
for (i in 3:200) {
  idis <- i * lngCh
  
  Coefs[i, 1] <- idis 
  distanceSub <- as.data.frame(closest_distancesAll[, grep(cellx, colnames(closest_distancesAll))])
  
  rownames(distanceSub) <- closest_distancesAll[, 1]
  
  # Define breaks using fixed values and the current radius
  breaks <- c(5.230785, idis, 707.867107)
  
  # Create a new column based on quantile breaks
  distanceSub$Score_Label <- cut(distanceSub[, 1], breaks = breaks, labels = c("Low", "High"), include.lowest = TRUE)
  
  distanceSub$ID <- sub("-.*$", "", rownames(distanceSub))
  
  SubHigh <- subset(distanceSub, Score_Label == "High")
  SubHigh2 <- as.data.frame(table(SubHigh$ID))
  
  SubLow <- subset(distanceSub, Score_Label == "Low")
  SubLow2 <- as.data.frame(table(SubLow$ID))
  
  Samples$Close_Mac <- SubLow2$Freq[match(Samples$Mixture, SubLow2$Var1)]
  Samples$Dist_Mac <- SubHigh2$Freq[match(Samples$Mixture, SubHigh2$Var1)]
  
  Samples$Close_Mac <- Samples$Close_Mac / (max(Samples$Close_Mac, na.rm = TRUE))
  Samples$Dist_Mac <- Samples$Dist_Mac / (max(Samples$Dist_Mac, na.rm = TRUE))
  
  model <- coxph(Surv(RFSurvival, RFStatus) ~ Close_Mac, data = Samples, x = TRUE)
  Coefs[i, 2] <- model$coefficients
  Coefs[i, 3] <- nrow(SubLow2)
  
  model <- coxph(Surv(RFSurvival, RFStatus) ~ Dist_Mac, data = Samples, x = TRUE)
  Coefs[i, 4] <- model$coefficients
  Coefs[i, 5] <- nrow(SubHigh2)
}

CoefsLumA <- Coefs

# Repeat the process for LumB

# Read and process the minimum distances file for LumB
raw_text <- readLines("~LumB/Min_distancesAll.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
closest_distancesAll <- read_csv(combined_text)
closest_distancesAll <- as.data.frame(closest_distancesAll)
closest_distancesAll$ImNum <- sub("-.*$", "", closest_distancesAll[, 1])

# Read and process the LumB samples file
raw_text <- readLines("~/METABRIC/AllInc.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
Samples <- read_csv(combined_text)
Samples <- subset(Samples, Samples$PAM50 == "LumB")

# Reuse the cell type pattern and escape special characters
celltx <- list("MHC I & II^{hi}")
celltx <- gsub("\\^", "\\\\\\^", celltx)
celltx <- gsub("\\{", "\\\\\\{", celltx)
celltx <- gsub("\\}", "\\\\\\}", celltx)
celltx <- gsub("\\+", "\\\\\\+", celltx)
cellx <- celltx

# Initialize a dataframe to store coefficients
Coefs <- data.frame(matrix(ncol = 6))
colnames(Coefs) <- c("radi", "close", "nclose", "dist", "ndist", "data")

# Loop through different radii to calculate coefficients for LumB
for (i in 3:175) {
  idis <- i * lngCh
  
  Coefs[i, 1] <- idis 
  distanceSub <- as.data.frame(closest_distancesAll[, grep(cellx, colnames(closest_distancesAll))])
  
  rownames(distanceSub) <- closest_distancesAll[, 1]
  
  # Define breaks using fixed values and the current radius
  breaks <- c(5.230785, idis, 707.867107)
  
  # Create a new column based on quantile breaks
  distanceSub$Score_Label <- cut(distanceSub[, 1], breaks = breaks, labels = c("Low", "High"), include.lowest = TRUE)
  
  distanceSub$ID <- sub("-.*$", "", rownames(distanceSub))
  
  SubHigh <- subset(distanceSub, Score_Label == "High")
  SubHigh2 <- as.data.frame(table(SubHigh$ID))
  
  SubLow <- subset(distanceSub, Score_Label == "Low")
  SubLow2 <- as.data.frame(table(SubLow$ID))
  
  Samples$Close_Mac <- SubLow2$Freq[match(Samples$ImageNumber, SubLow2$Var1)]
  Samples$Dist_Mac <- SubHigh2$Freq[match(Samples$ImageNumber, SubHigh2$Var1)]
  
  Samples$Close_Mac <- Samples$Close_Mac / (max(Samples$Close_Mac, na.rm = TRUE))
  Samples$Dist_Mac <- Samples$Dist_Mac / (max(Samples$Dist_Mac, na.rm = TRUE))
  
  model <- coxph(Surv(RFSurvival, RFStatus) ~ Close_Mac, data = Samples, x = TRUE)
  Coefs[i, 2] <- model$coefficients
  Coefs[i, 3] <- nrow(SubLow2)
  
  model <- coxph(Surv(RFSurvival, RFStatus) ~ Dist_Mac, data = Samples, x = TRUE)
  Coefs[i, 4] <- model$coefficients
  Coefs[i, 5] <- nrow(SubHigh2)
}

CoefsLumB <- Coefs

# Combine LumA and LumB coefficients
CoefsLumA$data <- "LumA"
CoefsLumB$data <- "LumB"
Coefs <- na.omit(rbind(CoefsLumA, CoefsLumB))

# Plot the combined coefficients
pdf("~/Fig4l.pdf", width = 8, height = 5)
ggplot(Coefs, aes(x = radi, y = -dist, col = data)) +
  geom_point(size = 3) + 
  theme_bw() +
  geom_line() +
  xlim(0, 175) +
  ylim(-1.6, 1.95006) +
  scale_color_manual(values = c("darkblue", "cyan")) +
  xlab("r_T") +
  ylab("SC for macrophages with r_m > r_T") +
  theme(text = element_text(size = 18))
dev.off()
