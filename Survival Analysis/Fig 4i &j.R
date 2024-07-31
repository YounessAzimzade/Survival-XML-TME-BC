# Analysis of association of two groups of Macrophages (divided by their minimum distance to MHC I II)
# and RFS in LumA Subtype  
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(survival)  
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

# Define the cell types of interest
celltx <- list("MHC I & II^{hi}") 

# Escape special characters for pattern matching
celltx <- gsub("\\^", "\\\\\\^", celltx)      
celltx <- gsub("\\{", "\\\\\\{", celltx)  
celltx <- gsub("\\}", "\\\\\\}", celltx)  
celltx <- gsub("\\+", "\\\\\\+", celltx)

# Subset data based on the specified cell types
cellx <-  celltx 
distanceSub <- as.data.frame(matrix(ncol = 0, nrow = nrow(closest_distancesAll)))

for (pattern in cellx) {
  matched_cols <- grep(pattern, colnames(closest_distancesAll), value = TRUE)
  distanceSub <- cbind(distanceSub, closest_distancesAll[, matched_cols])
}

# Set row names and calculate minimum distances
rownames(distanceSub) <-  closest_distancesAll[,1]
distanceSub$Min <- apply(distanceSub, 1, function(row) min(row, na.rm = TRUE))

# Filter finite values for break calculation
filtered_values <- distanceSub$Min[is.finite(distanceSub$Min)]

# Define breaks using quantile on the filtered values
breaks <- quantile(filtered_values, probs = c(0, 1/2, 1))

# Create a new column based on quantile breaks
distanceSub$Score_Label <- cut(distanceSub$Min, breaks = breaks, labels = c("Low",  "High"), include.lowest = TRUE)

# Extract ID from row names
distanceSub$ID <-  sub("-.*$", "", rownames(distanceSub))

# Subset high and low groups
SubHigh <- subset(distanceSub, Score_Label == "High")
SubHigh2 <- as.data.frame(table(SubHigh$ID))

SubLow <- subset(distanceSub, Score_Label == "Low")
SubLow2 <- as.data.frame(table(SubLow$ID))

# Add frequency counts to Samples data
Samples$`Close to MHC I & II^{hi}` <- SubLow2$Freq[match(Samples$ImageNumber, SubLow2$Var1)]
Samples$`Distant from MHC I & II^{hi}` <- SubHigh2$Freq[match(Samples$ImageNumber, SubHigh2$Var1)]

# Replace NA values with 0
Samples[is.na(Samples)] <- 0

# Generate survival plot for 'Close to MHC I & II^{hi}'
pdf("~/LumA/Fig4i.pdf", width = 7, height = 5)
model <- coxph(Surv(RFSurvival, RFStatus) ~ `Close to MHC I & II ^{hi}`, data = Samples, x = TRUE)

p <- plot_surv_area(time = "RFSurvival",
                    status = "RFStatus",
                    variable = "Close to MHC I & II ^{hi}",
                    data = Samples,
                    model = model,
                    start_color = "white", end_color = "purple") + 
  ylim(0, 1) +  
  labs(x = "\n Survival Time (Months)", y = "Survival Probabilities\n") +
  theme(legend.position = c(0.05, 0.05), # Position at bottom left
        legend.justification = c(0, 0), # Anchor point of legend
        legend.box.just = "left", # Justify the legend box
        legend.margin = margin(0, 0, 0, 0), # Remove margin around legend
        legend.box.margin = margin(0, 0, 0, 0)) + # Remove margin around legend box
  theme(text = element_text(size = 18))

print(p)
dev.off()

# Generate survival plot for 'Distant from MHC I & II^{hi}'
pdf("~/LumA/Fig4j.pdf", width = 7, height = 5)
model <- coxph(Surv(RFSurvival, RFStatus) ~ `Distant from MHC I & II ^{hi}`, data = Samples, x = TRUE)

p <- plot_surv_area(time = "RFSurvival",
                    status = "RFStatus",
                    variable = "Distant from MHC I & II ^{hi}",
                    data = Samples,
                    model = model,
                    start_color = "white", end_color = "yellow2") + 
  ylim(0, 1) +  
  labs(x = "\n Survival Time (Months)", y = "Survival Probabilities\n") +
  theme(legend.position = c(0.05, 0.05), # Position at bottom left
        legend.justification = c(0, 0), # Anchor point of legend
        legend.box.just = "left", # Justify the legend box
        legend.margin = margin(0, 0, 0, 0), # Remove margin around legend
        legend.box.margin = margin(0, 0, 0, 0)) + # Remove margin around legend box
  theme(text = element_text(size = 18))

print(p)
dev.off()
