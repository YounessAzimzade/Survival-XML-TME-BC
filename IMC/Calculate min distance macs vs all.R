# Calculate minimum distance between macrophages and the rest of cell types
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024

# Load necessary libraries
library(Matrix)
library(proxy)
library(readr)
library(stringi)

# Read and encode the text file
raw_text <- readLines("D:/Dropbox (UiO)/SpTr/METABRIC/AllInc.csv", encoding = "UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")

# Combine the text into a single string
combined_text <- paste(correct_text, collapse = "\n")

# Read the data from the text connection
Samples <- read_csv(combined_text)
Samples <- Samples[Samples$PAM50=="LumA",]

METAB$Score <- METAB$Macrophage # - METAB$myCAFs #+  0.462005713* METAB$GenMod2  -0.34896937*METAB$Macrophage   # #  #0.98462520* METAB$GenMod5 #

breaks <- quantile(METAB$Score, probs = c(0, 1/3, 2/3, 1))

# Create a new column based on quantile breaks
METAB$Score_Label <- cut(METAB$Score, breaks = breaks, labels = c("Low/Int", "Low/Int", "High"), include.lowest = TRUE)


# Load spatial omics data including cell types and coordinates
SingleCells <- read.csv("D:/Dropbox (UiO)/SpTr/METABRIC/SingleCells.csv")

# Replace dashes with dots in metabric_id
SingleCells$metabric_id <- gsub('-', '.', SingleCells$metabric_id)

# Subset SingleCells data to include only relevant samples
SingleCells <- subset(SingleCells, SingleCells$ImageNumber %in% Samples$ImageNumber)

# Select relevant columns from SingleCells
SingleCells <- SingleCells[, c(1:11, 51:53)]

# Extract unique cell types
celltypes <- as.data.frame(unique(SingleCells$cellPhenotype))

# Subset cell types to include only macrophages
celltypes2 <- subset(celltypes, celltypes$`unique(SingleCells$cellPhenotype)` == "Macrophages")

# Get the list of sample image numbers
SamplesList <- Samples$ImageNumber

# Create an empty dataframe to store the closest distances for all samples
col_names <- celltypes[, 1]
Min_Distance_All <- data.frame(matrix(NA, ncol = nrow(celltypes), nrow = 0))
colnames(Min_Distance_All) <- col_names

# Loop through each sample in the SamplesList
for (k in SamplesList) {
  print(k)
  # Subset SingleCells data for the current sample
  Sample <- SingleCells[SingleCells$ImageNumber == k, c(2, 4, 12, 13)]
  
  # Modify cellPhenotype to include ObjectNumber and sample ImageNumber
  Sample$cellPhenotype <-paste(Sample$ObjectNumber,Sample$cellPhenotype,"_",sep="_")
  Sample$cellPhenotype <- paste(k,Sample$cellPhenotype,sep="-")
  
  # Extract cell type and coordinates
  AllCells <- Sample[, 2:4]
  colnames(AllCells) <- c("Cell Type", "x", "y")
  
  # Compute the distance matrix
  distance_matrix <- as.matrix(proxy::dist(AllCells[, c("x", "y")]))
  
  # Replace diagonal elements (distance to itself) with NA
  diag(distance_matrix) <- NA
  distance_matrix <- as.data.frame(distance_matrix)
  
  # Set column and row names of the distance matrix
  colnames(distance_matrix) <- AllCells$`Cell Type`
  rownames(distance_matrix) <- AllCells$`Cell Type`
  
  # Subset the distance matrix to include only macrophages
  celltx <- paste0("Macrophages", "_")
  distanceSub <- as.data.frame(distance_matrix[grep(celltx, rownames(distance_matrix)), ])
  
  # Create an empty dataframe to store the closest distances for the current sample
  if (nrow(distanceSub) > 0) {
    closest_distances <- data.frame(matrix(NA, ncol = nrow(celltypes), nrow = nrow(distanceSub)))
    colnames(closest_distances) <- col_names
    rownames(closest_distances) <- rownames(distanceSub)
    
    # Calculate the minimum distance for each cell type
    for (j in celltypes$`unique(SingleCells$cellPhenotype)`) {
      cellt <- gsub("\\^", "\\\\\\^", j)      # Escape '^'
      cellt <- gsub("\\{", "\\\\\\{", cellt)  # Escape '{'
      cellt <- gsub("\\}", "\\\\\\}", cellt)  # Escape '}'
      cellt <- gsub("\\+", "\\\\\\+", cellt)  # Escape '+'
      cellt <- paste0(cellt, "_")
      distanceSub2 <- as.data.frame(distanceSub[, grep(cellt, colnames(distanceSub))])
      closest_distances[, j] <- apply(distanceSub2, 1, function(row) min(row, na.rm = TRUE))
    }
    # Append the closest distances for the current sample to the overall dataframe
    Min_Distance_All <- rbind(Min_Distance_All, closest_distances)
  }
}

# Write the closest distances to a CSV file
write.csv(Min_Distance_All, file = "~/LumA/Min_Distance_All.csv",
          sep = "\t", row.names = TRUE, quote = FALSE, fileEncoding = "UTF-8")
