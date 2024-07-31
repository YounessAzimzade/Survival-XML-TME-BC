# DEA and PEA in HLA-ABC^{lo} vs. HLA-ABC^{hi} in scRNA seq
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 31.07.2024

 
# Load necessary libraries for the analysis
library(plyr)
library(dplyr) 
library(clusterProfiler)
library(fgsea)
library(ggplot2)
library(ggpubr)  
library(signatureSearch)
library(ggrepel) 
library(reshape2)   
library(Seurat)

# Read the Seurat object from the specified file
seurat_obj <- readRDS(file = "D:/Dropbox/scRNA/Seurat_All.rds")

# Subset the Seurat object to include only Macrophage cells
seurat_obj <- subset(seurat_obj, Cell_Type_All == "Macrophage")

# Read the metadata for Macrophages from the specified file
Macs <- read.delim("D:/Dropbox/scRNA/MetadataMacs.txt")
Macs$Cell_ID <- gsub("-",".", Macs$Cell_ID) # Replace '-' with '.' in Cell_IDs

# Add HLA information to the Seurat object's metadata
seurat_obj@meta.data$HLA <- Macs$HLA.ABC.Level[match(seurat_obj@meta.data$Cell_ID, Macs$Cell_ID)]

# Extract metadata from the Seurat object
metadata <- seurat_obj@meta.data

#############################################################################
##############################################################################

## SEURAT

# Normalize the data in the Seurat object
seurat_obj <- NormalizeData(seurat_obj)

# Identify the 2000 most highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# Perform PCA (Principal Component Analysis)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Find neighbors based on the PCA results
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Find clusters in the data
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Identify markers between two conditions, HLA-ABC low and high
markers <- FindMarkers(seurat_obj, ident.1 = "HLA-ABC^{lo}", ident.2 = "HLA-ABC^{hi}", group.by = 'HLA')

# Print the top markers
head(markers)

# Add a column to indicate significant markers based on adjusted p-value and log fold change
markers$significant <- with(markers, ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 2, "Significant", "Not Significant"))

# Create a volcano plot of the markers
padj_threshold <- 0.05
logFC_threshold <- 2
markers$GeneName <- rownames(markers)

pdf("D:/Dropbox/scRNA/Fig5g.pdf", width = 14, height =10)
ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(size=3, aes(color = ifelse(abs(avg_log2FC) > logFC_threshold & p_val_adj < padj_threshold, "violet", "grey")), alpha = 0.6) +
  geom_text_repel(
    aes(label = ifelse(abs(avg_log2FC) > logFC_threshold & p_val_adj < padj_threshold, GeneName, "")),
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.2,
    segment.color = "grey50",
    size=7
  ) +
  scale_color_identity() +
  theme_bw() + 
  labs(
    x = "log2 Fold Change",
    y = "-log10(padj)",
    color = "Significant" 
  ) +  ylim(0,80)+ 
  theme(
    text = element_text(size = 28),
    axis.title.x = element_text(size = 28), # Set x-axis label size
    axis.title.y = element_text(size = 28)  # Set y-axis label size
  ) +
  annotate("text", x = -logFC_threshold - 0.95, y = 75, label = "HLA-ABC^{lo}", size = 9, color = "blue") +
  annotate("text", x = logFC_threshold + 0.5, y = 75, label = "HLA-ABC^{hi}", size = 9, color = "red")
dev.off()

#######################################################################################################################
#######################################################################################################################
#### PEA

# Load the required libraries for pathway enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db) 
library(msigdbr) 

# Extract gene names from markers
markers$gene <- rownames(markers)

# Select upregulated genes for HLA-ABC low
upregulated_ER_n <- markers %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  pull(gene)

# Select upregulated genes for HLA-ABC high
upregulated_ER_p <- markers %>% 
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>% 
  pull(gene)

# Function to run enrichment analysis
run_enrichment <- function(genes) {
  entrez_ids <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  enrichment_results <- enricher(entrez_ids, TERM2GENE = m_t2g)
  return(enrichment_results)
}

# Run enrichment for upregulated genes in HLA-ABC low
enrichment_results_ER_n <- run_enrichment(upregulated_ER_n)

# Run enrichment for upregulated genes in HLA-ABC high
enrichment_results_ER_p <- run_enrichment(upregulated_ER_p)

# Function to plot the enrichment results
plot_enrichment <- function(enrichment, condition) {
  barplot(enrichment, showCategory = 10, title = paste("Hallmark Pathway Enrichment -", condition))
  dotplot(enrichment, showCategory = 10, title = paste("Hallmark Pathway Enrichment -", condition))
}

# Convert enrichment results to data frames
df_ER_n <- as.data.frame(enrichment_results_ER_n)
df_ER_p <- as.data.frame(enrichment_results_ER_p)

# Add a condition column
df_ER_n$Condition <- "HLA-ABC^{lo}"
df_ER_p$Condition <- "HLA-ABC^{hi}"

# Combine the data frames
combined_results <- rbind(df_ER_n, df_ER_p)

# Select the top enriched pathways
top_terms <- combined_results %>%
  group_by(Condition) %>%
  top_n(n = 12, wt = -p.adjust) %>%
  pull(Description) %>%
  unique()

# Filter for the top enriched pathways
combined_results <- combined_results %>%
  filter(Description %in% top_terms)

# Create a bar plot to compare pathway enrichment
pdf("D:/Dropbox/scRNA/Fig5h.pdf", width = 11, height = 7)
ggplot(combined_results, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +  # Adjust the width parameter here
  coord_flip() +
  labs(x = "Pathway", y = "-log10(p.adjust)") +
  theme_minimal() +
  scale_fill_manual(values = c("HLA-ABC^{lo}" = "blue", "HLA-ABC^{hi}" = "red")) +
  theme(legend.position = c(0.5, 0.5)) +
  theme(text = element_text(size = 19)) 
dev.off()
