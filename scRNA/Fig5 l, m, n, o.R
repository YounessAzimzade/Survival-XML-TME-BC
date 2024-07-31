# Cell cell communcation assay 
# Based on  https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md
# Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 05.07.2024
# Load necessary libraries
library(nichenetr)
library(tidyverse)
library(circlize)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(Matrix)
library(readr)
library(stringi)
library(ggpubr)  
library(bigmemory)
library(Seurat)

# Load the list of genes of interest
Genes_OI <- read.delim("~/T cells genes.txt")   
Genes_OI <- unique(Genes_OI)

# Define the path to the data
path <- "D:/DropBox/scRNA/" 

# Load scRNA counts for macrophages and T cells
Counts_Macs <- read.delim(paste0(path, "scMacsAll.txt"))
Counts_T_cells <- read.delim2(paste0(path, "scTcellsAll.txt"))

# Load metadata for macrophages and T cells
MacsID <- read.delim(paste0(path, "MetadataMacs.txt"))
T_cellsID <- read.delim(paste0(path, "MetadataTcells.txt"))

# Filter macrophages with low HLA-ABC level and T cells of CD8+ type
Macs_ids <- MacsID %>% filter(`HLA.ABC.Level` == "HLA-ABC^{lo}") %>% pull(Cell_ID)
T_cells_ids <- T_cellsID %>% filter(`celltype_minor` == "T cells CD8+") %>% pull(Cell_ID)

# Process macrophages counts to find expressed genes with a threshold
EG_Mac_hi <- Counts_Macs[,colnames(Counts_Macs) %in% Macs_ids] %>% 
  apply(1, function(x) {10 * (2 ** x - 1)}) %>% 
  apply(2, function(x) {log2(mean(x) + 1)}) %>% 
  .[. >= 4] %>% 
  names()

# Process T cells counts to find expressed genes with a threshold
EG_CD4T <- Counts_T_cells[,colnames(Counts_T_cells) %in% T_cells_ids] %>% 
  apply(1, function(x) {10 * (2 ** x - 1)}) %>% 
  apply(2, function(x) {log2(mean(x) + 1)}) %>% 
  .[. >= 4] %>% 
  names()

# Load the ligand-target matrix and ligand-receptor network
ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
lr_network <- readRDS("lr_network_human_21122021.rds")

# Identify expressed genes in sender and receiver cells
expressed_genes_sender <- EG_Mac_hi
expressed_genes_receiver <- EG_CD4T

# Define genes of interest and background expressed genes
geneset_oi <- Genes_OI$Gene.Symbol
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# Identify ligands and receptors in the network
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands, expressed_genes_sender)

receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

# Filter the ligand-receptor network for expressed ligands and receptors
lr_network_expressed <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

# Identify potential ligands
potential_ligands <- lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

# Predict ligand activities
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi, 
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix, 
  potential_ligands = potential_ligands
)

# Arrange ligand activities and select the top ligands
ligand_activities %>% arrange(-aupr_corrected) 
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
head(best_upstream_ligands)

# Show histogram of ligand activity scores
p_hist_lig_activity <- ggplot(ligand_activities, aes(x = aupr_corrected)) + 
  geom_histogram(color = "black", fill = "darkorange") + 
  geom_vline(aes(xintercept = min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))), color = "red", linetype = "dashed", size = 1) + 
  labs(x = "ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

# Get weighted ligand-target links for the top ligands
active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, 
                                                                  geneset = geneset_oi,
                                                                  ligand_target_matrix = ligand_target_matrix, 
                                                                  n = 250) %>% 
  bind_rows()

# Remove NA values from the data frame
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df, 
  ligand_target_matrix = ligand_target_matrix, 
  cutoff = 0.25
)

# Order ligands and targets for visualization
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique()
vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()

# Load the scales library
library(scales)

# Create heatmap of ligand-target interactions
p_ligand_target_network <- vis_ligand_target %>%
  make_heatmap_ggplot("Prioritized ligands in HLA-ABC{lo} macrophages", "genes in CD8 T cells", color = "#00BFC4",
                      legend_position = "top", x_axis_position = "top", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke", high = "#00BFC4",
                       breaks = c(0, 0.05, 0.1),
                       limits = c(0, 0.1),  # Set the limits for the color scale
                       oob = squish,  # Use squish to cap values outside the limits
                       guide = guide_colorbar(order = 1, title.position = "top", title.hjust = 0.5)) +
  theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
        legend.text = element_text(size = 8))

# Save the plot as a PDF file
pdf(paste0(path, "/Fig5l.pdf"), width = 6, height = 5.5)
p_ligand_target_network
dev.off() 
