## transfer annotation 

library(Seurat)
library(Matrix)
library(readr)
library(stringi)
library(ggpubr) 
library(dplyr) 
library(bigmemory)


## set dimension for feature transfer function 
dimss <- 60 
 
NatGen_MetaData<-read.csv("D:/Dropbox/scRNA/NatGen21/metadata.csv")
NatGen_Counts<-readMM("D:/Dropbox/scRNA/NatGen21/matrix.mtx") 
NatGen_Genes<-read.table("D:/Dropbox/scRNA/NatGen21/genes.tsv")
rownames(NatGen_Counts) <- NatGen_Genes$V1
colnames(NatGen_Counts) <- NatGen_MetaData$X

# Create Seurat object
NatGen_SO <- CreateSeuratObject(counts=NatGen_Counts)

NatGen_SO@meta.data$Cell_ID <- NatGen_MetaData$X
NatGen_SO@meta.data$BC_Sub <- NatGen_MetaData$subtype
NatGen_SO@meta.data$celltype_major <- NatGen_MetaData$celltype_major
NatGen_SO@meta.data$celltype_minor <- NatGen_MetaData$celltype_minor
NatGen_SO@meta.data$celltype_subset <- NatGen_MetaData$celltype_subset 
 
NatGen_SO <- subset(NatGen_SO, celltype_major=="T-cells")


NatGen_SO <- NormalizeData(NatGen_SO)
NatGen_SO <- FindVariableFeatures(NatGen_SO)
NatGen_SO <- ScaleData(NatGen_SO)


remove(NatGen_Counts,NatGen_Genes,NatGen_MetaData)


path <- "D:/DropBox/scRNA/NatMed20/Data/"
NatMed_Counts<-  readRDS(paste0(path, "counts_cells_cohort1.rds"))
NatMed_SO <- CreateSeuratObject(counts = NatMed_Counts)

raw_text<-readLines(paste0(path, "metaData_cohort1_web.csv"),encoding="UTF-8")
correct_text <- stri_encode(raw_text, from = "UTF-8", to = "UTF-8")
combined_text <- paste(correct_text, collapse = "\n")
NatMed_MetaData  <- read_csv(combined_text)


rownames(NatMed_MetaData) <- NatMed_MetaData$Cell

#NatMed_MetaData <- subset(NatMed_MetaData, timepoint=="Pre" )
#table(NatMed_MetaData$cellType)

NatMed_SO@meta.data$orig.ident  <-  NatMed_MetaData$patient_id


NatMed_SO@meta.data$Cell_ID <- NatMed_MetaData$Cell
NatMed_SO@meta.data$BC_Sub <- NatMed_MetaData$BC_type

 
NatMed_SO@meta.data$cellType <-  NatMed_MetaData$cellType
NatMed_SO@meta.data$Phase <-  NatMed_MetaData$timepoint
 
 
NatMed_SO <- subset(NatMed_SO,cellType=="T_cell" & Phase=="Pre" )
 
remove(NatMed_MetaData, NatMed_Counts, raw_text, correct_text)


gc() 
# For the new dataset
NatMed_SO <- NormalizeData(NatMed_SO)
NatMed_SO <- FindVariableFeatures(NatMed_SO)
NatMed_SO <- ScaleData(NatMed_SO)


gc() 
 
# Get gene lists for the reference Seurat object
reference_genes <- rownames(GetAssayData(NatGen_SO, assay = "RNA", slot = "counts"))
# Get gene lists for the query Seurat object
query_genes <- rownames(GetAssayData(NatMed_SO, assay = "RNA", slot = "counts"))
# Find common genes
common_genes <- intersect(reference_genes, query_genes)

gc()
anchors <- FindTransferAnchors(reference=NatGen_SO, query=NatMed_SO, dims=1:dimss, features=common_genes)

# Transfer cell type annotations
predictions <- TransferData(anchorset = anchors, refdata = NatGen_SO$celltype_minor, dims = 1:dimss)

# Add predictions to the new Seurat object
NatMed_SO <- AddMetaData(NatMed_SO, metadata = predictions$predicted.id, col.name = "celltype_minor")

predictions2 <- TransferData(anchorset = anchors, refdata = NatGen_SO$celltype_subset, dims = 1:dimss)
NatMed_SO <- AddMetaData(NatMed_SO, metadata = predictions2$predicted.id, col.name = "celltype_subset")

 
# Merging Seurat objects

NatMed_SO <- subset(NatMed_SO, features = common_genes)
NatGen_SO <- subset(NatGen_SO, features = common_genes)
 
seurat_combined <- merge(NatGen_SO, y = NatMed_SO)
 
 
seurat_combined <- JoinLayers(seurat_combined)

metadata <- seurat_combined@meta.data
 
 
write.table(metadata,file="D:/Dropbox/scRNA/Metadata_T_cells.txt"
            ,sep="\t",row.names = TRUE, col.names = TRUE, quote = FALSE)

Counts   <- GetAssayData(object = seurat_combined, assay = "RNA",slot = "counts")


write.table(Counts,file="D:/Dropbox/scRNA/Counts_T_cells.txt"
            ,sep="\t",row.names = TRUE, col.names = TRUE, quote = FALSE)


 