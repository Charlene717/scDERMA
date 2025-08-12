## Ref: https://chatgpt.com/share/683f81af-51f8-8002-8328-e3a39fccd626
## Ref: https://chatgpt.com/c/683e93ce-9dd0-8002-9948-e6a467075a76

##### Presetting ######
rm(list = ls()) # Clean all variables ##* Comment out if running the entire script

## Speed up (Old Version)
memory.limit(150000)


####
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('tidyverse')) {install.packages('tidyverse'); library(tidyverse)}


if(!require('scMRMA')) {devtools::install_github("JiaLiVUMC/scMRMA"); library(scMRMA)}

####
Set_scMRMA_P <- 0.05


####
## Keloid
seuratObject_Sample <- readRDS(file = "C:/Charlene/Dataset_KGD_Lab/#_Seurat_Object/Keloid_Jojie/TNtype.combined-37.rds")
seuratObject_Sample$Type <- seuratObject_Sample$orig.ident1
seuratObject_Sample$Cell_Type_KGD <- Idents(seuratObject_Sample)

## Extract Fibroblasts
seuratObject_Sub_Fib <- seuratObject_Sample[, grepl("fibroblasts", seuratObject_Sample$Cell_Type_KGD)]




################################################################################
#### scMRMA for Multiple levels of  marker_df_scMRMA_KGD_ChatGPT ####

Set_Run_DotPlots <- FALSE
source("KGD_CTAnnot_MarkerList_Charlene.R")


DefaultAssay(seuratObject_Sample) <- "RNA"
seuratObject_Sample$seurat_clusters <- as.factor(seuratObject_Sample$seurat_clusters)

# Call scMRMA
scMRMA_Step1 <- scMRMA(
  input = seuratObject_Sample,
  p = Set_scMRMA_P, # 0.05,
  normalizedData = F,
  selfDB = marker_df_scMRMA_KGD_ChatGPT,
  selfClusters = seuratObject_Sample$seurat_clusters,
  # selfClusters = NULL,
  k = Set_SubClust_k # 20
)


seuratObject_Sample@meta.data$label_scMRMA_Level1 <- scMRMA_Step1[["multiR"]][["annotationResult"]]$Level1
seuratObject_Sample@meta.data$label_scMRMA_Level2 <- scMRMA_Step1[["multiR"]][["annotationResult"]]$Level2
seuratObject_Sample@meta.data$label_scMRMA_Level3 <- scMRMA_Step1[["multiR"]][["annotationResult"]]$Level3
seuratObject_Sample@meta.data$label_scMRMA_Level4 <- scMRMA_Step1[["multiR"]][["annotationResult"]]$Level4
DimPlot(seuratObject_Sample, reduction = "umap", label = T , group.by = "Cell_Type_KGD") +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_Level4", reduction = "umap", label = T)

DimPlot(seuratObject_Sample, group.by = "label_scMRMA_Level1", reduction = "umap", label = T) +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_Level2", reduction = "umap", label = T) +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_Level3", reduction = "umap", label = T) +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_Level4", reduction = "umap", label = T)

DimPlot(seuratObject_Sample, group.by = "label_scMRMA_Level4", reduction = "umap", label = T)


################################################################################
#### scMRMA for Single levels of marker_df_scMRMA_KGD_ChatGPT ####
Set_Run_DotPlots <- FALSE
source("KGD_CTAnnot_MarkerList_Charlene.R")

marker_df_scMRMA_KGD_ChatGPT_M <- marker_df_scMRMA_KGD_ChatGPT[,c(1,5)]
colnames(marker_df_scMRMA_KGD_ChatGPT_M)[2] <- "LevelS"


DefaultAssay(seuratObject_Sample) <- "RNA"
seuratObject_Sample$seurat_clusters <- as.factor(seuratObject_Sample$seurat_clusters)

# Call scMRMA
scMRMA_S <- scMRMA(
  input = seuratObject_Sample,
  p = Set_scMRMA_P, # 0.05,
  normalizedData = F,
  selfDB = marker_df_scMRMA_KGD_ChatGPT_M,
  selfClusters = seuratObject_Sample$seurat_clusters,
  # selfClusters = NULL,
  k = Set_SubClust_k # 20
)

seuratObject_Sample@meta.data$label_scMRMA_LevelS <- scMRMA_S[["multiR"]][["annotationResult"]]

Set_Com_broad_cell_clusters <- TRUE
if (Set_Com_broad_cell_clusters) {
  source("FUN_Annotate_broad_cell_clusters.R")
  seuratObject_Sample <- annotate_broad_cell_clusters(
    seurat_object      = seuratObject_Sample,
    broad_cell_type_col = "label_scMRMA_LevelS",
    seurat_cluster_col  = "seurat_clusters"
  )
}
seuratObject_Sample$label_scMRMA_LevelS_SeuratClusters <- seuratObject_Sample$BroadCellTypeAnnot_SeuratClusters


DimPlot(seuratObject_Sample, reduction = "umap", label = T , group.by = "Cell_Type_KGD") +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_LevelS", reduction = "umap", label = T)

DimPlot(seuratObject_Sample, reduction = "umap",label.size = 3,  label = T , group.by = "Cell_Type_KGD") +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_LevelS_SeuratClusters", reduction = "umap",label.size = 3,  label = T)

DimPlot(seuratObject_Sample, reduction = "umap", label = T , group.by = "Cell_Type_KGD") +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_LevelS", reduction = "umap", label = T) +
  DimPlot(seuratObject_Sample, group.by = "seurat_clusters", reduction = "umap", label = T)


################################################################################
#### DotPlots for marker_sets_KGD_ChatGPT ####

#### Visualization: DotPlots ####
marker_sets <- marker_sets_KGD_ChatGPT 

Set_Idents <- "seurat_clusters"
Set_Marker_DotPlot_ncol <- 4
Set_Marker_DotPlot_width <- 30
Set_Marker_DotPlot_height <- 60
source("RUNPlot_Marker_DotPlots.R")


Set_Idents <- "Cell_Type_KGD"
Set_Marker_DotPlot_ncol <- 4
Set_Marker_DotPlot_width <- 30
Set_Marker_DotPlot_height <- 80
source("RUNPlot_Marker_DotPlots.R")


Set_Idents <- "label_scMRMA_LevelS"
Set_Marker_DotPlot_ncol <- 4
Set_Marker_DotPlot_width <- 30
Set_Marker_DotPlot_height <- 80
source("RUNPlot_Marker_DotPlots.R")


Set_Idents <- "label_scMRMA_LevelS_SeuratClusters"
Set_Marker_DotPlot_ncol <- 4
Set_Marker_DotPlot_width <- 30
Set_Marker_DotPlot_height <- 80
source("RUNPlot_Marker_DotPlots.R")




################################################################################
#### scMRMA for Single levels of marker_df_scMRMA_KGD ####

source("KGD_CTAnnot_MarkerList.R")

colnames(marker_df_scMRMA_KGD)[2] <- "KGD"


DefaultAssay(seuratObject_Sample) <- "RNA"
seuratObject_Sample$seurat_clusters <- as.factor(seuratObject_Sample$seurat_clusters)

# Call scMRMA
scMRMA_S <- scMRMA(
  input = seuratObject_Sample,
  p = Set_scMRMA_P, # 0.05,
  normalizedData = F,
  selfDB = marker_df_scMRMA_KGD,
  selfClusters = seuratObject_Sample$seurat_clusters,
  # selfClusters = NULL,
  k = Set_SubClust_k # 20
)

seuratObject_Sample@meta.data$label_scMRMA_KGD <- scMRMA_S[["multiR"]][["annotationResult"]]

Set_Com_broad_cell_clusters <- TRUE
if (Set_Com_broad_cell_clusters) {
  source("FUN_Annotate_broad_cell_clusters.R")
  seuratObject_Sample <- annotate_broad_cell_clusters(
    seurat_object      = seuratObject_Sample,
    broad_cell_type_col = "label_scMRMA_KGD",
    seurat_cluster_col  = "seurat_clusters"
  )
}
seuratObject_Sample$label_scMRMA_KGD_SeuratClusters <- seuratObject_Sample$BroadCellTypeAnnot_SeuratClusters


DimPlot(seuratObject_Sample, reduction = "umap", label = T , group.by = "Cell_Type_KGD") +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_KGD", reduction = "umap", label = T)

DimPlot(seuratObject_Sample, reduction = "umap",label.size = 3,  label = T , group.by = "Cell_Type_KGD") +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_KGD_SeuratClusters", reduction = "umap",label.size = 3,  label = T)

DimPlot(seuratObject_Sample, reduction = "umap", label = T , group.by = "Cell_Type_KGD") +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_KGD", reduction = "umap", label = T) +
  DimPlot(seuratObject_Sample, group.by = "seurat_clusters", reduction = "umap", label = T)

