# ----------------------------------------------------------------------------
# Preprocessing | scVI integration and Louvain clustering of CD8+ T cells
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/1.2.cd8Init.R
# Input    : 1.Data/cd8Set.RDS
# Output   : 1.Data/cd8SetScvi.RDS
# Upstream : 01_data_prep.R
# Note     : Requires a conda env with scvi-tools
#            set SCVI_ENV below
# ----------------------------------------------------------------------------

library(Seurat)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(tidyverse)

options(future.globals.maxSize = 1e9)

# t <- readRDS('1.Data/cillo_et_al_figure_2_data.rds')

# Conda environment that provides scvi-tools -- edit to match your setup.
SCVI_ENV <- 'scvi-env'

obj <- readRDS('1.Data/cd8Set.RDS')

# obj <- subset(obj, cells = t$cd8_scrnaseq$cb)

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$patient_id)

obj

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
# 
# obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
# obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
# 
# obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
# DimPlot(obj, reduction = "umap.unintegrated", group.by = c("patient_id"))

obj <- IntegrateLayers(
  object = obj, 
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = SCVI_ENV,
  verbose = TRUE
)

obj <- FindNeighbors(obj, reduction = "scviX", dims = 1:10)
obj <- FindClusters(obj, resolution = 0.3, cluster.name = "scvi_clusters")

obj <- RunUMAP(obj, 
               reduction = "scviX", 
               dims = 1:10, 
               reduction.name = "umap.scvi")
DimPlot(
  obj,
  reduction = "umap.scvi",
  group.by = c("patient_id", "scvi_clusters"),
  combine = T, label.size = 2
) +
  coord_fixed()

FeaturePlot(obj, reduction = "umap.scvi", features = 'LAG3') +
  coord_fixed()

saveRDS(obj, '1.Data/cd8SetScvi.RDS')

# marker gene 
obj2 <- JoinLayers(obj)

DimPlot(obj2, group.by = c("scvi_clusters", 'final_clusters'), reduction = "umap.scvi",
        label = T) +
  coord_fixed()

Idents(obj2) <- "scvi_clusters"
mkG <- COSG::cosg(obj2)

write.csv(mkG, '3.Output/cd8MkG.csv')

obj2 <- AddMetaData(obj2, 
                    t$cd8_scrnaseq %>% 
                      filter(cb %in% Cells(obj2)) %>% 
                      select(final_clusters))


t$cd8_scrnaseq %>%
  ggplot(aes(UMAP_1, UMAP_2, col = final_clusters)) +
  geom_point() +
  coord_fixed()

# cluster 0: effector memory T cells. The combination of T cell-specific markers (GZMK, LTB, TCF7), memory-associated genes (CD44, IL7R, GPR183)
# 1: exhausted CD8+ T cells with cells that are not completely dysfunctional (IFNG, CCL3, CCL4, PRF1), The cells show a clear exhaustion signature while retaining some effector capabilities, which is typical of exhausted T cells that are not terminally differentiated. 
# 2: NK-like CD8+ T cells
# 3: innate-like T cells, specifically MAIT (Mucosal Associated Invariant T) cells.
# 4: NK-like CD8+ T cells
# 5: naïve or memory CD8+ 
# 6: Cell Cycle and Proliferation
# 7: progenitor exhausted T 

FeaturePlot(obj2, features = 'LAG3', reduction = 'umap.scvi') +
  coord_fixed()
