# ----------------------------------------------------------------------------
# Data preparation | Subset CD4+ Tconv and CD8+ T cells from the published trial dataset
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/1.1.datPrep.R
# Input    : 18071 overall Seurat object (GSE225063, Cillo et al. Cell 2024)
# Output   : 1.Data/cd8Set.RDS
#            1.Data/cd4Set.RDS
# Note     : Edit INPUT_OBJ below to point at the downloaded dataset
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(COSG)

# Published dataset from Cillo et al. Cell 2024 (GEO: GSE225063) -- see README.
# Download it into 1.Data/ and point INPUT_OBJ at the file.
INPUT_OBJ <- '1.Data/18071_overall_refined_cell_types_ser_obj_230702.rds'

obj <- readRDS(INPUT_OBJ)

objSub <- subset(obj, 
                 timepoint %in% c('Baseline', 'Week4') &
                 scvi_cell_types_refined %in% c("CD4+ Tconv", "CD8+ T cells"))

objSub$treatment <- ifelse(objSub$timepoint == 'Baseline',
                           'Baseline',
                           objSub$treatment)
objSub <- subset(objSub, treatment == 'nivo', invert = T)

DimPlot(objSub, split.by = 'scvi_cell_types_refined') +
  coord_fixed()

saveRDS(subset(objSub, scvi_cell_types_refined == "CD8+ T cells"),
        '1.Data/cd8Set.RDS')

saveRDS(subset(objSub, scvi_cell_types_refined == "CD4+ Tconv"),
        '1.Data/cd4Set.RDS')

# target gene CBL CBLB

objcbl <- AddModuleScore(subset(objSub, scvi_cell_types_refined == "CD8+ T cells"), 
                         features = list(CBLs = c('CBL', 'CBLB')), name = 'CBLs')

FeaturePlot(objcbl, features = 'CBLs1') +
  coord_fixed() +
  scale_color_gradient2(high = 'darkred', low = 'blue', mid = 'white', 
                        midpoint = 0)

