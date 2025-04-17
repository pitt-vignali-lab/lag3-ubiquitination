# ----------------------------------------------------------------------------
# Fig S8 | GESECA functional enrichment before and after relatlimab, by patient group
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/4.2.funcEnrich.R
# Input    : 1.Data/cd8SetScvi.RDS
#            3.Output/workingDf_pid_cat.csv
# Output   : 3.Output/sigPtwCD8.csv
#            3.Output/funcEnrich_new.pdf
# Upstream : 08_cluster_distribution_by_treatment.R
#            utils_plot_helpers.R
# Note     : KEGG gene sets via msigdbr
#            heatmap colours are Z-scaled enrichment scores
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(limma)
library(ComplexHeatmap)

obj <- readRDS('1.Data/cd8SetScvi.RDS')

objTil <- subset(obj, sample_type == 'TIL' & 
                     treatment %in% c('Baseline', 'rela'))

pidCat <- read.csv('3.Output/workingDf_pid_cat.csv') 

# Add metadata and filter out cells with missing categories
objTil <- objTil %>%
  AddMetaData(metadata = pidCat %>% 
    select(pid, cat)  %>% 
    merge(data.frame(pid = objTil$patient_id) %>% rownames_to_column('cell'), by = 'pid') %>% 
    select(cell, cat) %>% 
    column_to_rownames('cell'))

objTil <- subset(objTil, patient_id %in% pidCat$pid)

objTil$condi <- paste(objTil$treatment, objTil$cat, sep = '_')

table(objTil$condi)

es <- AggregateExpression(objTil, group.by = 'condi', assays = 'RNA')[[1]]
es <- es[rowSums(es) != 0, ]
es2 <- normalizeBetweenArrays(log2(as.matrix(es)), method="quantile")
es2[is.infinite(es2)] <- 0
es2[is.na(es2)] <- 0

pathwaysDF <- msigdbr("human", category="C7", subcategory = 'IMMUNESIGDB') 

pathways <- split(as.character(pathwaysDF$gene_symbol), pathwaysDF$gs_name)

gesecaRes <- geseca(pathways, es2, minSize = 15, maxSize = 500)

gesecaResCD8 <- gesecaRes %>% 
  filter(grepl('CD8', pathway),
         padj < 0.05)

plotGesecaTable(gesecaResCD8, pathways, E=es2)

write.csv(gesecaResCD8, '3.Output/sigPtwCD8.csv')

source('2.Script/plotInfo.r')

selectedPtw <- read.csv('1.Data/Relevant_CD8_Pathways_with_Detailed_Categories.csv') %>% 
    select(pathway, detailed_category) %>% 
    filter(detailed_category != 'Other')


t <- generate_plot_data(gesecaResCD8 %>% filter(pathway %in% selectedPtw$pathway), 
                        pathways, 
                        E=es2) 

ptw <- t$plot_data %>% 
  filter(abs(value) > 1.7) %>% 
  pull(pathway) %>% 
  unique()

mat <- t$plot_data %>% 
  pivot_wider(names_from = sample, values_from = value) %>% 
  # filter(pathway %in% ptw) %>% 
  column_to_rownames('pathway') %>% 
  as.matrix() 

pdf('3.Output/funcEnrich_new.pdf', width = 8, height = 9)
draw(Heatmap(mat[selectedPtw$pathway, ],
             row_names_max_width = max_text_width(
               rownames(mat), 
               gp = gpar(fontsize = 8)
             ),
             row_split = selectedPtw$detailed_category,
             row_title_rot = 0,
             row_title_gp = gpar(fontsize = 10, fontface = "bold"),
             column_title_gp = gpar(fontsize = 10, fontface = "bold"),  # Optional: Bold for column titles
             column_names_rot = 45,  # Fixes to column rotation parameter
             row_names_gp = gpar(fontsize = 8, fontface = "bold"),
             column_names_gp = gpar(fontsize = 8, fontface = "bold"),
             rect_gp = gpar(col = "white", lwd = 1),
             column_dend_height = unit(15, "mm"),
             row_dend_width = unit(15, "mm"), 
             width = ncol(mat) * unit(5, "mm"),
             height = nrow(mat) * unit(5, "mm"),
             heatmap_legend_param = list(title = "Normalized Expression Values",
                                         direction = "horizontal",
                                         title_position = "topcenter",
                                         legend_width = ncol(mat) * unit(5, "mm"))
), 
heatmap_legend_side = "bottom")
dev.off()

data_dist <- dist(t(mat), method = "manhattan")
# Hierarchical clustering based on the distance matrix
hc <- hclust(data_dist)

# Plot the dendrogram
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "Samples", sub = "", cex = 0.8)

# Convert the distance object to a matrix for visualization
data_dist_matrix <- as.matrix(data_dist)

# Load the pheatmap library for heatmap visualization
if(!require(pheatmap)) install.packages("pheatmap")
library(pheatmap)

# Plot the heatmap
pheatmap(data_dist_matrix,
         main = "Distance Matrix Heatmap",
         color = colorRampPalette(c("#008cff", "white", "#ea00ff25"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         border_color = NA)
