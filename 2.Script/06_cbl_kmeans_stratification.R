# ----------------------------------------------------------------------------
# Fig 7J | K-means (k=2) stratification of CBL-high vs CBL-low cells
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/cblhilowExp.R
# Input    : 1.Data/cd8SetScvi.RDS
#            3.Output/workingDf.csv
# Output   : 3.Output/cblsExpHiLowUMAP.png
#            3.Output/cblsExpHiLowRatio.png
# Upstream : 03_cd8_annotation_module_scores.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

obj <- readRDS('1.Data/cd8SetScvi.RDS')

annoDf <- read.csv('3.Output/workingDf.csv') %>% 
  select(cell, anno, group, lag3, cbl) %>% 
  column_to_rownames('cell') %>% 
  mutate(anno = factor(anno,
                       levels = c('Centory memory', 'Effect memory', 'NK-like', 'MAIT', 
                                  'Progenitor exhausted', 'Exhausted', 'Cycling')),
         group = factor(group, levels = c('Low', 'High')))


obj <- AddMetaData(obj, metadata = annoDf)

DimPlot(obj, split.by = 'group', group.by = 'anno', reduction = 'umap.scvi') +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom',
        text = element_text(face = 2, size = 8, color = 'black'))

ggsave('3.Output/cblsExpHiLowUMAP.png', width = 8, height = 4.5)

annoDf %>% 
  group_by(anno, group) %>% 
  summarise(freq = n()) %>% 
  group_by(anno) %>% 
  mutate(ratio = freq/sum(freq)) %>% 
  ggplot(aes(x = anno, y = ratio, fill = anno, alpha = group)) +
  geom_bar(stat = 'identity') +
  scale_alpha_manual(breaks = c('High', 'Low'),
                     values = c(0.8, 0.2)) +
  labs(x = '',
       y = 'Distribution Ratio',
       alpha = 'CBLs expression',
       fill = 'Cluster') +
  guides(fill = guide_legend(order = 1),  # Set order for fill legend
         alpha = guide_legend(order = 2)) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.box = 'vertical',
        legend.title.align = 0.5,
        legend.title.position = "top",
        aspect.ratio = 1, 
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        text = element_text(face = 2, size = 8, color = 'black'))
  
ggsave('3.Output/cblsExpHiLowRatio.png', width = 5, height = 7)
