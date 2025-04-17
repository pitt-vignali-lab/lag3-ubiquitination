# ----------------------------------------------------------------------------
# Fig 7L | CD8+ functional cluster distribution at baseline and after LAG3 blockade
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/4.1.LAGCBLCDinTIL.R
# Input    : 1.Data/cd8SetScvi.RDS
# Output   : 3.Output/workingDf_pid_cat.csv
#            3.Output/CD8_population_distribution_by_treatment_LAGCBL_group*.png
# Upstream : 03_cd8_annotation_module_scores.R
# Note     : Group proportions compared by chi-square test
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ggpubr)

obj <- readRDS( '1.Data/cd8SetScvi.RDS')

# Subset the object for TIL samples at baseline
objTilBsln <- subset(obj, sample_type == 'TIL' & 
                     treatment == 'Baseline')

pidRela <- unique(subset(obj, sample_type == 'TIL' & 
                    treatment == 'rela')$patient_id)

# Add LAG3 and CBLs scores
objTilBsln <- AddModuleScore(
  JoinLayers(objTilBsln), 
  features = list(c('LAG3', 'CBL', 'CBLC'),
                  'LAG3',
                  c('CBL', 'CBLC')),
  name = c('LAGCBL', 'LAG3', 'CBLs')
)

# spliting the overall patinet_id into two groups, with hgiher and lower LAGCBL scores with k meas method
# Splitting the overall patient_id into two groups, with higher and lower LAGCBL scores using k-means method
df <- rbind(data.frame(pid = objTilBsln$patient_id,
                 exp = objTilBsln$LAG32,
                 cat = 'LAG3'),
                 data.frame(pid = objTilBsln$patient_id,
                 exp = objTilBsln$CBLs3,
                 cat = 'CBLs')) 

plotDf <- df %>% 
  group_by(pid, cat) %>% 
  summarise(mean_exp = mean(exp), .groups = 'drop') %>%  # Calculate mean expression by pid and category
  pivot_wider(names_from = cat, values_from = mean_exp) %>%  # Reshape data to wide format with categories as columns
  mutate(cat = case_when(LAG3 > 0 & CBLs > 0 ~ 'LAG3_CBLs_High',  # Simplified category names
                         LAG3 > 0 & CBLs < 0 ~ 'LAG3_CBLs_Low',
                         TRUE ~ 'DN')) 
  
ggplot(plotDf, aes(x = LAG3, y = CBLs, col = cat)) +
  geom_point() +  # Scatter plot points
  geom_point(data = plotDf %>% 
               filter(pid %in% pidRela), shape = 1, size = 3, show.legend = F) +
  ggrepel::geom_text_repel(aes(label = pid), hjust = -0.2, vjust = 0.5,
                           show.legend = FALSE,
                           size = 8/2.8, fontface = 2) +  # Add labels for each point
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Vertical line at x = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Horizontal line at y = 0
  labs(
    title = 'LAG3 vs CBLs Scores',
    x = 'Mean LAG3 Score',
    y = 'Mean CBLs Score',
    color = bquote(LAG3^{"+"} ~ " and " ~ CBLs^{High/Low})  # Legend title with superscript
  ) +
  theme_minimal() +  # Use minimal theme for a clean look
  theme(aspect.ratio = 1,
        text = element_text(size = 8, face = "bold", color = 'black'),
        legend.text = element_text(size = 8, face = "bold", color = 'black'),
        legend.title = element_blank()) +  # Set aspect ratio and text size
  scale_color_discrete(labels = c(bquote(LAG3^{"+"} ~ CBLs^{High}), 
                                bquote(LAG3^{"+"} ~ CBLs^{Low}), 
                                "DN")) 

df %>% 
  group_by(pid, cat) %>% 
  summarise(mean_exp = mean(exp), .groups = 'drop') %>%  # Calculate mean expression by pid and category
  pivot_wider(names_from = cat, values_from = mean_exp) %>%  # Reshape data to wide format with categories as columns
  mutate(cat = case_when(LAG3 > 0 & CBLs > 0 ~ 'LAG3_CBLs_High',  # Simplified category names
                         LAG3 > 0 & CBLs < 0 ~ 'LAG3_CBLs_Low',
                         TRUE ~ 'DN')) %>% 
  write.csv('3.Output/workingDf_pid_cat.csv', row.names = F)

# objTilBsln <- AddMetaData(objTilBsln, 
#                     data.frame(clst = objTilBsln$scvi_clusters) %>% 
#                       rownames_to_column('cell') %>% 
#                       merge(data.frame(clst = 0:7,
#                    anno = factor(c('Effect memory', 'Exhausted', 'NK-like', 'MAIT', 'NK-like',
#                                    'Centory memory', 'Cycling', 'Progenitor exhausted'),
#                                  levels = c('Centory memory', 'Effect memory', 'NK-like', 'MAIT', 
#                                             'Progenitor exhausted', 'Exhausted', 'Cycling'))),
#         by = 'clst') %>% 
#                       select(cell, anno) %>% 
#                       column_to_rownames('cell') )

# objTilBsln <- AddMetaData(objTilBsln, 
#                           metadata = data.frame(pid = objTilBsln$patient_id) %>% 
#                            rownames_to_column('cell') %>%
#                             merge(df %>% 
#   group_by(pid, cat) %>% 
#   summarise(mean_exp = mean(exp), .groups = 'drop') %>%  # Calculate mean expression by pid and category
#   pivot_wider(names_from = cat, values_from = mean_exp) %>%  # Reshape data to wide format with categories as columns
#   mutate(cat = case_when(LAG3 > 0 & CBLs > 0 ~ 'LAG3_CBLs_High',  # Simplified category names
#                          LAG3 > 0 & CBLs < 0 ~ 'LAG3_CBLs_Low',
#                          TRUE ~ 'DN')), by = 'pid') %>%
#                          select(cell, cat) %>%
#                           column_to_rownames('cell')) 


# data.frame(Embeddings(objTilBsln, reduction = 'umap.scvi'),
#            cat = objTilBsln$cat,
#            anno = objTilBsln$anno) %>%  
#   mutate(cat = factor(cat, levels = c('DN', 'LAG3_CBLs_Low', 'LAG3_CBLs_High'))) %>%
#   ggplot(aes(x = umapscvi_1, y = umapscvi_2, col = anno)) + 
#   geom_point() + 
#   facet_wrap(~cat) +
#   coord_fixed() +
#   theme_void() +
#   NoLegend()
# ggsave('3.Output/umap_scvi_LAG3_CBLs_til_bsln_cat_split.png', width = 9, height = 3)

# table(objTilBsln$cat, objTilBsln$anno) %>% 
#   data.frame() %>% 
#   group_by(Var1) %>% 
#   mutate(ratio = Freq * 100 / sum(Freq),
#          Var1 = factor(Var1, levels = c('LAG3_CBLs_High', 'LAG3_CBLs_Low', 'DN'))) %>%
#   ggplot(aes(y = Var1, x = ratio, fill = Var2)) +
#   geom_bar(stat = 'identity', col = 'gray20') +
#   labs(x = 'Percentage (%)',
#        y = 'Tissue',
#        fill = 'Cluster') +
#   coord_fixed(10) +
#   theme_bw() +
#   theme(legend.position = 'bottom',
#         text = element_text(face = 2, size = 8, color = 'black')) 
# ggsave('3.Output/Functional_clusters_distribution_in_TIL_samples.png', width = 8, height = 5)

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

plotDf <- df %>% 
  group_by(pid, cat) %>% 
  summarise(mean_exp = mean(exp), .groups = 'drop') %>%  # Calculate mean expression by pid and category
  pivot_wider(names_from = cat, values_from = mean_exp) %>%  # Reshape data to wide format with categories as columns
  mutate(cat = case_when(LAG3 > 0 & CBLs > 0 ~ 'LAG3+CBLs+',  # Classify groups based on CBLs and LAG3 values
                         LAG3 > 0 & CBLs < 0 ~ 'LAG3+CBLs-',
                         TRUE ~ 'DN'))  

ggplot(plotDf, aes(x = LAG3, y = CBLs, col = cat)) +
    geom_point() +  # Scatter plot points
  geom_point(data = plotDf %>% 
               filter(pid %in% pidRela), shape = 1, size = 3, show.legend = F) +
  geom_point(data = plotDf %>% 
               filter(pid %in% pidRela), shape = 3, size = 3, show.legend = F) +
    ggrepel::geom_text_repel(aes(label = pid), hjust = -0.2, vjust = 0.5,
                             show.legend = FALSE,
                             size = 8/2.8, fontface = 2) +  # Add labels for each point
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Vertical line at x = 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Horizontal line at y = 0
    scale_color_discrete(labels = c('LAG3+CBLs+' = bquote(LAG3^{"+"} ~ CBLs^{High}), 
                                  'LAG3+CBLs-' = bquote(LAG3^{"+"} ~ CBLs^{Low}), 
                                  'DN' = "DN")) +
    labs(
      title = 'LAG3 vs CBLs Scores',
      x = 'Mean LAG3 Score',
      y = 'Mean CBLs Score'
    ) +
    theme_minimal() +  # Use minimal theme for a clean look
    theme(legend.position = 'bottom',  # Set legend position
          aspect.ratio = 1,
          text = element_text(size = 8, face = "bold", color = 'black'))  # Set aspect ratio and text size

ggsave('3.Output/LAG3_vs_CBLs.pdf', width = 8, height = 8)

# according to exp columen in different catagories, separte the pateint into several groups
dfSum <- df %>%
  group_by(pid, cat) %>% 
  summarise(mean_exp = mean(exp), .groups = 'drop') %>%  # Calculate mean expression by pid and category
  pivot_wider(names_from = cat, values_from = mean_exp) %>%  # Reshape data to wide format with categories as columns
  mutate(
    group = ifelse(CBLs > 0 & LAG3 > 0, 'LAG3+CBLs+',  # Classify groups based on CBLs and LAG3 values
                   ifelse(CBLs < 0 & LAG3 > 0, 'LAG3+CBLs-', 'DN')))  # Set factor levels for group

# Merge the group information back to the original object
objTilBsln$ExpGroup <- dfSum$group[match(objTilBsln$patient_id, dfSum$pid)]

data.frame(pid = objTilBsln$patient_id, 
           LAGCBL = objTilBsln$LAGCBL1, 
           group = objTilBsln$ExpGroup) %>%
  group_by(pid, group) %>%
  mutate(LAGCBLmedian = median(LAGCBL), .groups = 'drop',
         group = factor(group, levels = c('DN', 'LAG3+CBLs-', 'LAG3+CBLs+'))) %>%
  ggplot(aes(x = reorder(pid, LAGCBLmedian), y = LAGCBL, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    labs(title = 'Distribution of LAG3 and CBLs Score by Group (Patient-wise)',
         x = 'Patient ID',
         y = expression(Score[LAG3~"&"~CBLs])) +
    facet_grid(~group, scales = 'free_x', space = 'free_x') +
    theme_minimal() +
    theme(legend.position = 'none',
          text = element_text(size = 8, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave( '3.Output/LAGCBL_scores_by_group_patient_wised.png', width = 6, height = 6)

# evluate the functional cluster shif in TIL samples under different treat ment alon with LAGCBL scores
annoDf <- read.csv('3.Output/workingDf.csv') %>% 
  select(cell, anno) %>% 
  filter(cell %in% Cells(obj)) %>%
  column_to_rownames('cell') %>% 
  mutate(anno = factor(anno,
                       levels = c('Centory memory', 'Effect memory', 'NK-like', 'MAIT', 
                                  'Progenitor exhausted', 'Exhausted', 'Cycling')))

obj$ExpGroup <- dfSum$group[match(obj$patient_id, dfSum$pid)]
table(obj$ExpGroup)

obj <- AddMetaData(obj, metadata = annoDf)

data.frame(Embeddings(obj, reduction = 'umap.scvi'),
           treat = obj$treatment,
           anno = obj$anno,
           cat = obj$ExpGroup) %>% 
  filter(treat == 'Baseline',
         !is.na(cat)) %>% 
  ggplot(aes(umapscvi_1, umapscvi_2, col = anno)) +
  geom_point() +
  facet_grid(~cat) +
  coord_fixed() +
  theme_void() +
    NoLegend()
  
ggsave('3.Output/umap_scvi_LAG3_CBLs_til_bsln_cat_split.png', width = 9, height = 3)


data.frame(treat = obj$treatment,
           anno = obj$anno,
           cat = obj$ExpGroup) %>% 
  filter(treat == 'Baseline',
         !is.na(cat)) %>%
  group_by(anno, cat) %>% 
  summarise(Freq = n()) %>% 
  # data.frame() %>%
  group_by(cat) %>%
  mutate(ratio = Freq * 100 / sum(Freq),
         cat = factor(cat, levels = c('LAG3+CBLs+', 'LAG3+CBLs-', 'DN'))) %>%
  ggplot(aes(y = cat, x = ratio, fill = anno)) +
  geom_bar(stat = 'identity', col = 'gray20') +
  labs(x = 'Percentage (%)',
       y = 'Tissue',
       fill = 'Cluster') +
  coord_fixed(10) +
  theme_bw() +
  theme(legend.position = 'bottom',
        text = element_text(face = 2, size = 8, color = 'black'))

ggsave('3.Output/Functional_clusters_distribution_in_TIL_samples.png', width = 8, height = 5)

# table(obj$ExpGroup, obj$anno, obj$treatment) %>% 
#   as.data.frame()  %>% 
#   group_by(Var1, Var3) %>% 
#   mutate(ratio = Freq * 100 / sum(Freq),
#          Var1 = factor(Var1, levels = c('DN', 'LAG3+CBLs-', 'LAG3+CBLs+'))) %>% 
#   filter(Freq != 0,
#          Var3 != 'combo') %>%  
#   ggplot(aes(x = Var1, y = ratio, fill = Var2)) +
#   geom_bar(stat = 'identity') +
#   # geom_point() +
#   # geom_line(aes(group = Var1)) +
#   labs(title = 'Functional Clusters in TIL Samples',
#        x = 'LAG3 and CBLs Groups',
#        y = 'Proportion (%)') +
#   facet_grid(Var3~Var2) +
#   theme_bw() +
#   theme(aspect.ratio = 1,
#         legend.position = 'none',
#         strip.background = element_blank(),
#         text = element_text(size = 8, face = "bold"))
# 
# ggsave('3.Output/Functional_clusters_in_TIL_samples.png', width = 14, height = 7)

# Assuming `obj$ExpGroup`, `obj$anno`, and `obj$treatment` are inputs to create a table
table(obj$ExpGroup, obj$anno, obj$treatment) %>%
  as.data.frame() %>%
  group_by(Var1, Var3) %>%
  mutate(
    ratio = Freq * 100 / sum(Freq),  # Calculate ratio as percentage
    Var1 = factor(Var1, levels = c('DN', 'LAG3+CBLs-', 'LAG3+CBLs+'))  # Order levels
  ) %>%
  filter(
    Freq != 0,        # Exclude rows where Freq is 0
    Var3 != 'combo'   # Exclude rows where Var3 is 'combo'
  ) %>%
  ggplot(aes(x = Var1, y = ratio, fill = Var2, alpha = Var3)) +
  geom_bar(stat = 'identity', position = position_dodge()) +  # Position bars side by side
  scale_alpha_manual(values = c(0.5, 1)) +  # Set alpha values for fill and pattern
  labs(
    title = 'Functional Clusters in TIL Samples',
    x = 'LAG3 and CBLs Groups',
    y = 'Proportion (%)',
    alpha = 'Treatment'
  ) +
  guides(fill = "none") +
  facet_grid( ~ Var2) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = 'bottom',
    strip.background = element_blank(),
    text = element_text(size = 8, face = "bold")
  )

ggsave('3.Output/Functional_clusters_in_TIL_samples.png', width = 14, height = 4)

# evluate the distribution of CD8 popultation distribuion in TIL samples under different treat ment alon with LAGCBL scores
objAll <- readRDS('1.Data/18071_overall_refined_cell_types_ser_obj_230702.rds')

objSub <- subset(objAll, 
                 sample_type == 'TIL' &
                 timepoint %in% c('Baseline', 'Week4'))

objSub$treatment <- ifelse(objSub$timepoint == 'Baseline',
                           'Baseline',
                           objSub$treatment)

objSub <- subset(objSub, treatment == 'nivo', invert = T)

objSub$ExpGroup <- dfSum$group[match(objSub$patient_id, dfSum$pid)]

# Plotting the distribution of CD8 population in TIL samples under different treatments along with LAGCBL scores
# add stat_summary to show the mean
library(ggpubr)

table(treat = objSub$treatment, pid = objSub$patient_id) %>% 
  data.frame() %>% 
  filter(treat == 'rela', Freq != 0) %>% 
  pull(pid) -> relapids
  
table(treat = objSub$treatment, pid = objSub$patient_id) %>% 
  data.frame() %>% 
  filter(treat == 'combo', Freq != 0) %>% 
  pull(pid) -> combopids

table(objSub$treatment, objSub$ExpGroup, objSub$scvi_cell_types_refined, objSub$patient_id) %>% 
  as.data.frame() %>% 
   filter(Var1 %in% c('Baseline', 'combo'), # Corrected filtering syntax
         Var4 %in% combopids) %>%
  group_by(Var1, Var2) %>% 
  mutate(ratio = Freq * 100 / sum(Freq)) %>% 
  filter(Var3 == 'CD8+ T cells') %>% 
  ggplot(aes(x = Var1, y = ratio)) +
  # geom_boxplot() +
  geom_violin() +
  geom_line(aes(group = Var4)) +
  geom_point() +
  # stat_compare_means() +  # Add this line to show p-values
  facet_grid(~Var2) +
  labs(x = '',
       y = 'Percentage (%)') +
  theme_minimal() +
  theme(legend.position = 'bottom',
        text = element_text(size = 8, face = "bold"))

ggsave( '3.Output/CD8_population_distribution_by_treatment_LAGCBL_group in_combo.png', width = 4, height = 4)

table(objSub$treatment, objSub$ExpGroup, objSub$scvi_cell_types_refined, objSub$patient_id) %>% 
  as.data.frame() %>% 
  filter(Var1 %in% c('Baseline', 'rela'), # Corrected filtering syntax
         Var4 %in% relapids) %>%
  group_by(Var1, Var2) %>%
  mutate(ratio = Freq * 100/sum(Freq)) %>%
  filter(Var3 == 'CD8+ T cells') %>%
  ggplot(aes(x = Var1, y = ratio)) +
  # geom_boxplot() +
  geom_violin() + 
  geom_line(aes(group = Var4)) +
  geom_point() +
  # stat_compare_means() + # Add this line to show p-values
  facet_grid(~Var2) +
  labs(x = ' ',
       y = 'Percentage (%)') +
  theme_minimal() +
  theme(legend.position = 'bottom',
        text = element_text(size = 8, face = "bold"))

ggsave( '3.Output/CD8_population_distribution_by_treatment_LAGCBL_group in_rela.png', width = 4, height = 4)
