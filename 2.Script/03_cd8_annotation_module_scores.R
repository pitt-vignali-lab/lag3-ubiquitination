# ----------------------------------------------------------------------------
# Fig 7G | CD8+ TIL cluster annotation; LAG3 and CBL module scores on UMAP
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/2.1.cd8deepDive.R
# Input    : 1.Data/cd8SetScvi.RDS
# Output   : 3.Output/annotedUMAP*
#            3.Output/cd8MkG.pdf
#            3.Output/lag3ExpressUMAP.pdf
#            3.Output/cblExpressUMAP.pdf
#            3.Output/workingDf.csv
# Upstream : 02_cd8_integration_clustering.R
# Note     : Cluster markers identified with COSG
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

obj <- readRDS('1.Data/cd8SetScvi.RDS')

obj <- subset(obj, treatment == 'combo', invert = T)

DimPlot(obj, 
        group.by = "scvi_clusters", 
        reduction = "umap.scvi",
        label = T) +
  coord_fixed()

obj2 <- JoinLayers(obj)

obj2 <- AddModuleScore(obj2, 
                       features = list(lag3 = c('LAG3'),
                                       cbl = c('CBL', 'CBLB')),
                       name = c('lag3', 'cbl'))

# Create the plotDf dataframe
plotDf <- data.frame(
  Embeddings(obj2, reduction = "umap.scvi"),
  clst = obj2$scvi_clusters,
  treat = obj2$treatment,
  tissue = obj2$sample_type,
  lag3 = obj2$lag31,
  cbl = obj2$cbl2,
  pid = obj2$patient_id
) %>% 
  rownames_to_column('cell')

# Function to simplify merging process
merge_response <- function(treatment_type, response_col_name) {
  obj2@meta.data %>%
    filter(treatment == treatment_type) %>%
    select(patient_id, clinical_response_wk4) %>%
    mutate(clinical_response_wk4 = ifelse(clinical_response_wk4 == 'NA', NA, clinical_response_wk4)) %>% 
    unique() %>%
    rename(!!response_col_name := clinical_response_wk4)
}

# Merge with 'rela' treatment
plotDf <- plotDf %>% 
  merge(merge_response('rela', 'relaRes'), by.x = 'pid', by.y = 'patient_id', all = TRUE)

# Merge with 'combo' treatment
# plotDf <- plotDf %>% 
#   merge(merge_response('combo', 'comboRes'), by.x = 'pid', by.y = 'patient_id', all = TRUE)

# cluster 0: effector memory T cells
# 1: exhausted CD8+ T cells with cells that are not completely dysfunctional (IFNG, CCL3, CCL4, PRF1), The cells show a clear exhaustion signature while retaining some effector capabilities, which is typical of exhausted T cells that are not terminally differentiated. 
# 2: NK-like CD8+ T cells
# 3: innate-like T cells, specifically MAIT (Mucosal Associated Invariant T) cells.
# 4: NK-like CD8+ T cells
# 5: naïve or memory CD8+ 
# 6: Cell Cycle and Proliferation
# 7: progenitor exhausted T 

# Step 1: Process `mkG`
mkG <- rbind(data.frame(clst = 0, gene = c('GZMK', 'CD44', 'CXCR4', 'ITM2C', 'IL7R')),
             data.frame(clst = 1, gene = c('PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'CXCL13')),
             data.frame(clst = 2, gene = c('GZMB', 'GNLY', 'FCGR3A', 'KLRD1', 'NKG7')),
             data.frame(clst = 3, gene = c('IL7R', 'KLRB1', 'RORA', 'CXCR4', 'KLRG1')),
             data.frame(clst = 4, gene = c('GZMB', 'GNLY', 'FCGR3A', 'KLRD1', 'NKG7')),
             data.frame(clst = 5, gene = c('LEF1', 'CCR7', 'TCF7', 'SELL', 'IL7R')),
             data.frame(clst = 6, gene = c('MKI67', 'CDK1', 'CCNB1', 'TOP2A', 'BIRC5')),
             data.frame(clst = 7, gene = c('TOX', 'IKZF3', 'BCL2L11', 'TFRC', 'CCR5'))) %>% 
  # pivot_longer(cols = everything(), names_to = 'clst', values_to = 'gene') %>% 
  # mutate(clst = gsub('names.', '', clst)) %>% 
  mutate(clst = as.numeric(clst)) %>%  # Ensure 'clst' is numeric for merging
  merge(data.frame(clst = 0:7,
                   anno = factor(c('Effector memory', 'Exhausted', 'NK-like', 'Effector-like\nTissue-Adapted', 'NK-like',
                                   'Central memory', 'Cycling', 'Progenitor exhausted'),
                                 levels = c('Central memory', 'Effector memory', 'NK-like', 'Effector-like\nTissue-Adapted', 
                                            'Progenitor exhausted', 'Exhausted', 'Cycling'))),
        by = 'clst')

# Step 2: Add Metadata to `obj2`
obj3 <- AddMetaData(obj2, 
                    data.frame(clst = obj2$scvi_clusters) %>% 
                      rownames_to_column('cell') %>% 
                      merge(data.frame(clst = 0:7,
                   anno = factor(c('C2: Effector memory', 
                                   'C6: Exhausted', 
                                   'C3: NK-like', 
                                   'C4: Effector-like Tissue-Adapted', 
                                   'C3: NK-like',
                                   'C1: Central memory', 
                                   'C7: Cycling', 
                                   'C5: Progenitor exhausted'),
                                 levels = c('C1: Central memory', 
                                            'C2: Effector memory', 
                                            'C3: NK-like', 
                                            'C4: Effector-like Tissue-Adapted', 
                                            'C5: Progenitor exhausted', 
                                            'C6: Exhausted', 
                                            'C7: Cycling'))),
        by = 'clst') %>% 
                      select(cell, anno) %>% 
                      column_to_rownames('cell') )

# Step 3: Dimensional Plot
DimPlot(obj3, 
        group.by = 'anno', 
        reduction = 'umap.scvi',
        label = TRUE) +
  coord_fixed()

# Step 4: Dot Plot with Additional Customization
dot_data <- DotPlot(obj3, 
                    features = unique(mkG$gene),
                    group.by = 'anno',
                    dot.scale = 2)$data %>% 
  merge(mkG, by.x = 'features.plot', by.y = 'gene') 

ggplot(dot_data, aes(y = features.plot, x = gsub(':.*', '', id), 
                     fill = avg.exp.scaled, size = pct.exp)) +
  geom_point(shape = 21) + 
  scale_fill_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
  facet_grid(anno~., scales = 'free', space = 'free') +
  labs(x = '', y = '', fill = 'Average Expression', size = 'Percent Expression') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.box = "vertical",
        text = element_text(face = 2, size = 8, color = 'black'),
        strip.text.y = element_text(angle = 0, face = 2, size = 8, color = 'black', hjust = 0),
        strip.background = element_blank())


ggsave('3.Output/cd8MkG.pdf', width = 4, height = 10)

plotDf <- plotDf %>%
  merge(data.frame(clst = 0:7,
                   anno = factor(c('C2: Effector memory', 
                                   'C6: Exhausted', 
                                   'C3: NK-like', 
                                   'C4: Effector-like\nTissue-Adapted', 
                                   'C3: NK-like',
                                   'C1: Central memory', 
                                   'C7: Cycling', 
                                   'C5: Progenitor exhausted'),
                                 levels = c('C1: Central memory', 
                                            'C2: Effector memory', 
                                            'C3: NK-like', 
                                            'C4: Effector-like\nTissue-Adapted', 
                                            'C5: Progenitor exhausted', 
                                            'C6: Exhausted', 
                                            'C7: Cycling'))),
        by = 'clst')

plotDf %>% 
  ggplot(aes(umapscvi_1, umapscvi_2, col = anno)) +
  geom_point() +
  # facet_grid(~tissue) +
  geom_text(data = plotDf %>% group_by(anno) %>% 
              summarise(umapscvi_1 = median(umapscvi_1),
                        umapscvi_2 = median(umapscvi_2)),
            aes(label = anno), col = 'black',
            size = 12/2.8, fontface = 2) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'none')

ggsave('3.Output/annotedUMAP_labeled.pdf', width = 6, height = 6.5)

p <- ggplot()

for (i in c('C1: Central memory', 
            'C2: Effector memory', 
            'C3: NK-like', 
            'C4: Effector-like\nTissue-Adapted', 
            'C5: Progenitor exhausted', 
            'C6: Exhausted', 
            'C7: Cycling')) {
  p <- p + 
    geom_point(data = plotDf %>% filter(anno == i), aes(umapscvi_1, umapscvi_2, fill = anno), shape = 21, size = 2, stroke = 1) + 
    geom_point(data = plotDf %>% filter(anno == i), aes(umapscvi_1, umapscvi_2, col = anno), size = 1)
}


p +
  geom_label(data = plotDf %>% group_by(anno) %>% 
              summarise(umapscvi_1 = median(umapscvi_1),
                        umapscvi_2 = median(umapscvi_2)),
            aes(umapscvi_1, umapscvi_2, label = anno), col = 'black',
            size = 12/2.8, fontface = 2) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'none')

ggsave('3.Output/annotedUMAP_labeled2.pdf', width = 6, height = 6.5)

p +
  # geom_label(data = plotDf %>% group_by(anno) %>% 
  #              summarise(umapscvi_1 = median(umapscvi_1),
  #                        umapscvi_2 = median(umapscvi_2)),
  #            aes(umapscvi_1, umapscvi_2, label = anno), col = 'black',
  #            size = 12/2.8, fontface = 2) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'none')

ggsave('3.Output/annotedUMAP_labeled3.pdf', width = 6, height = 6.5)

plotDf %>% 
  ggplot(aes(umapscvi_1, umapscvi_2, col = anno)) +
  geom_point() +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'none')

ggsave('3.Output/annotedUMA.png', width = 6, height = 6.5)

table(plotDf$anno, plotDf$tissue) %>% 
  data.frame() %>% 
  group_by(Var2) %>% 
  mutate(percentage = (Freq/sum(Freq)) *100) %>%
  ggplot(aes(y = Var2, x = percentage, fill = Var1)) +
  geom_bar(stat = 'identity', col = 'gray20') +
  labs(x = 'Percentage (%)',
       y = 'Tissue',
       fill = 'Cluster') +
  coord_fixed(10) +
  theme_bw() +
  theme(legend.position = 'bottom',
        text = element_text(face = 2, size = 8, color = 'black'))
ggsave('3.Output/clusterTissueDistribution.pdf', width = 8, height = 3)

# 1. patient wised lag3 cbls expression correlation test ----
corRes <- plotDf %>%
  filter(tissue == 'TIL') %>% 
  group_by(pid, tissue, treat) %>%
  summarise(
    cor_value = cor.test(lag3, cbl)$estimate,
    p_value = cor.test(lag3, cbl)$p.value
  ) %>% 
  mutate(fdr = p.adjust(p_value)) # %>%
  # merge(plotDf %>% 
  #         select(pid, relaRes, comboRes) %>% 
  #         unique()) 
  # na.omit() %>% 
  # filter(fdr < 0.05) %>% 
corRes %>% 
  # mutate(relaRes = replace_na(relaRes, "NA"),
  #        comboRes = replace_na(comboRes, "NA")) %>% 
  ggplot(aes(x = reorder(pid, cor_value), y = cor_value)) +
  geom_point(aes(col = fdr < 0.05),
             size = 2) +
  facet_grid(~treat, scales = 'free', , space = "free_x") +
  scale_color_manual(breaks = c(TRUE, FALSE),
                     values = c('red', 'gray')) +
  labs(x = 'Patient',
       y = expression(Cor[italic(LAG3 ~ v.s. ~ CBLs)]),
       col = 'Significant (FDR < 0.05)') +
  theme_bw() +
  theme(# aspect.ratio = 1,
        text = element_text(face = 2, size = 12, color = 'black'),
        axis.text.x = element_text(angle = 90, 
                                   face = 2,
                                   hjust = 1, 
                                   vjust = 0.5),
        legend.position = 'bottom')

   ggsave('3.Output/pre_post_treated_patient_wised_cor_lag3_cbls_TIL_rela_only.pdf', width = 10, height = 5.5)


# 2. lag3 and cbls expression umap ----
  cblDf <- plotDf %>% 
    mutate(cbl = round(cbl, 1))

  p <- ggplot()
  
  # Plotting points for each unique cbl value
  for (i in sort(unique(abs(cblDf$cbl)))) {
    p <- p + geom_point(data = cblDf %>% filter(abs(cbl) == i),
                        aes(umapscvi_1, umapscvi_2, col = cbl))
  }
  
  # Adding color scales and theme
  p + 
    scale_fill_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
    scale_color_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = 'bottom',
          text = element_text(size = 12, face = 2),
          legend.text = element_text(size = 12, face = 2))
  
  ggsave('3.Output/cblExpressUMAP.pdf', width = 4, height = 4.5)
  
  # LAG3 
  cblDf <- plotDf %>% 
    mutate(lag3 = round(lag3, 1))
  
  p <- ggplot()
  
  # Plotting points for each unique cbl value
  for (i in sort(unique(abs(cblDf$lag3)))) {
    p <- p + geom_point(data = cblDf %>% filter(abs(lag3) == i),
                        aes(umapscvi_1, umapscvi_2, col = lag3))
  }
  
  # Adding color scales and theme
  p + 
    scale_fill_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
    scale_color_gradient2(high = 'red', low = 'blue', mid = 'white', midpoint = 0) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = 'bottom',
          text = element_text(size = 12, face = 2),
          legend.text = element_text(size = 12, face = 2))
  
  ggsave('3.Output/lag3ExpressUMAP.pdf', width = 4, height = 4.5)
  
# 2. CBLs expression high/low groups split ----
# kmeans_result <- kmeans(plotDf$cbl, centers = 3)
# gene_df <- plotDf %>%
#     mutate(group = ifelse(kmeans_result$cluster == 1, "Low", "High"))

# test_result <- wilcox.test(cbl ~ group, data = gene_df)

# ggplot(gene_df, aes(cbl, fill = group)) +
#   geom_density(alpha = 0.5) +
#   annotate('text', x = Inf, y = Inf, 
#            label = 'High#: 20745, Low#: 36415 \n p-value < 2.2e-16', 
#            vjust = 2, hjust = 1.2, col = 'red',
#            size = 8/2.8, fontface = 2) +
#   labs(x = 'CBLs Expression',
#        y = 'Density') +
#   theme_bw() +
#   theme(legend.position = 'bottom',
#         aspect.ratio = 1,
#         text = element_text(face = 2, size = 8, color = 'black')) 

# gene_df  %>% 
#   ggplot(aes(x = group, y = lag3, fill = group)) +
#   ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, 
#                        justification = -.3, point_colour = NA, alpha = 0.5) + 
#   geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5, show.legend = F) +
#   labs(y = 'LAG3 Expression',
#        x = 'CBLs Expression Group') +
#   theme_bw() +
#   theme(legend.position = 'bottom',
#         aspect.ratio = 1,
#         text = element_text(face = 2, size = 8, color = 'black'))

# ggsave('3.Output/lag3ExpWithCBLs.png', width = 4, height = 4.5)

# ggsave('3.Output/cblHiLowSplit.png', width = 4, height = 4.5)

# write.csv(gene_df, '3.Output/workingDf.csv')
