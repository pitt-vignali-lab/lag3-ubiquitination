# ----------------------------------------------------------------------------
# Fig 7H | Patient-wise Pearson correlation between LAG3 and CBL expression
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/2.2.corDotRelaCombo.R
# Input    : 1.Data/cd8SetScvi.RDS
# Output   : 3.Output/pre_post_treated_patient_wised_cor_lag3_cbls_TIL_rela_combo.pdf
# Upstream : 02_cd8_integration_clustering.R
# Note     : P values adjusted by FDR
#            red dots mark FDR < 0.05
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

obj <- readRDS('1.Data/cd8SetScvi.RDS')

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

# 1. patient wised lag3 cbls expression correlation test ----
corRes <- plotDf %>%
  filter(tissue == 'TIL') %>% 
  group_by(pid, tissue, treat) %>%
  summarise(
    cor_value = cor.test(lag3, cbl)$estimate,
    p_value = cor.test(lag3, cbl)$p.value
  ) %>% 
  mutate(fdr = p.adjust(p_value),
         treat = factor(treat, levels = c('Baseline', 'rela', 'combo'))) # %>%
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

ggsave('3.Output/pre_post_treated_patient_wised_cor_lag3_cbls_TIL_rela_combo.pdf', width = 10, height = 5.5)
