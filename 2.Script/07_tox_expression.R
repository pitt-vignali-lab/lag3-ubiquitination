# ----------------------------------------------------------------------------
# Fig 7K | TOX expression across LAG3 / CBL patient groups in TOX-positive TILs
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/ToxExp.R
# Input    : 1.Data/cd8SetScvi.RDS
#            3.Output/workingDf_pid_cat.csv
# Output   : 3.Output/TOX_violin.png
# Upstream : 08_cluster_distribution_by_treatment.R
# Note     : TOX_violin.png is the TOX-positive panel (tox > 0) used in the figure
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(ggpubr)

obj <- readRDS('1.Data/cd8SetScvi.RDS')

objTil <- subset(obj, sample_type == 'TIL' & 
                     treatment %in% c('Baseline'))

pidCat <- read.csv('3.Output/workingDf_pid_cat.csv') 

# Add metadata and filter out cells with missing categories
objTil <- objTil %>%
  AddMetaData(metadata = pidCat %>% 
    select(pid, cat)  %>% 
    merge(data.frame(pid = objTil$patient_id) %>% rownames_to_column('cell'), by = 'pid') %>% 
    select(cell, cat) %>% 
    column_to_rownames('cell'))

objTil <- AddModuleScore(JoinLayers(objTil), 
                         features = 'TOX', 
                         name = 'TOX')

my_comparisons <- list( c('DN', 'LAG3_CBLs_Low'), 
                        c('DN', 'LAG3_CBLs_High'), 
                        c('LAG3_CBLs_Low', 'LAG3_CBLs_High') )

remove_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25)
  Q3 <- quantile(data[[column]], 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  data[data[[column]] >= lower_bound & data[[column]] <= upper_bound, ]
}

df = data.frame(cat = objTil$cat, tox = objTil$TOX1) %>% 
  filter(tox > 0) %>%
  mutate(cat = factor(cat, levels = c('DN', 'LAG3_CBLs_Low', 'LAG3_CBLs_High'))) 
  
df %>%
  ggplot(aes(x = cat, y = tox)) +
  geom_violin(aes(fill = cat), alpha = 0.5) +
  geom_boxplot(width = 0.1) + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_classic() +
  theme(aspect.ratio = 1,
        text = element_text(size = 8, face = 2)) 

ggsave('3.Output/TOX_violin.png', width = 6, height = 6)

df = data.frame(cat = objTil$cat, tox = objTil$TOX1) %>% 
  # filter(tox > 0) %>%
  mutate(cat = factor(cat, levels = c('DN', 'LAG3_CBLs_Low', 'LAG3_CBLs_High'))) 
  
df %>%
  ggplot(aes(x = cat, y = tox)) +
  geom_violin(aes(fill = cat), alpha = 0.5) +
  geom_boxplot(width = 0.01) + 
  geom_point(data = df %>% 
                   group_by(cat) %>%
                     summarise(median = median(tox)), 
            aes(x = cat, y = median), col = 'red') +
  stat_compare_means(comparisons = my_comparisons) +
  theme_classic() +
  theme(aspect.ratio = 1,
        text = element_text(size = 8, face = 2)) 

ggsave('3.Output/TOX_violin2.png', width = 6, height = 6)
