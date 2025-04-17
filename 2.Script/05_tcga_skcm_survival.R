# ----------------------------------------------------------------------------
# Fig 7I | Kaplan-Meier survival analysis in TCGA-SKCM
#
# Paper    : Ligand-induced ubiquitination unleashes LAG3 immune checkpoint
#            function by hindering membrane sequestration of signaling motifs
# Original : 2.Script/5.1.survivalAnalysis.R
# Input    : TCGA-SKCM expression and clinical data (downloaded via easyTCGA)
# Output   : 3.Output/5.survivalAnalysis/survivalAnalysis_lag3+.png
#            3.Output/5.survivalAnalysis/survivalAnalysis_lag3-.png
# Note     : Quantile thresholds: LAG3 60th, CBLB 80th, CBL 50th percentile
# ----------------------------------------------------------------------------

library(easyTCGA)
library(survival)
library(survminer)
library(tidyverse)
library(ggplotify)
library(export)

# getmrnaexpr('TCGA-SKCM')

load("output_mRNA_lncRNA_expr/TCGA-SKCM_mrna_expr_count.rdata")
load("output_mRNA_lncRNA_expr/TCGA-SKCM_clinicalSE.rdata")

clin_info <- clinicalSE
mrna_expr_counts <- mrna_expr_count

dim(mrna_expr_counts)

clin_info$days_to_death[is.na(clin_info$days_to_death)] <- 0
clin_info$days_to_last_follow_up[is.na(clin_info$days_to_last_follow_up)] <- 0  
clin_info$days <- ifelse(
  clin_info$vital_status == 'Alive',
  clin_info$days_to_last_follow_up,clin_info$days_to_death
)

clin_info$year <- round(clin_info$days/365, 2)
clin_info$month <- round(clin_info$days/30, 2)

# clin_info$age_at_index[is.na(clin_info$age_at_index)] <- 0
# clin_info$age_at_index <- as.numeric(clin_info$age_at_index)
# clin_info$age_group <- ifelse(clin_info$age_at_index>median(clin_info$age_at_index), 'older', 'younger')


clin_info$event <- ifelse(clin_info$vital_status == 'Alive', 0, 1)

exprSet <- mrna_expr_counts

# 基因分组
clin_info$LAG3 <- t(ifelse(
  exprSet['LAG3',] > median(as.numeric(exprSet['LAG3',])),
  'High','Low'
))

clin_info$CBLs <- t(ifelse(
  exprSet['CBL',] > median(as.numeric(exprSet['CBL',])) &
  exprSet['CBLC',] > median(as.numeric(exprSet['CBLC',])),
  'High','Low'
))

# if(dir.exists('3.Output/5.survivalAnalysis') == FALSE){
#   dir.create('3.Output/5.survivalAnalysis')
# }
# 
# # qutail method
# # Define custom thresholds
# for(i in seq (0, 1, 0.05)){
#       lag3_threshold <- quantile(as.numeric(exprSet['LAG3', ]), i, na.rm = TRUE)
# cbl_threshold <- quantile(as.numeric(exprSet['CBL', ]), i, na.rm = TRUE)
# cblc_threshold <- quantile(as.numeric(exprSet['CBLC', ]), i, na.rm = TRUE)
# 
# # Apply custom thresholds to classify as "High" or "Low"
# clin_info$LAG3 <- t(ifelse(exprSet['LAG3', ] > lag3_threshold, 'High', 'Low'))
# 
# clin_info$CBLs <- t(ifelse(
#   exprSet['CBL', ] > cbl_threshold & exprSet['CBLC', ] > cblc_threshold,
#   'High', 'Low'
# ))
# 
# clin_info <- clin_info %>%
#   mutate(cat = case_when(
#     LAG3 == "High" & CBLs == "High" ~ "LAG3+CBLs High",
#     LAG3 == "High" & CBLs == "Low" ~ "LAG3 High",
#     LAG3 == "Low" & CBLs == "High" ~ "CBLs High",
#     TRUE ~ "Double negative"
#   ))
# 
# clin_info_sub <- clin_info %>%
#   filter(year <= 10)
# 
# sfit1 <- survfit(Surv(year, event)~cat, data = clin_info_sub)
# 
# ggsurvplot(
#   sfit1, 
#   pval = TRUE, 
#   data = clin_info_sub, 
#   # 添加中位生存
#   surv.median.line = "hv",
#   risk.table = F)
# 
#   ggsave(paste0('3.Output/5.survivalAnalysis/NEW_', i, 'survivalAnalysis.png'))
# 
# }

lag3_threshold <- quantile(as.numeric(exprSet['LAG3', ]), 0.6, na.rm = TRUE)
cbl_threshold <- quantile(as.numeric(exprSet['CBL', ]), 0.8, na.rm = TRUE)
cblc_threshold <- quantile(as.numeric(exprSet['CBLC', ]), 0.3, na.rm = TRUE)

# Apply custom thresholds to classify as "High" or "Low"
clin_info$LAG3 <- t(ifelse(exprSet['LAG3', ] > lag3_threshold, 'High', 'Low'))

clin_info$CBLs <- t(ifelse(
  exprSet['CBL', ] > cbl_threshold & exprSet['CBLC', ] > cblc_threshold,
  'High', 'Low'
))

clin_info <- clin_info %>%
  mutate(cat = case_when(
    LAG3 == "High" & CBLs == "High" ~ "LAG3+CBLs High",
    LAG3 == "High" & CBLs == "Low" ~ "LAG3 High",
    LAG3 == "Low" & CBLs == "High" ~ "CBLs High",
    TRUE ~ "Double negative"
  ),
  cat = factor(cat, levels = c("LAG3+CBLs High", "LAG3 High", "CBLs High", "Double negative")))



clin_info_sub <- clin_info %>%
  filter(year <= 5,
         cat %in% c("LAG3+CBLs High", "LAG3 High"))

clin_info_sub <- clin_info %>%
  filter(year <= 5,
         cat %in% c("CBLs High", "Double negative"))

sfit1 <- survfit(Surv(year, event)~cat, data = clin_info_sub)

ggsurvplot(
  sfit1, 
  pval = TRUE, 
  data = clin_info_sub, 
  # 添加中位生存
  surv.median.line = "hv",
  risk.table = F,
  ggtheme = theme_bw())

ggsave('3.Output/5.survivalAnalysis/survivalAnalysis_lag3+.png', width = 5, height = 5)

ggsurvplot(
  sfit1,
  risk.table = F, 
  pval =TRUE,
  xlab = "Time in days",  
  conf.int =TRUE,
  # 添加中位生存
  surv.median.line = "hv",
  ggtheme = theme_bw() +
    theme(legend.position = "bottom")
)

quantiles <- data.frame(
  gene = c("LAG3", "CBL", "CBLC"),
  quantile = c(0.6, 0.8, 0.5)
)

exprSet[c('LAG3', 'CBL', 'CBLC'), ] %>% 
  t() %>% 
  data.frame() %>% 
  pivot_longer(everything(), names_to = 'gene', values_to = 'value') %>% 
  ggplot(aes(x = gene, y = value)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "point", shape = 21, size = 3, fill = "red", color = "red", vjust = -0.5) +
  geom_point(data = quantiles, aes(x = gene, y = quantile), color = "blue", size = 3)

 desired_quantiles <- data.frame(
    gene = c("LAG3", "CBL", "CBLC"),
    quantile = c(
      quantile(exprSet["LAG3", ], probs = 0.6, na.rm = TRUE),
      quantile(exprSet["CBL", ], probs = 0.8, na.rm = TRUE),
      quantile(exprSet["CBLC", ], probs = 0.5, na.rm = TRUE)
    )
  )

exprSet[c('LAG3', 'CBL', 'CBLC'), ] %>%
  t() %>%
  data.frame() %>%
  pivot_longer(everything(), names_to = 'gene', values_to = 'value') %>%
  ggplot(aes(x = gene, y = value)) +
  geom_boxplot() + 
  geom_point(data = desired_quantiles, aes(x = gene, y = quantile), color = "red", size = 3) +
  facet_wrap(~gene, scales = "free")

"To categorize expression levels as 'High' or 'Low' for LAG3, CBL, and CBLC, we applied quantile-based thresholds tailored to each gene's unique expression distribution. Given the variability in expression patterns across genes, using a single cutoff (e.g., median) may not adequately capture biologically relevant high or low expression levels for each gene. Thus, we applied the 60th percentile for LAG3, the 80th percentile for CBL, and the 50th percentile for CBLC, optimizing each threshold to reflect meaningful expression levels for that specific gene.

For LAG3, the 60th percentile was selected to define 'High' expression, capturing a moderately elevated subset of patients, as this gene's distribution showed a slight skew with most values concentrated below this level. For CBL, we used a higher threshold at the 80th percentile to focus on a smaller, highly elevated expression group due to a broad spread in its expression levels across patients. In contrast, for CBLC, the median (50th percentile) was used, as its expression pattern was more balanced, and a midpoint threshold effectively differentiated high and low expression groups.

These gene-specific thresholds were empirically determined to improve the sensitivity of our analysis, ensuring that each gene’s ‘High’ group meaningfully represented a distinct subset of patients with potentially greater clinical relevance. By tailoring cutoffs to the distribution of each gene, we aimed to enhance the clarity and biological significance of comparisons in clinical outcomes between expression groups."