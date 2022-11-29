#!/usr/bin/env Rscript


#### PREPARE ENVIRONMENT ####
rm(list = ls())
library(tidyverse)
# Set global variables
P_SIGNIFICANT = 0.05


#### IMPORT DATA ####
deseq_data_file = "/path/to/DataFiles/DATE_ProjectName_DESeq2Results.txt"  # EDIT: file path and file name for your experiment's DESeq2 results
deseq_data = read_tsv(deseq_data_file)


#### GET TRUE POSITIVE AND TRUE NEGATIVE GENES ####
true_positive_genes = c("GENEID1", "GENEID2", "ETC.")  # EDIT: vector of Ensembl_Gene IDs for genes expected to be significantly enriched
true_negative_genes = c("GENEID3", "GENEID4", "ETC.")  # EDIT: vector of Ensembl_Gene IDs for genes expected to NOT be significantly enriched


#### PREPARE DATA ####
# Pivot longer
deseq_data_fc = deseq_data %>%
  select(c(Ensembl_Gene, ends_with("C"))) %>%
  pivot_longer(ends_with("C"), names_to = "Condition", values_to = "Fold_Change") %>%
  mutate("Condition" = str_sub(Condition, end = -4))
deseq_data_p = deseq_data %>%
  select(c("Ensembl_Gene", ends_with("P"))) %>%
  pivot_longer(ends_with("P"), names_to = "Condition", values_to = "P_Value") %>%
  mutate("Condition" = str_sub(Condition, end = -3))
deseq_data = inner_join(deseq_data_fc, deseq_data_p, by = c("Ensembl_Gene", "Condition"))
# Label true positive and true negative genes
deseq_data = deseq_data %>%
  filter(Ensembl_Gene %in% true_positive_genes | Ensembl_Gene %in% true_negative_genes, 
         P_Value <= P_SIGNIFICANT) %>%
  mutate("True_Positive" = if_else(Ensembl_Gene %in% true_positive_genes, TRUE, FALSE), 
         "True_Negative" = !True_Positive) %>%
  group_by(Condition) %>%
  arrange(desc(Fold_Change)) %>%
  drop_na() %>%
  # Calculate True Positive Rate and False Positive Rate
  mutate("TPR" = cumsum(True_Positive) / sum(True_Positive), 
         "FPR" = cumsum(True_Negative) / sum(True_Negative)) %>%
  # Calculate AUC
  mutate("FPR_Increase" = FPR - c(0, FPR[1:(length(FPR) - 1)]), 
         "Partition_Area" = TPR * FPR_Increase,  # TPR is height of rectangle, FPR_Increase is width of rectangle 
         "Cumulative_Sum_Of_Partitions" = cumsum(Partition_Area)) %>%
  # Calculate cutoff
  mutate("TPR-FPR" = TPR - FPR)

cutoffs = deseq_data %>%
  #filter(FPR <= 0.2) %>%
  summarize("Max_TPR-FPR" = max(`TPR-FPR`)) %>%
  add_column("Fold_Change" = 0, 
             "Num_TP" = 0, 
             "Num_FP" = 0, 
             "TPR" = 0, 
             "FPR" = 0)
for (i in 1:length(cutoffs$Condition)) {
  condition = cutoffs$Condition[i]
  cutoff = cutoffs %>%
    filter(Condition == condition) %>%
    pull(`Max_TPR-FPR`)
  roc_data_at_cutoff = roc_data %>%
    filter(Condition == condition, 
           `TPR-FPR` == cutoff)
  cutoffs$Fold_Change[i] = roc_data_at_cutoff$Fold_Change
  cutoffs$Num_TP[i] = roc_data_at_cutoff$Num_TP
  cutoffs$Num_FP[i] = roc_data_at_cutoff$Num_FP
  cutoffs$TPR[i] = roc_data_at_cutoff$TPR
  cutoffs$FPR[i] = roc_data_at_cutoff$FPR
}

areas_under_curve = deseq_data %>%
  summarize("AUC" = max(Cumulative_Sum_Of_Partitions))


#### PLOT ROC ####
ggplot(deseq_data, aes(x = FPR, y = TPR, color = Condition)) +
  geom_step(size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "False Positive Rate", 
       y = "True Positive Rate") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())
