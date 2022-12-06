#!/usr/bin/env Rscript


#### PREPARE ENVIRONMENT ####
rm(list = ls())
library(tidyverse)
# Set global variables
P_SIGNIFICANT = 0.05
MIN_TRANSCRIPT_LENGTH = 100
MIN_COUNTS = 1


#### IMPORT DATA ####
deseq_data_file = "/path/to/DataFiles/DATE_ProjectName_DESeq2Results.txt"  # EDIT: file path and file name for your experiment's DESeq2 results
deseq_data = read_tsv(deseq_data_file)
htseq_counts_file = "/path/to/DataFiles/DATE_ProjectName_HTSeqCounts.txt"  # EDIT: file path and file name for your experiment's HTSeq counts
htseq_counts = read_tsv(htseq_counts_file)
transcript_lengths_file = "/path/to/reference/CommonGeneNamesAndTranscriptLengths/TranscriptLengths.txt"  # EDIT: file path and file name for the transcript lengths reference table
transcript_lengths = read_tsv(transcript_lengths_file)
true_positive_genes_file = "/path/to/reference/LocalizedGenes/LOCATION_Localized_Genes_Table.txt"  # EDIT: file path and file name for the reference table containing the genes you expect to be significantly enriched in your experiment
true_positive_genes = read_tsv(true_positive_genes_file) %>% pull(Ensembl_Gene)
true_negative_genes = c("GENEID3", "GENEID4", "ETC.")  # EDIT: vector of Ensembl_Gene IDs for genes expected to NOT be significantly enriched


#### PREPARE DATA ####

# Pivot data table
deseq_data_fc = deseq_data %>%
  select(c(Ensembl_Gene, ends_with("C"))) %>%
  pivot_longer(ends_with("C"), names_to = "Condition", values_to = "Fold_Change") %>%
  mutate("Condition" = str_sub(Condition, end = -4))
deseq_data_p = deseq_data %>%
  select(c("Ensembl_Gene", ends_with("P"))) %>%
  pivot_longer(ends_with("P"), names_to = "Condition", values_to = "P_Value") %>%
  mutate("Condition" = str_sub(Condition, end = -3))
roc_data = inner_join(deseq_data_fc, deseq_data_p, by = c("Ensembl_Gene", "Condition"))

# Filter out genes with short transcript length
genes_filtered_min_transcript_length = deseq_data %>%
  filter(transcript_lengths$Transcript_Length >= MIN_TRANSCRIPT_LENGTH) %>%
  pull(Ensembl_Gene)
roc_data = roc_data %>%
  filter(Ensembl_Gene %in% genes_filtered_min_transcript_length)

# Filter out genes with dropoff or no counts in any condition
genes_filtered_min_counts = htseq_counts %>%
  filter(if_all(contains("-"), ~ . >= MIN_COUNTS)) %>%
  pull(Ensembl_Gene)
roc_data = roc_data %>%
  filter(Ensembl_Gene %in% genes_filtered_min_counts)

# Label true positive and true negative genes
roc_data = roc_data %>%
  filter(Ensembl_Gene %in% true_positive_genes | Ensembl_Gene %in% true_negative_genes, 
         P_Value <= P_SIGNIFICANT) %>%
  mutate("True_Positive" = if_else(Ensembl_Gene %in% true_positive_genes, TRUE, FALSE), 
         "False_Positive" = !True_Positive) %>%
  group_by(Condition) %>%
  arrange(desc(Fold_Change)) %>%
  drop_na() %>%
  
  # Calculate True Positive Rate and False Positive Rate
  mutate("Num_TP" = cumsum(True_Positive), 
         "Num_FP" = cumsum(False_Positive), 
         "TPR" = Num_TP / sum(True_Positive), 
         "FPR" = Num_FP / sum(False_Positive)) %>%
  
  # Calculate AUC
  mutate("FPR_Increase" = FPR - c(0, FPR[1:(length(FPR) - 1)]), 
         "Partition_Area" = TPR * FPR_Increase,  # TPR is height of rectangle, FPR_Increase is width of rectangle 
         "Cumulative_Sum_Of_Partitions" = cumsum(Partition_Area)) %>%
  
  # Calculate cutoff
  mutate("TPR-FPR" = TPR - FPR)

# Make table for cutoff data
cutoffs = roc_data %>%
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

# Make table for AUC data
areas_under_curve = roc_data %>%
  summarize("AUC" = max(Cumulative_Sum_Of_Partitions))


#### ROC PLOT ####

# Make ROC
roc_plot = ggplot(deseq_data, aes(x = FPR, y = TPR, color = Condition)) +
  geom_step(size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "False Positive Rate", 
       y = "True Positive Rate") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())
print(roc_plot)

# Save ROC
pdf(file = "/path/to/Figures/DATE_ProjectName_ROC.pdf", width = 7, height = 5)  # EDIT: file path and file name for output ROC plot PDF
print(roc_plot)
dev.off()
png(file = "/path/to/Figures/DATE_ProjectName_ROC.png", width = 7, height = 5, units = "in", res = 300)  # EDIT: file path and file name for output ROC plot PNG
print(roc_plot)
dev.off()

