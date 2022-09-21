#!/usr/bin/env Rscript


#### PREPARE ENVIRONMENT ####
rm(list = ls())
library(tidyverse)
library(pROC)


#### GET TRUE POSITIVE AND TRUE NEGATIVE GENES ####
true_positive_genes = c("GENENAME1", "GENENAME2", "ETC.")  # EDIT: vector of Common_Gene names for genes expected to be significantly enriched
true_negative_genes = c("GENENAME3", "GENENAME4", "ETC.")  # EDIT: vector of Common_Gene names for genes expected to NOT be significantly enriched


#### IMPORT EXPERIMENT DATA FILES ####
deseq_data_file = "/path/to/DataFiles/DATE_ProjectName_DESeq2Results.txt"  # EDIT: file path and file name for your experiment's DESeq2 results
deseq_data = read_tsv(deseq_data_file) 
htseq_counts_file = "/path/to/DataFiles/DATE_ProjectName_HTSeqCounts.txt"  # EDIT: file path and file name for your experiment's HTSeq counts
htseq_counts = read_tsv(htseq_counts_file)


#### ROC CURVE ####

# Prepare data
roc_data = deseq_data %>%
  filter(Common_Gene %in% true_positive_genes | Common_Gene %in% true_negative_genes) %>%
  mutate("Actual" = ifelse(Common_Gene %in% true_positive_genes, 1, 0)) %>%
  drop_na(Actual, `YourSample-FC`, `YourSample-P`)  # EDIT: the column names in your DESeq2 results table for the FC and P values to be used in your ROC curve

# Calculate fit
fit = glm(roc_data$Actual ~ roc_data$`YourSample-FC` + roc_data$`YourSample-P`,  # EDIT: the column names in your DESeq2 results table for the FC and P values to be used in your ROC curve 
          family = binomial)

# Plot and save ROC curve
par(pty = "s", cex.axis = 1.2, cex.lab = 1.5, cex = 1.5)
pdf(file = "/path/to/Figures/DATE_ProjectName_ROCcurve.pdf", width = 4, height = 4)  # EDIT: file path and file name for output ROC curve
roc(roc_data$Actual, 
    fit$fitted.values, 
    plot = TRUE, 
    legacy.axes = TRUE, 
    xlab = "False Positive Rate", 
    ylab = "True Positive Rate", 
    col = "#377eb8", 
    lwd = 3, 
    print.auc = TRUE, 
    print.auc.y = 0.2)
dev.off()


#### VOLCANO PLOT ####

# Prepare data
deseq_data_fc = deseq_data %>%
  select(c("Ensembl_Gene", "Common_Gene", ends_with("-FC"))) %>%
  pivot_longer(ends_with("-FC"), names_to = "Condition", values_to = "Fold_Change") %>%
  mutate("Condition" = str_sub(Condition, end = -4))
deseq_data_p = deseq_data %>%
  select(c("Ensembl_Gene", "Common_Gene", ends_with("-P"))) %>%
  pivot_longer(ends_with("-P"), names_to = "Condition", values_to = "P_Value") %>%
  mutate("Condition" = str_sub(Condition, end = -3))
volcano_data = inner_join(deseq_data_fc, deseq_data_p, by = c("Ensembl_Gene", "Common_Gene", "Condition")) %>%
  drop_na()

# Label genes expected to be enriched (true positive)
volcano_data = mutate(volcano_data, "Expected_Positive" = Common_Gene %in% true_positive_genes)

# Filter out genes with low counts
MIN_COUNTS = 10
genes_filtered = htseq_counts %>%
  filter(if_all(contains("-"), ~ . > MIN_COUNTS)) %>%
  pull("Common_Gene")
volcano_data = volcano_data %>%
  filter(Common_Gene %in% genes_filtered)

# Make volcano plot
volcano_plot = ggplot() +
  geom_point(data = subset(volcano_data, Expected_Positive == FALSE), 
             aes(x = Fold_Change, y = -log10(P_Value), color = Expected_Positive), alpha = 1/4, size = 2) +
  geom_point(data = subset(volcano_data, Expected_Positive == TRUE), 
             aes(x = Fold_Change, y = -log10(P_Value), color = Expected_Positive), alpha = 1/4, size = 2) +
  facet_grid(cols = vars(Condition)) +
  labs(x = "Fold Change", y = bquote(-log[10](pvalue))) +
  scale_color_manual(values = c("grey50", "red4")) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 10)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        plot.margin = unit(c(3, 3, 3, 3), "lines"))
print(volcano_plot)

# Save volcano plot
pdf(file = "/path/to/Figures/DATE_ProjectName_VolcanoPlot.pdf", width = 7, height = 5)  # EDIT: file path and file name for output volcano plot
print(volcano_plot)
dev.off()
png(file = "/path/to/Figures/DATE_ProjectName_VolcanoPlot.png", width = 7, height = 5, units = "in", res = 300)  # EDIT: file path and file name for output volcano plot
print(volcano_plot)
dev.off()

