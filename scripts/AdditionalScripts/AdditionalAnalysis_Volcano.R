#!/usr/bin/env Rscript


#### PREPARE ENVIRONMENT ####
rm(list = ls())
library(tidyverse)
# Set global variables
MIN_TRANSCRIPT_LENGTH = 100
MIN_COUNTS = 1


#### IMPORT DATA ####
deseq_data_file = "/path/to/DataFiles/DATE_ProjectName_DESeq2Results.txt"  # EDIT: file path and file name for your experiment's DESeq2 results
deseq_data = read_tsv(deseq_data_file) 
htseq_counts_file = "/path/to/DataFiles/DATE_ProjectName_HTSeqCounts.txt"  # EDIT: file path and file name for your experiment's HTSeq counts
htseq_counts = read_tsv(htseq_counts_file)
transcript_lengths_file = "/path/to/reference/CommonGeneNamesAndTranscriptLengths/TranscriptLengths.txt"  # EDIT: file path and file name for the transcript lengths reference table
transcript_lengths = read_tsv(transcript_lengths_file)
true_positive_genes_file = "/path/to/reference/LocalizedGenes/LOCATION_Localized_Genes_Table.txt"  # EDIT: file path and file name for the reference table containing the genes you expect to be significantly enriched in your experiment (replace next line with "true_positive_genes = c()" to skip this part)
true_positive_genes = read_tsv(true_positive_genes_file) %>% pull(Common_Gene)


#### PREPARE DATA ####

# Pivot data table
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

# Filter out genes with short transcript length
genes_filtered_min_transcript_length = deseq_data %>%
  filter(transcript_lengths$Transcript_Length >= MIN_TRANSCRIPT_LENGTH) %>%
  pull(Ensembl_Gene)
volcano_data = volcano_data %>%
  filter(Ensembl_Gene %in% genes_filtered_min_transcript_length)

# Filter out genes with dropoff or no counts in any condition
genes_filtered_min_counts = htseq_counts %>%
  filter(if_all(contains("-"), ~ . >= MIN_COUNTS)) %>%
  pull(Ensembl_Gene)
volcano_data = volcano_data %>%
  filter(Ensembl_Gene %in% genes_filtered_min_counts)


#### VOLCANO PLOT ####

# Make volcano plot
volcano_plot = ggplot() +
  geom_point(data = subset(volcano_data, Expected_Positive == FALSE), 
             aes(x = Fold_Change, y = -log10(P_Value), color = Expected_Positive), alpha = 1/4, size = 2) +
  geom_point(data = subset(volcano_data, Expected_Positive == TRUE), 
             aes(x = Fold_Change, y = -log10(P_Value), color = Expected_Positive), alpha = 1/4, size = 2) +
  facet_grid(cols = vars(Condition)) +
  labs(x = "Fold Change", y = bquote(-log[10](pvalue))) +
  scale_color_manual(values = c("grey50", "darkred")) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 10)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        plot.margin = unit(c(3, 3, 3, 3), "lines"))
print(volcano_plot)

# Calculate file width
num_conditions = volcano_data %>%
  group_by(Condition) %>%
  summarize() %>%
  nrow()
w = if_else(num_conditions < 4, 7, num_conditions * 2)  # width is 7in unless there are 4+ conditions

# Save volcano plot
pdf(file = "/path/to/Figures/DATE_ProjectName_VolcanoPlot.pdf", width = w, height = 5)  # EDIT: file path and file name for output volcano plot PDF
print(volcano_plot)
dev.off()
png(file = "/path/to/Figures/DATE_ProjectName_VolcanoPlot.png", width = w, height = 5, units = "in", res = 300)  # EDIT: file path and file name for output volcano plot PNG
print(volcano_plot)
dev.off()

