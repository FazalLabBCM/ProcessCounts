#!/usr/bin/env Rscript

  
#### PREPARE ENVIRONMENT ####

# Clear environment
rm(list = ls())

# Retrieve command line arguments
args = commandArgs(trailingOnly = TRUE)
if (interactive()) {
  scriptdir = "."
  outputdir = "."
  project_name = "test"
} else {
  scriptdir = args[1]
  outputdir = args[2]
  project_name = args[3]
}

# Import libraries
library(tidyverse)
library(DESeq2)
library(corrplot)
library(Rtsne)


#### DEFINE VARIABLES ####

# Folders
data_folder = "DataFiles/"
figures_folder = "Figures/"

# Files
genes_file = paste(scriptdir, "AdditionalFiles", "CommonGeneNames.txt", sep = "/")
lengths_file = paste(scriptdir, "AdditionalFiles", "TranscriptLengths.txt", sep = "/")

# Date
date = format(Sys.Date(), "%Y%m%d")

# Additional file prefix
deseq_prefix = "DESeq2"


#### PERFORM DESEQ ANALYSIS ####

target = "target"
control = "control"

for (file in list.files(path = outputdir, pattern = "^[[:upper:]]+-.*[[:digit:]]+C\\.txt", recursive = TRUE, full.names = TRUE)) {
  print(paste0("working with file ", file))
  base_name = str_sub(file, end = str_locate(file, pattern = "C\\.txt")[1])
  abbreviation = str_sub(base_name, 
                         start = str_locate(base_name, pattern = "/[:upper:]+-")[1] + 1, 
                         end = str_locate(base_name, pattern = "-[[:digit:]]+C$")[1] - 1)
  
  # Read data file
  data = read_tsv(file)
  gene = data$Ensembl_Gene
  data = select(data, -Ensembl_Gene)
  
  # Define column types (condition is either "target" or "control")
  num_targets = ncol(select(data, matches(".+-T[[:digit:]]$")))
  num_controls = ncol(select(data, matches(".+-C[[:digit:]]$")))
  condition = factor(c(rep(target, num_targets), rep(control, num_controls)))
  colTypes = tibble(condition)
  
  # DESeq analysis
  DESeq_data_set = DESeqDataSetFromMatrix(data, colTypes, formula(~condition))
  DESeq_data_set = estimateSizeFactors(DESeq_data_set)
  DESeq_data_set = DESeq(DESeq_data_set)
  DESeq_result = results(DESeq_data_set, alpha = 0.05)
  
  # Prepare data for Pearson correlation plot
  counts = as_tibble(counts(DESeq_data_set))
  sizeFactors = as_tibble(sizeFactors(DESeq_data_set)) %>%
    add_column("name" = names(sizeFactors(DESeq_data_set)), .before = 1)
  # Add targets counts and size factors to tibbles
  if (!exists("targets_counts")) {
    targets_counts = select(counts, matches(".+-T[[:digit:]]$"))
  } else {
    targets_counts = add_column(targets_counts, select(counts, matches(".+-T[[:digit:]]$")))
  }
  new_targets_sizeFactors = sizeFactors %>%
    filter(grepl("-T", name)) %>%
    pull(value)
  names(new_targets_sizeFactors) = sizeFactors %>%
    filter(grepl("-T", name)) %>%
    pull(name)
  if (!exists("targets_sizeFactors")) {
    targets_sizeFactors = new_targets_sizeFactors
  } else {
    targets_sizeFactors = c(targets_sizeFactors, new_targets_sizeFactors)
  }
  # Add controls counts and size factors to tibbles
  if (!exists("controls_counts")) {
    controls_counts = select(counts, matches(".+-C[[:digit:]]$"))
  } else {
    controls_counts = add_column(controls_counts, select(counts, matches(".+-C[[:digit:]]$")))
  }
  new_controls_sizeFactors = sizeFactors %>%
    filter(grepl("-C", name)) %>%
    pull(value)
  names(new_controls_sizeFactors) = sizeFactors %>%
    filter(grepl("-C", name)) %>%
    pull(name)
  if (!exists("controls_sizeFactors")) {
    controls_sizeFactors = new_controls_sizeFactors
  } else {
    controls_sizeFactors = c(controls_sizeFactors, new_controls_sizeFactors)
  }
  
  # Save DESeq results
  result_data = select(data, contains(abbreviation))
  result_data = add_column(result_data, 
                           "baseMean" = DESeq_result$baseMean, 
                           "log2FoldChange" = DESeq_result$log2FoldChange, 
                           "pValueAdjusted" = DESeq_result$padj)
  result_data = add_column(result_data, "Ensembl_Gene" = gene, .before = 1)
  file_name = str_sub(base_name, start = str_locate(base_name, pattern = abbreviation)[1])
  write_tsv(result_data, file = paste0(data_folder, "DESeq/", paste(deseq_prefix, file_name, sep = "_"), ".txt"))
}


#### MERGE AND SAVE ALL RESULTS ####

# Merge DESeq data into one data table
for (file in list.files(path = outputdir, pattern = "DESeq2_[[:upper:]]+-.*[[:digit:]]+C\\.txt", recursive = TRUE, full.names = TRUE)) {
  print(paste0("working with file ", file))
  base_name = str_sub(file, end = str_locate(file, pattern = "C\\.txt")[1])
  
  # Read file
  DESeq_data = read_tsv(file)
  if (!exists("all_DESeq_data")) {
    all_DESeq_data = select(DESeq_data, Ensembl_Gene)
  }
  
  # Select and rename columns
  DESeq_data = select(DESeq_data, c(Ensembl_Gene, log2FoldChange, pValueAdjusted))
  name = str_extract(base_name, pattern = "[:upper:]+-.*[:digit:]+C")
  colNameFC = str_c(name, "-FC")
  colNameP = str_c(name, "-P")
  DESeq_data = setNames(DESeq_data, c("Ensembl_Gene", colNameFC, colNameP))
  
  # Join data
  all_DESeq_data = inner_join(all_DESeq_data, DESeq_data)
}
# Add common gene names
genes = read_tsv(genes_file)
genes_joined = left_join(all_DESeq_data, genes, by = "Ensembl_Gene")
all_DESeq_data = add_column(all_DESeq_data, "Common_Gene" = genes_joined$Common_Gene, .after = 1)
# Save all DESeq data to file
write_tsv(all_DESeq_data, file = paste0(data_folder, paste(date, project_name, "DESeqResults.txt", sep = "_")))


#### CREATE PEARSON CORRELATION PLOT ####

filter_counts_tibble = function(counts_tibble) {
  min_avg_count = (ncol(counts_tibble))*1000
  counts_mat = as.matrix(counts_tibble)
  counts_sums = as.matrix(rowSums(counts_mat))
  counts_mat_filtered = counts_mat[(which(counts_sums > min_avg_count)),]
  return (counts_mat_filtered)
}

prepare_cor_mat = function(counts_tibble, sizeFactors_vector) {
  counts_table = as.data.frame(filter_counts_tibble(counts_tibble))
  counts_norm = t(t(counts_table) / sizeFactors_vector)
  cor_mat = cor(counts_norm, method = "pearson")
  # Reorder using correlation as distance
  dd = as.dist((1 - cor_mat) / 2)
  hc = hclust(dd)
  cor_mat = cor_mat[hc$order, hc$order]
  return(cor_mat)
}

create_cor_plot = function(cor_mat, targets_or_controls, colors) {
  if (targets_or_controls == "Controls") {
    num_squares = NULL
  } else {
    num_squares = floor(ncol(cor_mat) / 3)
    if (num_squares == 0) {
      num_squares = 1
    }
  }
  corrplot(cor_mat, 
           method="circle", 
           order = "hclust", 
           col.lim = c(0, 1), 
           tl.col = "black", 
           hclust.method = "centroid", 
           addrect = num_squares, 
           col = colors(40))
}

create_and_save_cor_plot = function(cor_mat, targets_or_controls) {
  colors = colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F", "cyan", "#007FFF", 
                               "blue", "#00007F","#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", 
                               "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  # Save as PDF
  pdf(file = paste0(figures_folder, paste(date, project_name, targets_or_controls, "PearsonCorrelationPlot.pdf", sep = "_")), 
      width = 7, height = 7)
  create_cor_plot(cor_mat, targets_or_controls, colors)
  dev.off()
  # Save as PNG
  png(file = paste0(figures_folder, paste(date, project_name, targets_or_controls, "PearsonCorrelationPlot.png", sep = "_")), 
      width = 7, height = 7, units = "in", res = 300)
  create_cor_plot(cor_mat, targets_or_controls, colors)
  dev.off()
}

# Targets
targets_cor_mat = prepare_cor_mat(targets_counts, targets_sizeFactors)
create_and_save_cor_plot(targets_cor_mat, "Targets")
# Controls
controls_cor_mat = prepare_cor_mat(controls_counts, controls_sizeFactors)
create_and_save_cor_plot(controls_cor_mat, "Controls")
# Targets and controls
targets_and_controls_counts = cbind(targets_counts, controls_counts)
targets_and_controls_sizeFactors = c(targets_sizeFactors, controls_sizeFactors)
targets_and_controls_cor_mat = prepare_cor_mat(targets_and_controls_counts, targets_and_controls_sizeFactors)
create_and_save_cor_plot(targets_and_controls_cor_mat, "TargetsAndControls")


#### CREATE TSNE PLOT ####

create_and_save_tsne_plot = function(counts_tibble, targets_or_controls) {
  # Transform data
  vst_mat = rlogTransformation(filter_counts_tibble(counts_tibble))
  # Make PCA plot of samples using rlog transformed data
  vst_PCA = prcomp(t(vst_mat))
  
  # Perform tSNE analysis
  d = dist(t(vst_mat))
  set.seed(40)
  max_perplexity = floor((attributes(d)$Size - 1) / 3)
  if (max_perplexity < 1) {
    cat("Too few", targets_or_controls, "samples to create tSNE plot.", "\n")
    return(0)
  }
  tsne_out = Rtsne(d, is_distance = TRUE, perplexity = max_perplexity, verbose = TRUE) 
  combined_tsne = as.data.frame(cbind(tsne_out$Y, rownames(t(vst_mat))))
  colnames(combined_tsne) = c("SNE_1","SNE_2","name")
  
  combined_tsne = combined_tsne %>%
    mutate("SNE_1" = as.double(SNE_1), "SNE_2" = as.double(SNE_2)) %>%
    mutate(color = ifelse(grepl("-T[[:digit:]]$", combined_tsne$name), yes = "Target", no = "Control"))
  
  # Create tSNE plot
  tsne_plot = ggplot(data = combined_tsne, aes(x = SNE_1, y = SNE_2, label = name)) + 
    geom_point(size=5, alpha = 0.3, aes(color = factor(color))) + 
    geom_text(size = 4)  +
    theme_bw() +
    coord_cartesian(clip = "off") +
    scale_color_manual(guide = "none", values = c("Target" = "red4", "Control" = "dodgerblue4")) +
    labs(x = "tSNE 2", y = "tSNE 1") +
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"), 
          plot.margin = unit(c(3,3,1,1), "lines"))
  
  # Save as PDF
  pdf(file = paste0(figures_folder, paste(date, project_name, targets_or_controls, "tSNEPlot.pdf", sep = "_")), 
      width = 7, height = 7)
  print(tsne_plot)
  dev.off()
  # Save as PNG
  png(file = paste0(figures_folder, paste(date, project_name, targets_or_controls, "tSNEPlot.png", sep = "_")), 
      width = 7, height = 7, units = "in", res = 300)
  print(tsne_plot)
  dev.off()
}

# Targets
create_and_save_tsne_plot(targets_counts, "Targets")
# Controls
create_and_save_tsne_plot(controls_counts, "Controls")
# Targets and controls
create_and_save_tsne_plot(targets_and_controls_counts, "TargetsAndControls")


#### CLEAR ENVIRONMENT ####
rm(list = ls())
