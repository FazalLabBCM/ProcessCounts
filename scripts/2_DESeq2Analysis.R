#!/usr/bin/env Rscript

  
#### PREPARE ENVIRONMENT ####

# Clear environment
rm(list = ls())

# Retrieve command line arguments
args = commandArgs(trailingOnly = TRUE)
if (interactive()) {
  SCRIPTDIR = "."
  OUTPUTDIR = "."
  PROJECT_NAME = "test"
  COMBINE_CONTROLS = "FALSE"
} else {
  SCRIPTDIR = args[1]
  OUTPUTDIR = args[2]
  PROJECT_NAME = args[3]
  COMBINE_CONTROLS = args[4]
}

# Import libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Rtsne))


#### DEFINE VARIABLES ####

# Folders
data_folder = paste(OUTPUTDIR, "DataFiles/", sep = "/")
figures_folder = paste(OUTPUTDIR, "Figures/", sep = "/")

# Files
genes_file = paste(SCRIPTDIR, "AdditionalFiles", "CommonGeneNames.txt", sep = "/")
lengths_file = paste(SCRIPTDIR, "AdditionalFiles", "TranscriptLengths.txt", sep = "/")

# Date
date = format(Sys.Date(), "%Y%m%d")

# Additional file prefix
deseq_prefix = "DESeq2"


#### PERFORM DESEQ ANALYSIS ####

target = "target"
control = "control"

target_pattern = ".+-T[[:digit:]]$"
control_pattern = ".+-C[[:digit:]]$"

for (file in list.files(path = data_folder, pattern = "^[^/_]+-[^/_]*[[:digit:]]+C\\.txt", recursive = TRUE, full.names = TRUE)) {
  cat("Working with file ", file, "\n")
  base_name = str_sub(file, end = str_locate(file, pattern = "C\\.txt")[1])
  abbreviation = str_sub(base_name, 
                         start = str_locate(base_name, pattern = "/[^/_]+-")[1] + 1, 
                         end = str_locate(base_name, pattern = "-[[:digit:]]+C$")[1] - 1)
  
  # Read data file
  data = read_tsv(file, show_col_types = FALSE)
  gene = data$Ensembl_Gene
  data = select(data, -Ensembl_Gene)
  
  # Define column types (condition is either "target" or "control")
  num_targets = ncol(select(data, matches(target_pattern)))
  num_controls = ncol(select(data, matches(control_pattern)))
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
    targets_counts = select(counts, matches(target_pattern))
  } else {
    targets_counts = add_column(targets_counts, select(counts, matches(target_pattern)))
  }
  new_targets_sizeFactors = sizeFactors %>%
    filter(grepl(target_pattern, name)) %>%
    pull(value)
  names(new_targets_sizeFactors) = sizeFactors %>%
    filter(grepl(target_pattern, name)) %>%
    pull(name)
  if (!exists("targets_sizeFactors")) {
    targets_sizeFactors = new_targets_sizeFactors
  } else {
    targets_sizeFactors = c(targets_sizeFactors, new_targets_sizeFactors)
  }
  if (COMBINE_CONTROLS == "FALSE") {
    # Add controls counts and size factors to tibbles
    if (!exists("controls_counts")) {
      controls_counts = select(counts, matches(control_pattern))
    } else {
      controls_counts = add_column(controls_counts, select(counts, matches(control_pattern)))
    }
    new_controls_sizeFactors = sizeFactors %>%
      filter(grepl(control_pattern, name)) %>%
      pull(value)
    names(new_controls_sizeFactors) = sizeFactors %>%
      filter(grepl(control_pattern, name)) %>%
      pull(name)
    if (!exists("controls_sizeFactors")) {
      controls_sizeFactors = new_controls_sizeFactors
    } else {
      controls_sizeFactors = c(controls_sizeFactors, new_controls_sizeFactors)
    }
  }
  
  # Save DESeq results
  result_data = select(data, contains(abbreviation))
  result_data = add_column(result_data, 
                           "baseMean" = DESeq_result$baseMean, 
                           "log2FoldChange" = DESeq_result$log2FoldChange, 
                           "pValueAdjusted" = DESeq_result$padj)
  result_data = add_column(result_data, "Ensembl_Gene" = gene, .before = 1)
  file_name = str_sub(base_name, start = str_locate(base_name, pattern = abbreviation)[1])
  write_tsv(result_data, file = paste0(data_folder, "DESeq2/", paste(file_name, deseq_prefix, sep = "_"), ".txt"))
}

if (COMBINE_CONTROLS == "TRUE") {
  controls_counts = select(counts, matches(control_pattern))
  controls_sizeFactors = sizeFactors %>%	
    filter(grepl("-C", name)) %>%	
    pull(value)	
  names(controls_sizeFactors) = sizeFactors %>%	
    filter(grepl("-C", name)) %>%	
    pull(name)
}


#### MERGE AND SAVE ALL RESULTS ####

# Merge DESeq data into one data table
cat("Merging all DESeq2 results", "\n")
for (file in list.files(path = data_folder, pattern = "[^/_]+-[^/_]*[[:digit:]]+C_DESeq2\\.txt", recursive = TRUE, full.names = TRUE)) {
  cat("  merging file ", file, "\n")
  base_name = str_sub(file, end = str_locate(file, pattern = "C_DESeq2\\.txt")[1])
  
  # Read file
  DESeq_data = read_tsv(file, show_col_types = FALSE)
  if (!exists("all_DESeq_data")) {
    all_DESeq_data = select(DESeq_data, Ensembl_Gene)
  }
  
  # Select and rename columns
  DESeq_data = select(DESeq_data, c(Ensembl_Gene, log2FoldChange, pValueAdjusted))
  name = str_extract(base_name, pattern = "[^/_]+-[^/_]*[:digit:]+C")
  colNameFC = str_c(name, "-FC")
  colNameP = str_c(name, "-P")
  DESeq_data = setNames(DESeq_data, c("Ensembl_Gene", colNameFC, colNameP))
  
  # Join data
  all_DESeq_data = inner_join(all_DESeq_data, DESeq_data, by = "Ensembl_Gene")
}
# Add common gene names
genes = read_tsv(genes_file, show_col_types = FALSE)
genes_joined = left_join(all_DESeq_data, genes, by = "Ensembl_Gene")
all_DESeq_data = add_column(all_DESeq_data, "Common_Gene" = genes_joined$Common_Gene, .after = 1)
# Save all DESeq data to file
write_tsv(all_DESeq_data, file = paste0(data_folder, paste(date, PROJECT_NAME, "DESeq2Results.txt", sep = "_")))


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
  text_label_size = 28.8 / ncol(cor_mat)
  corrplot(cor_mat, 
           method="circle", 
           order = "hclust", 
           col.lim = c(0, 1), 
           is.corr = FALSE,
           tl.cex = ifelse(text_label_size < 1, text_label_size, 1),
           tl.col = "black", 
           hclust.method = "centroid", 
           addrect = num_squares, 
           col = colors(150))}

create_and_save_cor_plot = function(cor_mat, targets_or_controls) {
  colors = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", 
                              "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  # Save as PDF
  pdf(file = paste0(figures_folder, paste(date, PROJECT_NAME, targets_or_controls, "PearsonCorrelationPlot.pdf", sep = "_")), 
      width = 7, height = 7)
  create_cor_plot(cor_mat, targets_or_controls, colors)
  invisible(dev.off())
  # Save as PNG
  png(file = paste0(figures_folder, paste(date, PROJECT_NAME, targets_or_controls, "PearsonCorrelationPlot.png", sep = "_")), 
      width = 7, height = 7, units = "in", res = 300)
  create_cor_plot(cor_mat, targets_or_controls, colors)
  invisible(dev.off())
}

# Targets
cat("Creating Pearson correlation plot for targets", "\n")
targets_cor_mat = prepare_cor_mat(targets_counts, targets_sizeFactors)
create_and_save_cor_plot(targets_cor_mat, "Targets")
# Controls
cat("Creating Pearson correlation plot for controls", "\n")
controls_cor_mat = prepare_cor_mat(controls_counts, controls_sizeFactors)
create_and_save_cor_plot(controls_cor_mat, "Controls")
# Targets and controls
cat("Creating Pearson correlation plot for targets and controls", "\n")
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
  text_size = 100 / ncol(counts_tibble)
  tsne_plot = ggplot(data = combined_tsne, aes(x = SNE_1, y = SNE_2, label = name)) + 
    geom_point(size=5, alpha = 0.3, aes(color = factor(color))) + 
    geom_text_repel(size = ifelse(text_size < 4, text_size, 4), point.padding = -10, force = 0.1, force_pull = 100, max.overlaps = 20)  +
    theme_bw() +
    coord_cartesian(clip = "off") +
    scale_color_manual(guide = "none", values = c("Target" = "red4", "Control" = "dodgerblue4")) +
    labs(x = "tSNE 2", y = "tSNE 1") +
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"), 
          plot.margin = unit(c(3,3,1,1), "lines"))
  
  # Save as PDF
  pdf(file = paste0(figures_folder, paste(date, PROJECT_NAME, targets_or_controls, "tSNEPlot.pdf", sep = "_")), 
      width = 7, height = 7)
  print(tsne_plot)
  invisible(dev.off())
  # Save as PNG
  png(file = paste0(figures_folder, paste(date, PROJECT_NAME, targets_or_controls, "tSNEPlot.png", sep = "_")), 
      width = 7, height = 7, units = "in", res = 300)
  print(tsne_plot)
  invisible(dev.off())
}

# Targets
cat("Creating tSNE plot for targets", "\n")
create_and_save_tsne_plot(targets_counts, "Targets")
# Controls
cat("Creating tSNE plot for controls", "\n")
create_and_save_tsne_plot(controls_counts, "Controls")
# Targets and controls
cat("Creating tSNE plot for targets and controls", "\n")
create_and_save_tsne_plot(targets_and_controls_counts, "TargetsAndControls")


#### CLEAR ENVIRONMENT ####
rm(list = ls())
