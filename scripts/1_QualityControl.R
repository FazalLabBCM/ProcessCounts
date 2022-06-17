#!/usr/bin/env Rscript


#### PREPARE ENVIRONMENT ####

# Clear environment
rm(list = ls())

# Retrieve command line arguments
args = commandArgs(trailingOnly = TRUE)
if (interactive()) {
  SCRIPTDIR = "."
  DATADIR = "./data"
  OUTPUTDIR = "."
  PROJECT_NAME = "test"
  COMBINE_CONTROLS = "FALSE"
} else {
  SCRIPTDIR = args[1]
  DATADIR = args[2]
  OUTPUTDIR = args[3]
  PROJECT_NAME = args[4]
  COMBINE_CONTROLS = args[5]
}

# Import libraries
library(tidyverse)
library(colorspace)
library(grid)


#### DEFINE VARIABLES ####

# Create folders for data and figures
data_folder = paste(OUTPUTDIR, "DataFiles", sep = "/")
figures_folder = paste(OUTPUTDIR, "Figures", sep = "/")
dir.create(data_folder)
dir.create(figures_folder)

# Folders
data_folder = paste0(data_folder, "/")
figures_folder = paste0(figures_folder, "/")

# Files
genes_file = paste(SCRIPTDIR, "AdditionalFiles", "CommonGeneNames.txt", sep = "/")
lengths_file = paste(SCRIPTDIR, "AdditionalFiles", "TranscriptLengths.txt", sep = "/")

# Date
date = format(Sys.Date(), "%Y%m%d")

# Store sample location abbreviations for later use
locations_conditions = c()


#### COLLECT DATA FOR EACH SAMPLE ####

for (file in list.files(path = DATADIR, pattern = ".*_aligned\\.txt", recursive = TRUE, full.names = TRUE)) {
  cat("Working with file", file, "\n")
  base_name = str_sub(file, end = str_locate(file, pattern = "_aligned\\.txt")[1])
  file_name = str_sub(base_name, start = str_locate(file, pattern = "/[^/_]+_[^/_]+_[:alpha:]+_[:digit:]_aligned\\.txt")[1])
  
  # Extract sample name (location + experimental condition + target/control + number)
  # Extract location
  location = str_extract(file_name, pattern = "/[^/_]+_")
  location = str_sub(location, 2, -2)
  # Extract experimental condition
  condition_start = str_locate(file_name, pattern = location)[2] + 2
  condition = str_sub(file_name, start = condition_start, end = -1)
  condition = str_sub(condition, end = str_locate(condition, pattern = "^[^/_]+")[2])
  # Join location and condition
  if (condition == "none") {
    location_condition = location
  } else {
    location_condition = str_c(location, "-", condition)
  }
  if (!(location_condition %in% locations_conditions)) {
    locations_conditions = c(locations_conditions, location_condition)
  }
  # Extract target/control
  if (str_detect(base_name, pattern = "_[Tt]arget_")) {
    type = "T"
  } else {
    type = "C"
  }
  # Extract number
  number = str_extract(base_name, pattern = "_[:digit:]_")
  number = str_sub(number, 2, -2)
  # Join sample name
  sample_name = str_c(location_condition, "-", type, number)
  cat("  sample name:", sample_name, "\n")

  # Read counts file
  data = read_tsv(file, col_names = c("Ensembl_Gene", sample_name), show_col_types = FALSE)
  data = dplyr::slice(data, 1:(length(data$Ensembl_Gene) - 5))
  # Join data
  if (!exists("all_data")) {
    all_data = data
  } else {
    all_data = inner_join(all_data, data, by = "Ensembl_Gene")
  }
  
  # Read STAR final log out file
  log_file = str_c(base_name, "trim_starLog.final.out")
  log = read_tsv(log_file, col_names = c("type_of_info", "value"), show_col_types = FALSE)
  # Filter/select important information
  sample = rep(sample_name, 5)
  if (!exists("reads")) {
    reads_vector = c("Reads mapped uniquely", 
                     "Reads mapped to multiple loci", 
                     "Reads mapped to too many loci", 
                     "Reads unmapped: too short", 
                     "Reads unmapped: other")
    reads = factor(reads_vector, levels = rev(reads_vector))
  }
  number = c(log$value[8], 
             log$value[23], 
             log$value[25], 
             log$value[30], 
             log$value[32])
  percent = c(str_sub(log$value[9], end = -2),
              str_sub(log$value[24], end = -2), 
              str_sub(log$value[26], end = -2), 
              str_sub(log$value[31], end = -2), 
              str_sub(log$value[33], end = -2))
  log_final = tibble("sample" = sample, "reads" = reads, "number" = number, "percent" = percent)
  log_final = mutate(log_final, number = as.double(number), percent = as.double(percent))
  # Join data
  if (!exists("all_log_final")) {
    all_log_final = log_final
  } else {
    all_log_final = add_row(all_log_final, 
                            sample = log_final$sample,
                            reads = log_final$reads,
                            number = log_final$number,
                            percent = log_final$percent)
  }
}
# Order all counts data columns alphabetically
all_data = select(all_data, "Ensembl_Gene", order(colnames(all_data)[2:length(colnames(all_data))]) + 1)
# Add common gene names
genes = read_tsv(genes_file, show_col_types = FALSE)
genes_joined = left_join(all_data, genes, by = "Ensembl_Gene")
all_data = add_column(all_data, "Common_Gene" = genes_joined$Common_Gene, .after = 1)
# Calculate total number of reads per sample
all_log_totals = all_log_final %>%
  group_by(sample) %>%
  summarize(total = sum(number))


#### SAVE DATA TABLES and BAR GRAPH ####

# Save all counts data and log data to files
cat("Saving data files", "\n")
write_tsv(all_data, file = paste0(data_folder, paste(date, PROJECT_NAME, "HTSeqCountsData.txt", sep = "_")))
write_tsv(all_log_final, file = paste0(data_folder, paste(date, PROJECT_NAME, "STARLogData.txt", sep = "_")))
write_tsv(all_log_totals, file = paste0(data_folder, paste(date, PROJECT_NAME, "STARTotalReadsData.txt", sep = "_")))

# Plot percentage reads mapped bar graph
cat("Making percentage reads mapped bar graph", "\n")
percentage_reads_mapped_bar_graph = ggplot(all_log_final, aes(x = sample, y = percent, fill = reads)) +
  geom_bar(stat = "identity") +
  labs(x = "Samples", y = "Percentage (%)", fill = "Reads") +
  scale_fill_discrete_diverging(palette = "Blue-Red", nmax = 9, order = c(9, 8, 7, 6, 2)) +
  coord_flip(clip = "off") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", 
        legend.direction = "vertical", 
        plot.margin = unit(c(3,5,1,1), "lines"))
# Add total number of reads per sample to bar graph
percentage_reads_mapped_bar_graph = percentage_reads_mapped_bar_graph +
  annotation_custom(grob = textGrob("  Reads\n(Millions)", hjust = 0, vjust = 0, gp = gpar(col = "grey30", fontsize = 10, lineheight = 0.9)), 
                    ymin = 110, ymax = 110, xmin = nrow(all_log_totals) + 0.5, xmax = nrow(all_log_totals) + 1)
for (i in 1:length(all_log_totals$total)) {
  percentage_reads_mapped_bar_graph = percentage_reads_mapped_bar_graph +
    annotation_custom(grob = textGrob(label = round(all_log_totals$total[i], -6)/1000000, hjust = 0),
                      ymin = 110, ymax = 110, xmin = i, xmax = i)
}

# Save as PDF
pdf(width = 6, height = 10, file = paste0(figures_folder, paste(date, PROJECT_NAME, "PercentageReadsMappedBarGraph.pdf", sep = "_")))
print(percentage_reads_mapped_bar_graph)
dev.off()
# Save as PNG
png(width = 6, height = 10, units = "in", res = 300, file = paste0(figures_folder, paste(date, PROJECT_NAME, "PercentageReadsMappedBarGraph.png", sep = "_")))
print(percentage_reads_mapped_bar_graph)
dev.off()


#### PREPARE DATA FOR DESEQ ANALYSIS ####

cat("Preparing data for DESeq analysis", "\n")

# Separate targets and controls
targets = select(all_data, matches(".+-T[[:digit:]]$"))
controls = select(all_data, matches(".+-C[[:digit:]]$"))

# Create new folder for DESeq data files
deseq_folder = paste0(data_folder, "DESeq")
dir.create(deseq_folder)
deseq_folder = paste0(deseq_folder, "/")

for (i in 1:length(locations_conditions)) {
  location_condition = locations_conditions[i]
  
  # Join location-condition targets with location-condition controls
  location_condition_targets = select(targets, matches(str_c(location_condition, "-T[[:digit:]]$")))
  if (COMBINE_CONTROLS == "TRUE") {
    location_condition_controls = controls
  } else {
    location_condition_controls = select(controls, matches(str_c(location_condition, "-C[[:digit:]]$")))
  }
  DESeq_data = add_column(location_condition_controls, location_condition_targets, .before = 1)
  DESeq_data = add_column(DESeq_data, "Ensembl_Gene" = all_data$Ensembl_Gene, .before = 1)
  
  # Save file
  file_name = str_c(location_condition, "-", ncol(location_condition_controls), "C")
  write_tsv(DESeq_data, file = paste0(deseq_folder, file_name, ".txt"))
}


#### CREATE CORRELATION PLOTS ####

generate_correlation_plots = function(data, color) {
  for (location_condition in locations_conditions) {
    location_condition_data = select(data, contains(location_condition))
    colnames = colnames(location_condition_data)
    if (length(colnames) < 2) {
      cat("  not enough replicates to create correlation plot for condition", location_condition, "\n")
      next
    }
    for (i in 1:(length(colnames) - 1)) {
      for (j in (i + 1):length(colnames)) {
        col_x = colnames[i]
        col_y = colnames[j]
        # Get Spearman correlation value
        x = pull(location_condition_data, !!sym(col_x))
        y = pull(location_condition_data, !!sym(col_y))
        r = round(cor(x, y, method = "spearman"), 3)
        # Create correlation plot
        correlation_plot = ggplot(data, aes(x = !!sym(col_x), y = !!sym(col_y))) +
          geom_point(color = color, alpha = 1/10) +
          scale_x_log10(labels = function(x) format(x, scientific = FALSE)) +
          scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
          annotate(geom="text", x=10, y=10000, label=bquote(italic(r) == .(r))) +
          theme_bw()
        # Save as PDF
        pdf(file = paste0(figures_folder, paste(date, PROJECT_NAME, col_x, col_y, "CorrelationPlot.pdf", sep = "_")))
        print(correlation_plot)
        dev.off()
        # Save as PNG
        png(width = 5, height = 5, units = "in", res = 300, file = paste0(figures_folder, paste(date, PROJECT_NAME, col_x, col_y, "CorrelationPlot.png", sep = "_")))
        print(correlation_plot)
        dev.off()
      }
    }
  }
}
cat("Creating correlation plots for targets", "\n")
generate_correlation_plots(targets, "red4")
cat("Creating correlation plots for controls", "\n")
generate_correlation_plots(controls, "dodgerblue4")


#### CLEAR ENVIRONMENT ####
rm(list = ls())
