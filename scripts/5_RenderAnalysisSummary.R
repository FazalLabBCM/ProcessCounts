#!/usr/bin/env Rscript


## PREPARE ENVIRONMENT ##

# Clear environment
rm(list = ls())

# Retrieve command line arguments
args = commandArgs(trailingOnly = TRUE)
if (interactive()) {
  projectdir = "."
  outputdir = "."
  file_prefix = "test"
} else {
  projectdir = args[1]
  outputdir = args[2]
  file_prefix = args[3]
}

# Import libraries
library(rmarkdown)

# Date
date = format(Sys.Date(), "%Y%m%d")


## RENDER OUTPUT SUMMARY ##
render_script_path = paste(projectdir, "AdditionalScripts", "AnalysisSummary.Rmd", sep = "/")
output_file_path = paste(outputdir, paste(date, file_prefix, "AnalysisSummary.html", sep = "_"), sep = "/")
render(render_script_path, 
       params = list(projectdir = projectdir, outputdir = outputdir), 
       output_file = output_file_path, 
       intermediates_dir = outputdir)
