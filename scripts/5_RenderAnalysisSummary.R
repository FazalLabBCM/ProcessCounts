#!/usr/bin/env Rscript


## PREPARE ENVIRONMENT ##

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
library(rmarkdown)

# Date
date = format(Sys.Date(), "%Y%m%d")


## RENDER OUTPUT SUMMARY ##
render_script_path = paste(scriptdir, "AdditionalScripts", "AnalysisSummary.Rmd", sep = "/")
output_file_path = paste(outputdir, paste(date, project_name, "AnalysisSummary.html", sep = "_"), sep = "/")
render(render_script_path, 
       params = list(scriptdir = scriptdir, outputdir = outputdir), 
       output_file = output_file_path, 
       intermediates_dir = outputdir)
