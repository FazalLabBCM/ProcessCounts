#!/usr/bin/env Rscript


## PREPARE ENVIRONMENT ##

# Clear environment
rm(list = ls())

# Retrieve command line arguments
args = commandArgs(trailingOnly = TRUE)
if (interactive()) {
  SCRIPTDIR = "."
  OUTPUTDIR = "."
  PROJECT_NAME = "test"
} else {
  SCRIPTDIR = args[1]
  OUTPUTDIR = args[2]
  PROJECT_NAME = args[3]
}

# Import libraries
suppressPackageStartupMessages(library(rmarkdown))

# Date
date = format(Sys.Date(), "%Y%m%d")


## RENDER OUTPUT SUMMARY ##
render_script_path = paste(SCRIPTDIR, "AdditionalScripts", "AnalysisSummary.Rmd", sep = "/")
output_file_path = paste(OUTPUTDIR, paste(date, PROJECT_NAME, "AnalysisSummary.html", sep = "_"), sep = "/")
suppressWarnings(render(render_script_path, 
                        params = list(SCRIPTDIR = SCRIPTDIR, OUTPUTDIR = OUTPUTDIR), 
                        output_file = output_file_path, 
                        intermediates_dir = OUTPUTDIR))
