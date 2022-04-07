#!/usr/bin/env bash

#SBATCH -n 1                        # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-12:00                  # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=16G                   # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name=ProcessCounts    # Short name for the job
#SBATCH --output=%j.ProcessCounts.log


# Define variables for project directories
SCRIPTDIR="${1}"
DATADIR="${2}"
OUTPUTDIR="${3}"
PROJECTNAME="${4}"


# RUN PROCESSCOUNTS PIPELINE
echo "________________________________________"
echo "STEP 1: QUALITY CONTROL"
Rscript --vanilla "${SCRIPTDIR}"/1_QualityControl.R "${DATADIR}" "${SCRIPTDIR}" "${PROJECTNAME}"
echo "________________________________________"
echo "STEP 2: DESEQ ANALYSIS"
Rscript --vanilla "${SCRIPTDIR}"/2_DESeqAnalysis.R "${SCRIPTDIR}" "${OUTPUTDIR}" "${PROJECTNAME}"
echo "________________________________________"
echo "FINAL STEP: RENDER ANALYSIS SUMMARY"
Rscript --vanilla "${SCRIPTDIR}"/5_RenderAnalysisSummary.R "${SCRIPTDIR}" "${OUTPUTDIR}" "${PROJECTNAME}"
echo "________________________________________"
echo "DONE"
echo ""

