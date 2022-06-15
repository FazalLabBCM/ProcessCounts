#!/usr/bin/env bash

#SBATCH -n 1                        # Number of cores
#SBATCH -N 1                        # Number of nodes
#SBATCH -t 0-08:00                  # Runtime in D-HH:MM
#SBATCH --mem=8G                    # Memory pool for all cores

set -e

# Define variables for project directories
SCRIPTDIR="${1}"
DATADIR="${2}"
OUTPUTDIR="${3}"
PROJECT_NAME="${4}"
COMBINE_CONTROLS="${5}"


# RUN PROCESSCOUNTS PIPELINE
echo "________________________________________"
echo "STEP 1: QUALITY CONTROL"
Rscript --vanilla "${SCRIPTDIR}"/1_QualityControl.R "${SCRIPTDIR}" "${DATADIR}" "${PROJECT_NAME}" "${COMBINE_CONTROLS}"
echo "________________________________________"
echo "STEP 2: DESEQ ANALYSIS"
Rscript --vanilla "${SCRIPTDIR}"/2_DESeqAnalysis.R "${SCRIPTDIR}" "${OUTPUTDIR}" "${PROJECT_NAME}" "${COMBINE_CONTROLS}"
echo "________________________________________"
echo "FINAL STEP: RENDER ANALYSIS SUMMARY"
Rscript --vanilla "${SCRIPTDIR}"/5_RenderAnalysisSummary.R "${SCRIPTDIR}" "${OUTPUTDIR}" "${PROJECT_NAME}"
echo "________________________________________"
echo "DONE"
echo ""

