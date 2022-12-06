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
MIN_TRANSCRIPT_LENGTH="${6}"


# RUN PROCESSCOUNTS PIPELINE
echo "________________________________________"
echo "STEP 1: QUALITY CONTROL"
Rscript --vanilla "${SCRIPTDIR}"/1_QualityControl.R "${SCRIPTDIR}" "${DATADIR}" "${OUTPUTDIR}" "${PROJECT_NAME}" "${COMBINE_CONTROLS}" "${MIN_TRANSCRIPT_LENGTH}"
echo "________________________________________"
echo "STEP 2: DESEQ2 ANALYSIS"
Rscript --vanilla "${SCRIPTDIR}"/2_DESeq2Analysis.R "${SCRIPTDIR}" "${OUTPUTDIR}" "${PROJECT_NAME}" "${COMBINE_CONTROLS}" "${MIN_TRANSCRIPT_LENGTH}"
echo "________________________________________"
echo "FINAL STEP: RENDER ANALYSIS SUMMARY"
Rscript --vanilla "${SCRIPTDIR}"/5_RenderAnalysisSummary.R "${SCRIPTDIR}" "${OUTPUTDIR}" "${PROJECT_NAME}"
echo "________________________________________"
echo "DONE"
echo ""

