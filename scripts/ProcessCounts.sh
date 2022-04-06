#!/usr/bin/env bash

#SBATCH -n 1                        # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-12:00                  # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=16G                   # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name=ProcessCounts    # Short name for the job
#SBATCH --output=%j.ProcessCounts.log


set -e
unset PYTHONPATH
unset R_LIBS_SITE

# Avoid dumping too much temporary data into system tmp directory
export TMPDIR=/storage/fazal/tmp
export TEMP=/storage/fazal/tmp
export TMP=/storage/fazal/tmp

# Activate environment
export PATH=/storage/fazal/software/R4.1.2/venv/bin:"${PATH}"

# Define variables for project directories
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DATADIR="${1}"
OUTPUTDIR="${PWD}"

# Get project file prefix from command arguments
filePrefix="${2}"


# RUN PROCESSCOUNTS PIPELINE
echo "________________________________________"
echo "STEP 1: QUALITY CONTROL"
Rscript --vanilla "${SCRIPTDIR}"/1_QualityControl.R "${DATADIR}" "${SCRIPTDIR}" "${filePrefix}"
echo "________________________________________"
echo "STEP 2: DESEQ ANALYSIS"
Rscript --vanilla "${SCRIPTDIR}"/2_DESeqAnalysis.R "${SCRIPTDIR}" "${OUTPUTDIR}" "${filePrefix}"
echo "________________________________________"
echo "FINAL STEP: RENDER ANALYSIS SUMMARY"
Rscript --vanilla "${SCRIPTDIR}"/5_RenderAnalysisSummary.R "${SCRIPTDIR}" "${OUTPUTDIR}" "${filePrefix}"
echo "________________________________________"
echo "DONE"
echo ""

