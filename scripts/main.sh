#!/usr/bin/env bash


set -e
unset PYTHONPATH
unset R_LIBS_SITE


# PARSE COMMAND LINE ARGUMENTS

while [[ $# -gt 0 ]]
do
  case "${1}" in
    -h|--help)
      echo ""
      echo "ProcessCounts"
      echo ""
      echo "Usage:"
      echo "  ProcessCounts/scripts/main.sh --option <argument>"
      echo ""
      echo "Options:"
      echo "  -h, --help                show this help message"
      echo ""
      echo "  -d, --data-dir DIR        path to directory containing processed data files"
      echo "                            (no default)"
      echo ""
      echo "  -o, --output-dir DIR      path to directory for output tables and figures"
      echo "                            (default is your current directory)"
      echo ""
      echo "  -n, --name NAME           name of project to be added to file names of tables and figures"
      echo "                            (default is ProjectName)"
      echo ""
      echo "  -e, --env-dir DIR         path to directory with conda virtual environment"
      echo "                            (default is /storage/fazal/software/R4.1.2/venv)"
      echo ""
      echo "  -t, --temp-dir DIR        path to directory for temporary files"
      echo "                            (default is /storage/fazal/tmp)"
      echo ""
      exit 0
      ;;
    -d|--data-dir)
      DATADIR="${2}"
      shift
      shift
      ;;
    -o|--output-dir)
      OUTPUTDIR="${2}"
      shift
      shift
      ;;
    -n|--name)
      PROJECTNAME="${2}"
      shift
      shift
      ;;
    -e|--env-dir)
      ENVDIR="${2}"
      shift
      shift
      ;;
    -t|--temp-dir)
      TEMPDIR="${2}"
      shift
      shift
      ;;
    -*|--*)
      echo "Unknown option: ${1}"
      exit 1
      ;;
    *)
      echo "Unknown argument: ${1}"
      exit 1
      ;;
  esac
done


# SET UNSPECIFIED ARGUMENTS TO DEFAULTS

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [[ -z "${DATADIR}" ]]
then
  echo "  Processed data directory not specified. Please supply path to processed data directory."
  exit 1
else
  echo "  Processed data directory set to ${DATADIR}."
fi

if [[ -z "${OUTPUTDIR}" ]]
then
  OUTPUTDIR="${PWD}"
  echo "  Output directory not specified. Defaulting to ${OUTPUTDIR}."
else
  echo "  Output directory set to ${OUTPUTDIR}."
fi

if [[ -z "${PROJECTNAME}" ]]
then
  PROJECTNAME=ProjectName
  echo "  Project name not specified. Defaulting to ${PROJECTNAME}."
else
  echo "  Project name set to ${PROJECTNAME}."
fi

if [[ -z "${ENVDIR}" ]]
then
  ENVDIR=/storage/fazal/software/R4.1.2/venv
  echo "  Environment directory not specified. Defaulting to ${ENVDIR}."
else
  echo "  Environment directory set to ${ENVDIR}."
fi
# Activate environment
export PATH="${ENVDIR}"/bin:"${PATH}"

if [[ -z "${TEMPDIR}" ]]
then
  TEMPDIR=/storage/fazal/tmp
  echo "  Temporary directory not specified. Defaulting to ${TEMPDIR}."
else
  echo "  Temporary directory set to ${TEMPDIR}."
fi
# Avoid dumping too much temporary data into system tmp directory
export TMPDIR="${TEMPDIR}"
export TEMP="${TEMPDIR}"
export TMP="${TEMPDIR}"


# CALL PROCESSCOUNTS

while true
do
  read -p "Do you wish to start the ProcessCounts pipeline? [y/n] " yn
  case $yn in
    [Yy]*)
      sbatch "${SCRIPTDIR}"/ProcessCounts.sh "${SCRIPTDIR}" "${DATADIR}" "${OUTPUTDIR}" "${PROJECTNAME}"
      break
      ;;
    [Nn]*)
      exit 0
      ;;
    *)
      echo "  Please answer yes or no."
      ;;
  esac
done

