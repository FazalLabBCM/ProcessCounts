# Download and Setup

This page is for users who do not have access to the MHGCP cluster at Baylor College of Medicine.
These instructions will help you download the necessary files and set up your computing environment 
in preparation for using the ProcessCounts pipeline.


### 1. Create virtual environment

[Install Conda or Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), 
make a folder for virtual environments, and add a folder inside called "DataAnalysis". Then use the code 
below to create a virtual environment for the ProcessCounts pipeline: 

```
conda create \
  -p /path/to/environments/DataAnalysis/venv \
  r-base=4.1.2 \
  r-essentials \
  bioconductor-deseq2 \
  r-colorspace \
  r-gridgraphics \
  r-rtsne \
  r-corrplot \
  -c conda-forge \
  -c bioconda \
  -c r
```


### 2. Download pipeline repository

Visit https://github.com/FazalLabBCM/ProcessCounts, and use the green ‘Code’ button to download the repository. 


### 3. Configure pipeline defaults

Edit the CONFIG file in the repository folder to include the absolute paths to your folders:

* ENVDIR is the virtual environment that you made in the first step: `path/to/environments/DataAnalysis/venv`.
  
* TEMPDIR is a folder for temporary files.
