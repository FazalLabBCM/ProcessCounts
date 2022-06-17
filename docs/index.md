# Data Analysis for RNA-Seq Experiments

The ProcessCounts pipeline is designed to process your FASTQ read counts with 
[DESeq2](https://github.com/mikelove/DESeq2) to obtain enrichment values and an analysis 
summary full of quality-control figures.


## Setup

If you are not a member of the Fazal Lab and don't have access to Baylor College 
of Medicine's MHGCP cluster, follow this link to 
[download and setup the pipeline](https://fazallabbcm.github.io/ProcessCounts/DownloadAndSetup) 
on your local computing environment.


## Running the Pipeline

1. From the command line, add the ProcessCounts scripts folder to your PATH environment variable. 
   For Fazal Lab members, this can be done with the following code:
   ```
   export PATH=/storage/fazal/pipelines/ProcessCounts/scripts:"${PATH}"
   ```
   
2. Run the following code:
   ```
   ProcessCounts -d /path/to/data -o /path/to/output -n ProjectName
   ```
   > **Note:**
   > Only run ProcessCounts once. This pipeline analyzes the data from every condition in your 
   > experiment at the same time. If your aligned data files are separated into subfolders, 
   > just supply the path to the folder containing the subfolders.

3. Make sure that the processed data and output file paths are correct. Then enter "y" to 
   start the pipeline.


## What next?

You can check the log file in your output folder to see the progress of your job as it runs.

If you encounter any errors, see our 
[troubleshooting](https://fazallabbcm.github.io/ProcessCounts/Troubleshooting) page for help.

Once [ProcessCounts](https://fazallabbcm.github.io/ProcessCounts) has finished...
