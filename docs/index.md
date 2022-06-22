# Data Analysis for RNA-Seq Experiments

The ProcessCounts pipeline is designed to process FASTQ read counts with 
[DESeq2](https://github.com/mikelove/DESeq2) to obtain enrichment values and an analysis 
summary full of quality-control figures.


## Setup

If you are not a member of the Fazal Lab and don't have access to Baylor College 
of Medicine's MHGCP cluster, follow this link to 
[download and setup the pipeline](https://fazallabbcm.github.io/ProcessCounts/DownloadAndSetup) 
on your local computing environment.

Make sure you run the [TrimMapCount](https://fazallabbcm.github.io/TrimMapCount) pipeline first.
Once it has finished successfully, you are ready to run the ProcessCounts pipeline.


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
   > **Note:** Only run ProcessCounts once. This pipeline analyzes the data from every condition 
   > in your experiment at the same time. If your aligned data files are separated into subfolders, 
   > just supply the path to the folder containing the subfolders.

3. Make sure that the processed data and output file paths are correct. Then enter "y" to 
   start the pipeline.


## What next?

You can check the log file in your output folder to see the progress of your job as it runs.

If you encounter any errors, see our 
[troubleshooting](https://fazallabbcm.github.io/ProcessCounts/Troubleshooting) page for help.

Once ProcessCounts has finished, open the Analysis Summary HTML file and look for the following things:

  1. Percentage of reads mapped uniquely is greater than 75% for each sample (remove any 
     samples that do not meet this criterion from future analyses) 

  2. Correlation between sample replicates is greater than 0.75 for the best correlated 
     replicates of each condition (repeat the experiment if the two best correlated replicates 
     for any condition have a correlation coefficient less than 0.7) 

  3. Examine the other plots, looking for correlation between replicates as well as differences 
     between experimental conditions

> Additional Analyses: If you wish to compare your results to the results of another experiment, 
> make a volcano plot or an ROC curve. (Code for a generic volcano plot and ROC curve can be 
> found [here]().)

Finally, if you haven't viewed the genome tracks for your experiment yet, use the 
[BamToBigWig](https://fazallabbcm.github.io/BamToBigWig) pipeline to visualize your data.
