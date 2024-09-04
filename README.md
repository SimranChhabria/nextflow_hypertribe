# Hypertribe nextflow pipeline:

This repo consists of the nextflow pipeline for HYPERTRIBE and deseq2 analysis in association with the nature paper.
https://www.nature.com/articles/s41467-020-15814-8#Sec8

## STEPS: 

- Pre-processing:
    - MutliQC
    - Qualimap
- Alignment:
    - Mapping the reads with reference genome

#### HYPERTRIBE STEPS:

- Post-Alignment Filtering:
    - SplitNCigarReads and reassign mapping qualities
    - Indel Realignment and Base Recalibration
    - Germline Variant calling and filtering
- Differential Analysis:
    - Base editing stats and comparisons
- Plots
 
#### RNA-seq expression STEPS:

- RNA seq counts from merged bams
- Differential gene expression using DESeq2
- Plots

## MAMBA ENVIRONMENT:
- Create a new mamba environment and installing tools needed for the pipeline
```
mamba env create -n nextflow_mamba -f env.yaml
```
## HOW TO RUN THE PIPELINE:
- INPUTS:
    - input_reads   =</path/of/the/directory/containing/fastq/files>
    - output_dir    =</path/of/the/directory/to/store/results>
    - genome_dir    =</path/of/the/directory/containing/reference/indexes>
    - sample_sheet  =</path/of/the/sample_sheet.csv>
          
     **Example of the sample_sheet.csv** 

     ```
     Condition,Control
     condition_samples_group,ctrl_sample_group
     MSI2-ADAR,MIG
     ```
       

- PARAMETERS:
    - diff_threshold_deseq2    = Value of threshold for differential analysis for RNA seq
    - diff_threshold_hyper     = Value of threshold for differential analysis for base edit
    - padj_threshold           = Value of padj threshold
    - workflow_HY              = "HYPERTRIBE" or "NULL"
    - workflow_DE              = "DESEQ2" or "NULL"
    - genome_build             = "hg38"

- COMMAND:

  Change all the required parameters in the bash script and run the following command

  ```
  bash run_hypertribe_pipeline.sh
  ```





