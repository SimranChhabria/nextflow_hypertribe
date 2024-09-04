#!/bin/bash
#BSUB -J 3_htseq_count[1-24]
#BSUB -W 8:00
#BSUB -R rusage[mem=20]
#BSUB -n 8
#BSUB -C 0
#BSUB -o logs/3_htseq_count_%J_%I.stdout
#BSUB -eo logs/3_htseq_count_%J_%I.stderr

# Set parameters
lab=ganeshk
project=Project_15162
genome=hg38_wtzfpadar_fszfpadar

# Initialize conda environment
source ~/.bashrc
conda activate genome
cd $LS_SUBCWD

# Get sample
samples=(../../projects/${lab}/${project}/output_data/deseq2/1_star_align/*/)
i=$((LSB_JOBINDEX-1))
sample=$(basename "${samples[$i]}")
echo $sample

# Get directories
align_folder=../../projects/${lab}/${project}/output_data/deseq2/1_star_align/${sample}/
gtf_path=../../projects/${lab}/${project}/genome_data/${genome}.gtf
counts_folder=../../projects/${lab}/${project}/output_data/deseq2/3_counts/

# Make directories
mkdir -p logs ${counts_folder}

# Htseq: count genes in reads
echo "Htseq: count genes in reads"
htseq-count \
-f bam \
-r pos \
-s reverse \
-a 10 \
-i gene_id \
-m union \
${align_folder}/${sample}.bam \
${gtf_path} > \
${counts_folder}/${sample}.counts \
