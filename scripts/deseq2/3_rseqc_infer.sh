#!/bin/bash
#BSUB -J 3_htseq_count
#BSUB -W 8:00
#BSUB -R rusage[mem=20]
#BSUB -n 8
#BSUB -C 0
#BSUB -o logs/3_rseqc_infer_%J.stdout
#BSUB -eo logs/3_rseqc_infer_%J.stderr

# Set parameters
lab=ganeshk
project=Project_15162
genome_base=hg38

# Initialize conda environment
source ~/.bashrc
conda activate genome
cd $LS_SUBCWD

# Get sample
samples=(../../projects/${lab}/${project}/output_data/deseq2/1_star_align/*/)
sample=$(basename "${samples[0]}")
echo $sample

# Get directories
align_folder=../../projects/${lab}/${project}/output_data/deseq2/1_star_align/${sample}/
bed_path=../../projects/${lab}/${project}/genome_data/${genome_base}_RefSeq.bed

# RSeQC: infer strandedness
echo "RSeQC: infer strandedness"
infer_experiment.py \
-i ${align_folder}/${sample}.bam \
-r ${bed_path}
