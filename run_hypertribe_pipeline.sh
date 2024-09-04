#!/bin/bash
#BSUB -J Nextflow_Hypertribe
#BSUB -W 48:00
#BSUB -n 1
#BSUB -C 0
#BSUB -o nextflow_hypertribe_%J.stdout
#BSUB -eo nextflow_hypertribe_%J.stderr

source ~/.bashrc
module load java/11.0.12
mamba activate nextflow_mamba
# module load samtools/1.13
# module load gatk/4.3.0.0
#module load R/R-4.2.0

genome_build="hg38"
genome_dir="/data/morrisq/baalii/nextflow_pipeline/genome/${genome_build}"
padj_threshold="0.05"
diff_threshold_deseq2="1"
diff_threshold_hyper="0.1"
project_name="Project_15619"
workflow_HY="NULL"
workflow_DE="DESEQ2"


genome="${genome_dir}/${genome_build}.fa"
genome_index="${genome_dir}/${genome_build}.fa.fai"
genome_dict="${genome_dir}/${genome_build}.dict"
genomeDir="${genome_dir}/${genome_build}_star_index/"
variants=${genome_dir}/"${genome_build}.vcf"
variants_index="${genome_dir}/${genome_build}.vcf.idx"
gtf="${genome_dir}/${genome_build}.gtf"
editor="ADAR"


# Set the input and output directory paths
input_reads="/data/morrisq/baalii/HyperTribe_MSI2_Wei/input_data/*15619*/*_L00{1,2,3,4}_R{1,2}_001.fastq.gz"
output_dir="/data//morrisq/simranch/nextflow/output_data/Project_15619/"
sample_sheet="/data//morrisq/simranch/nextflow/Test_Hypertribe/sample_sheet_test1.csv"



if [ ! -f "$genome" ]; then
    echo "Genome file does not exist: $genome"
    exit 1
fi

if [ ! -f "$genome_index" ]; then
    echo "Genome index file does not exist: $genome_index"
    exit 1
fi

if [ ! -f "$genome_dict" ]; then
    echo "Genome dictionary file does not exist: $genome_dict"
    exit 1
fi

if [ ! -d "$genomeDir" ]; then
    echo "Genome directory does not exist: $genomeDir"
    exit 1
fi

if [ ! -f "$variants" ]; then
    echo "Variants file does not exist: $variants"
    exit 1
fi

if [ ! -f "$variants_index" ]; then
    echo "Variants index file does not exist: $variants_index"
    exit 1
fi

if [ ! -f "$gtf" ]; then
    echo "GTF file does not exist: $gtf"
    exit 1
fi

mkdir -p $output_dir

# Run the Nextflow workflow in the background
# nextflow run prepare_genome.nf --genome_dir $genome_dir --genome $genome --variants $variants --denylist $denylist -resume -profile lsf &

nextflow run main_hypertribe.nf \
--genome_dir $genome_dir \
--genome $genome \
--genome_index $genome_index \
--genome_dict $genome_dict \
--genomeDir $genomeDir \
--variants $variants \
--variants_index $variants_index \
--gtf $gtf \
--genome_build $genome_build \
--editor $editor \
--input_dir $input_reads \
--output_dir $output_dir \
--sample_sheet $sample_sheet \
--padj_threshold $padj_threshold \
--diff_threshold_deseq2 $diff_threshold_deseq2 \
--diff_threshold_hyper $diff_threshold_hyper \
--project_name $project_name \
--workflow_HY $workflow_HY \
--workflow_DE $workflow_DE \
-profile lsf \
-resume

 