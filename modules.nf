
/*
 * Process 1A: Create a FASTA genome index (.fai) with samtools for GATK
 */

process PREPARE_GENOME_SAMTOOLS { 
  tag "$genome.baseName"
 
  input: 
    path genome_dir
    path genome
 
  output: 
    path "${genome}.fai"
  
  publishDir "${genome_dir}/", mode: 'copy'

  script:
  """
  samtools faidx ${genome}
  """
}


/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process PREPARE_GENOME_PICARD {
  tag "$genome.baseName"
  label 'mem_xlarge'
  
  

  input:
    path genome_dir
    path genome
  output:
    path "${genome.baseName}.dict"

  publishDir "${genome_dir}/", mode: 'copy'

  script:
  """
  gatk CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
  """
}


/*
 * Process 1C: Create STAR genome index file.
 */

process PREPARE_STAR_GENOME_INDEX {
  tag "$genome.baseName"

  

  input:
    path genome_dir
    path genome
  output:
    path "${genome.baseName}_index"

  publishDir "${genome_dir}/", mode: 'copy'
  script:
  """  
  mkdir ${genome.baseName}_index

  STAR --runMode genomeGenerate \
       --genomeDir ${genome.baseName}_index \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """
}

process PREPARE_VCF_FILE {
  tag "$variantsFile.baseName"

  input: 
    path genome_dir
    path variantsFile
    path denylisted

  output:
    tuple \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")
  
  publishDir "${genome_dir}/", mode: 'copy'
  
  script:  
  """
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${denylisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
  """
}




 /*
 * Process 2: Align RNA-Seq reads to the genome with STAR
 */

process RNASEQ_MAPPING_STAR {
  cpus = 8
  memory = 20.GB
  time = '48h'

  tag "$sample"

  input:
    path genome
    path genomeDir
    tuple val(sample), val(lane), path(reads)

  output:
    tuple \
      val(sample), \
      val(lane), \
      path("${lane}_Aligned.sortedByCoord.uniq.bam"), \
      path("${lane}_Aligned.sortedByCoord.uniq.bam.bai"), \
      path("qualimap_${sample}_${lane}")

  script:
  """
  # ngs-nf-dev Align reads to genome
  STAR --genomeDir $genomeDir \
       --readFilesIn $reads \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999

  # Run 2-pass mapping (improve alignmets using table of splice junctions and create a new index)
  STAR --genomeDir $genomeDir \
       --readFilesIn $reads \
       --runThreadN $task.cpus \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --sjdbFileChrStartEnd SJ.out.tab \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:${sample}_${lane} LB:library PL:illumina PU:machine SM:GM12878

  # Select only unique alignments, no multimaps
  (samtools view -H Aligned.sortedByCoord.out.bam; samtools view Aligned.sortedByCoord.out.bam| grep -w 'NH:i:1') \
  |samtools view -Sb - > ${lane}_Aligned.sortedByCoord.uniq.bam

  # Index the BAM file
  samtools index ${lane}_Aligned.sortedByCoord.uniq.bam

  qualimap bamqc --java-mem-size=16G -bam ${lane}_Aligned.sortedByCoord.uniq.bam -outdir qualimap_${sample}_${lane}
  """
}

process RUN_MULTIQC {
  tag 'multiqc'

  publishDir "${params.output_dir}/multiqc", mode: 'copy'

  input:
  path qualimap_results

  output:
  path '*'

  script:
  """
  multiqc $qualimap_results -o .
  """
}

process MERGE_BAMS {
    tag "$sample"

    cpus = 8
    memory = 20.GB
    time = '48h'

    input:
    tuple val(sample), path(bams)

    output:
    tuple val(sample), path("${sample}_merged.bam"), path("${sample}_merged.bam.bai")

    publishDir "${params.output_dir}/bam_files", mode: 'copy'

    script:
    """
    # Merge BAM files from different lanes
    samtools merge -f ${sample}_merged.bam ${bams.join(' ')}

    # Index the merged BAM file
    samtools index ${sample}_merged.bam
    """
}

/* 
* PART 1: BASE EDITING MODULES 
*/

process RNASEQ_GATK_SPLITNCIGAR {

  cpus = 8
  memory = 20.GB
  time = '10h'

  tag "$replicateId"
  label 'mem_large'
  
  input: 
    path genome
    path index
    path dict
    tuple val(replicateId), path(bam), path(index)

  output:
    tuple val(replicateId), path('1_split.bam'), path('1_split.bai')
  
  publishDir "${params.output_dir}/hypertribe/${sample}/", mode: 'copy'

  script:
  """
  # SplitNCigarReads and reassign mapping qualities
  gatk SplitNCigarReads \
            -R $genome \
            -I $bam \
            --refactor-cigar-string \
            -O 1_split.bam
  """
}

/*
 * Process 4: Base recalibrate to detect systematic errors in base quality scores, 
 *            select unique alignments and index
 *             
 */
process RNASEQ_GATK_RECALIBRATE {

  cpus = 8
  memory = 20.GB
  time = '10h'

  tag "$replicateId"
  label "mem_large"

  input: 
    path genome
    path index
    path dict
    tuple val(sample), path(bam), path(index)
    path variants_file
    path variants_file_index

  output:
    tuple \
      val(sample), \
      path("2_final.rnaseq.grp"), \
      path("3_recalibrated.bam"), \
      path("3_recalibrated.bam.bai")

  publishDir "${params.output_dir}/hypertribe/${sample}/", mode: 'copy'
  
  script: 
  """
  # Indel Realignment and Base Recalibration
  gatk BaseRecalibrator \
          -R $genome \
          -I $bam \
          --known-sites $variants_file \
          -O 2_final.rnaseq.grp 

  gatk ApplyBQSR \
          -R $genome -I $bam \
          --bqsr-recal-file 2_final.rnaseq.grp \
          -O 3_recalibrated.bam

  # Index BAM files
  samtools index 3_recalibrated.bam
  """
}

/*
 * Process 5: Call variants with GATK HaplotypeCaller.
 *            Calls SNPs and indels simultaneously via local de-novo assembly of 
 *            haplotypes in an active region.
 *            Filter called variants with GATK VariantFiltration.    
 */
process RNASEQ_CALL_VARIANTS {

  cpus = 8
  memory = 20.GB
  time = '24h'

  tag "$sampleId"
  label "mem_xlarge"

  input:
    path genome
    path index
    path dict
    tuple val(sampleId), path(rnaseq_grp), path(bam), path(bai)
 
  output: 
    tuple val(sampleId), \
    path('4_output.gatk.vcf.gz'), \
    path('5_final.vcf')

  publishDir "${params.output_dir}/hypertribe/${sampleId}/", mode: 'copy'

  script:
  def bam_params = bam.collect{ "-I $it" }.join(' ')
  """
  # fix absolute path in dict file
  sed -i 's@UR:file:.*${genome}@UR:file:${genome}@g' $dict
  
  # Variant calling
  gatk HaplotypeCaller \
          --native-pair-hmm-threads ${task.cpus} \
          --reference ${genome} \
          --output 4_output.gatk.vcf.gz \
          ${bam_params} \
          --standard-min-confidence-threshold-for-calling 20.0 \
          --dont-use-soft-clipped-bases 

  # Variant filtering
  gatk VariantFiltration \
          -R ${genome} -V 4_output.gatk.vcf.gz \
          --cluster-window-size 35 --cluster-size 3 \
          --filter-name FS --filter-expression "FS > 30.0" \
          --filter-name QD --filter-expression "QD < 2.0" \
          -O 5_final.vcf
  """
}

 /*
 * Process 6: Process to format the results
 *            the process runs an R script that gathers the files from 3_variant step
 */


process FORMAT_OUTPUT {

  cpus = 8
  memory = 20.GB
  time = '20h'


  publishDir "${params.output_dir}/", mode: 'copy'

  input:
    val gtf_path
    tuple val(samples), path(vcf_path_list, stageAs: "?/*"), path(bam_path_list, stageAs: "?/*"), path(bai_path_list, stageAs: "?/*")
  


  output:
    path '4_formatted_output.csv'
  //path '*.rds'

  script:
  """
  echo \$(pwd)

  Rscript $baseDir/scripts/base_editing/4_format_output.R \
  --working_dir \$(pwd) \
  --gtf_path $gtf_path \
  --samples $samples \
  --bam_path_list ${bam_path_list.join(',')} \
  --vcf_path_list ${vcf_path_list.join(',')} \
  --genome_build ${params.genome_build} \
  --editor ${params.editor} ${workflow.resume ? '--resume' : ''}
  """
}

/*
 * Process 7: Base editing stats and comparisons 
 *            the process runs a R script that does differential analysis based on sample sheet
 */

process BE_STATS {

  cpus = 8
  memory = 20.GB
  time = '20h'

  publishDir "${params.output_dir}/deseq2/base_edit/", mode: 'copy'
  tag "BE_STATS"

  input:
  tuple val(meta), path (formatted_file)
  

  output:
  tuple val(meta),path ('*csv')

  script:
  def ctrl      = "${meta.ctrl}"
  def condition = "${meta.cond}"
  """
  Rscript $baseDir/scripts/base_editing/5_test_differential.R \
     --formatted_file ${formatted_file} \
     --cond ${condition} \
     --ctrl ${ctrl}
  """

}

/* 
* PART 2: RNA-seq expression MODULES 
*/

/*
 * Process 1 : RNA seq counts from merged bams 
 *                 using htseq counts         
 */

process RNASEQ_COUNT {
  tag "$sample"
  label "mem_xlarge"

  cpus = 8
  memory = 20.GB
  time = '48h'

  publishDir "${params.output_dir}/deseq2/count/", mode: 'copy'

  input:
    path gtf_path
    tuple val(sample), path(bam), path(bai)

  output:
    tuple val(sample), path('*counts')

  script:
  """
  htseq-count \
  -f bam \
  -r pos \
  -s reverse \
  -a 10 \
  -i gene_id \
  -m union \
  ${bam} ${gtf_path} > "${sample}.counts"
  """
}

/*
 * Process 2 : Gathering all the sample counts from htseq  
 *             using R script into a single file      
 */

process FORMAT_OUTPUT_COUNTS {

  cpus = 8
  memory = 20.GB
  time = '24h'

  

  publishDir "${params.output_dir}/deseq2/", mode: 'copy'
  tag "FORMAT_OUTPUT_COUNTS"

  input:
    val gtf_path
    tuple val(samples), path(counts_path_list, stageAs: "?/*")
  

  output:
    path '4_formatted_output_counts.csv'
  //path '*.rds'

  script:
  """

  echo \$(pwd)

  Rscript $baseDir/scripts/deseq2/4_format_output.R \
  --gtf_path $gtf_path \
  --samples ${samples.join(',')} \
  --counts_path_list ${counts_path_list.join(',')} \
  --genome_build ${params.genome_build} 
  """

}

/*
 * Process 3 : Differential gene expression using DESeq2     
 */

process DESEQ2_RNA_SEQ {

  cpus = 8
  memory = 20.GB
  time = '12h'

  publishDir "${params.output_dir}/deseq2/DEG/", mode: 'copy'
  tag "DESEQ2"

  input:
  tuple val(meta), path (counts_file)

  output:
  tuple val(meta), path ('*csv')

  script:
  def ctrl      = "${meta.ctrl}"
  def condition = "${meta.cond}"
  """
  Rscript $baseDir/scripts/deseq2/5_test_differential.R \
     --counts_file ${counts_file} \
     --cond ${condition} \
     --ctrl ${ctrl}
  """
}




process PLOTS_RNA_SEQ {

  publishDir "${params.output_dir}/plots/", mode: 'copy'
  tag "PLOTS"

  cpus = 8
  memory = 20.GB
  time = '12h'

  input:
    val (threshold)
    tuple val(meta), path (deseq_file)
    val (formatted_file)
    val (analysis)

  output:
    tuple val(meta), path ('*jpeg') 

  script:
  def ctrl      = "${meta.ctrl}"
  def condition = "${meta.cond}"
  """
  Rscript $baseDir/scripts/plots_pipeline.R \
     --deseq_file ${deseq_file} \
     --formatted_file ${formatted_file} \
     --cond ${condition} \
     --ctrl ${ctrl} \
     --analysis ${analysis} \
     --padj_threshold ${params.padj_threshold} \
     --diff_threshold ${threshold} \
     --selection_name "NULL" \
     --project_name ${params.project_name}
  """
}
