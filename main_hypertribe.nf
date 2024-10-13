/*
 * Copyright (c) 2020, Seqera Labs.
 * Copyright (c) 2017-2019, Centre for Genomic Regulation (CRG).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

/*
 * 'CalliNGS-NF' - A Nextflow pipeline for variant calling with NGS data
 *
 * This pipeline that reproduces steps from the GATK best practics of SNP
 * calling with RNAseq data procedure:
 * https://software.broadinstitute.org/gatk/guide/article?id=3891
 *
 * Anna Vlasova
 * Emilio Palumbo
 * Paolo Di Tommaso
 * Evan Floden
 */

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */

params.genome_dir = "$baseDir/genome"
params.genome     = "${params.genome_dir}/genome.fa"
params.genome_index     = "${params.genome_dir}/genome.fa.fai"
params.genome_dict = "${params.genome_dir}/genome.dict"
params.genomeDir = "${params.genome_dir}/genome_index/"
params.variants   = "${params.genome_dir}/known_variants.filtered.recode.vcf.gz"
params.variants_index   = "${params.genome_dir}/known_variants.filtered.recode.vcf.gz.tbi"
//params.gtf        = "${params.genome_dir}/chr22_annotation.gtf"
//params.input_dir      = "$baseDir/Test_Hypertribe/input_data/*/*_L00{1,2}_R{1,2}.fq.gz"
params.genome_build  = 'hg38'
params.editor        = 'ADAR'
params.output_dir    = "$baseDir/Test_Hypertribe/output_data/"
params.sample_sheet  = "$baseDir/Test_Hypertribe/sample_sheet_deseq2.csv"


log.info """\
C A L L I N G S  -  N F    v 2.1
================================
  genome   : ${params.genome}
  genome_index : ${params.genome_index}
  variants : ${params.variants}
  gtf      : ${params.gtf}
  results  : ${params.output_dir}
  genome_build : ${params.genome_build}
  editor : ${params.editor}

}
"""

/*
 * Import modules
 */
include {
  RNASEQ_MAPPING_STAR;
  RUN_MULTIQC;
  MERGE_BAMS;
  } from './modules.nf'


include { HYPERTRIBE } from './workflows/HYPERTRIBE.nf'
include { DESEQ2 }     from './workflows/DESEQ2.nf'

/*
 * main pipeline logic
 */


def group_per_sample = { channel ->
  channel.groupTuple(by: [0])
}

workflow {
      
      // ** -- Read sampled from the folder
      // reads_ch = Channel.fromFilePairs(params.input_dir,
      //     size: 2,
      //     flat: false)
      //     { file -> file.name.split('_R')[0] }
      //     .map { sample, files ->
      //           def lane = files[0].name =~ /L00[1-4]/
      //           [sample.split('_')[0], lane[0], files.sort().collect { it.toString() }]
      //       }


      // --- Real part -- //

      // reads_ch = Channel.fromFilePairs(params.input_dir,
      //     size: 2,
      //     flat: false)
      //     { file -> file.name.split('_R')[0] }
      //     .map { sample, files ->
      //           def lane = files[0].name =~ /S.*_L00[1-4]/ 
      //           [sample.split('_')[0], lane[0], files.sort().collect { it.toString() }]
      //       }

      //reads_ch.view { "reads_ch: ${it}"}



      // ** -- Read the sample sheet
      deseq_ch = join_csv(file(params.sample_sheet))
      
      // ** -- Read the sample names sheet 
      samplesheet_ch = sample_names(file(params.sample_names, checkIfExists: true))
      samplesheet_ch.view()

      // PART 1: STAR RNA-Seq Mapping
      rnaseq_map_ch = RNASEQ_MAPPING_STAR(
            params.genome,
            params.genomeDir,
            samplesheet_ch)

       
      mapped_ch = RNASEQ_MAPPING_STAR.out
                  .map { id, sample_lane, files1, files2, files3 ->
                       def sample = sample_lane.split('_')[0]  // Splitting sample_lane to extract the sample
                       //def lane = sample_lane.split('_')[1]
                       def combined_id_sample = "${id}_${sample}"  // Combine id and sample
                       return [id, sample_lane,files1, files2, files3 ]}
      mapped_ch.view()

      // // PART 2: Run MultiQC on all qualimap results
      // //RUN_MULTIQC(rnaseq_map_ch.collect { it[4] })
  
      // Merge BAM files from different lanes
      merge_ch_input = group_per_sample(mapped_ch).map { s ->
        [s[0],s[2]]}
      

      merge_ch_input.view { "merge_ch_input: ${it}"}

      MERGE_BAMS(merge_ch_input)
      merged_bams_ch = MERGE_BAMS.out

      // *-------- BASE EDITING SCRIPTS WORKFLOW -------* //
      if (params.workflow_HY == "HYPERTRIBE") {
        HYPERTRIBE(merged_bams_ch,deseq_ch)
      }

      if (params.workflow_DE == "DESEQ2") {
        DESEQ2(merged_bams_ch,deseq_ch)
      }
      
  }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to extract information (meta data + file(s)) from csv file(s)
def join_csv(csv_file) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }

    Channel.of(csv_file).splitCsv(header: true)
        .map{ row ->
           def condition  = row['Condition']
           def control = row['Control']
           return [condition,control]
        }.map{condition, control ->
          def meta = [
                   "cond" : condition,
                   "ctrl" : control 
          ]
          return[meta]}
         

}


def sample_names(csv_file) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }

    Channel.of(csv_file).splitCsv(header: true)
        .map{ row ->
           def group_name  = row['Group']
           def sample_ID = row['Sample']
           def fastq_1 = "${params.input_dir_new}${row['Folder_name']}/${row['fastq_1']}"  // Combine input directory 
           def fastq_2 = "${params.input_dir_new}${row['Folder_name']}/${row['fastq_2']}" 
           return [group_name,sample_ID,fastq_1,fastq_2]
        }
         
}