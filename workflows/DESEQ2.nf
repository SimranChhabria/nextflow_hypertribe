/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

include {
  RNASEQ_COUNT;
  FORMAT_OUTPUT_COUNTS;
  DESEQ2_RNA_SEQ;
  PLOTS_RNA_SEQ as PLOTS_DESEQ2_SEQ;
  } from '../modules.nf'

workflow DESEQ2 {

    take:
    merge_bams_ch
    deseq_ch

    main:
    RNASEQ_COUNT(params.gtf,
                   merge_bams_ch)

    rna_ch = RNASEQ_COUNT.out

    names_ch = rna_ch.map { it[0] }.toList()
    files_ch = rna_ch.map { it[1] }.toList()
      
    combined_ch = names_ch.collect().toList()
                    .combine(files_ch.collect().toList())
  

    //-- Format the output --//

    FORMAT_OUTPUT_COUNTS(
         params.gtf,
         combined_ch
       )

    //-- Differential Gene Expression Analysis --//

    counts_ch = FORMAT_OUTPUT_COUNTS.out.view()
    deseq2_final_ch = deseq_ch
                      .combine(counts_ch).view()
    
    DESEQ2_RNA_SEQ(
        deseq2_final_ch 
      )

    deseq2_out_ch = DESEQ2_RNA_SEQ.out.view()

    //-- Generate the plots --//

    PLOTS_DESEQ2_SEQ(
         params.diff_threshold_deseq2,
         DESEQ2_RNA_SEQ.out,
         counts_ch.view().first(),
         "deseq2" 
       )

}