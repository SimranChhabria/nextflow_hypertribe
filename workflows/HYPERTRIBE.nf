/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

include {
  RNASEQ_GATK_SPLITNCIGAR;
  RNASEQ_GATK_RECALIBRATE;
  RNASEQ_CALL_VARIANTS;
  FORMAT_OUTPUT;
  BE_STATS;
  PLOTS_RNA_SEQ as PLOTS_HYPERTRIBE_SEQ;
  } from '../modules.nf'

workflow HYPERTRIBE {

    take:
    merge_bams_ch
    deseq_ch



    main:
    // *-------- BASE EDITING SCRIPTS -------* //

    // PART 3: GATK Prepare Mapped Reads
    RNASEQ_GATK_SPLITNCIGAR(
            params.genome,
            params.genome_index,
            params.genome_dict,
            merge_bams_ch )

    // PART 4: GATK Base Quality Score Recalibration Workflow
    RNASEQ_GATK_RECALIBRATE(
                  params.genome,
                  params.genome_index,
                  params.genome_dict,
                  RNASEQ_GATK_SPLITNCIGAR.out,
                  params.variants,
                  params.variants_index)

    // PART 5: GATK Variant Calling
    RNASEQ_CALL_VARIANTS(
            params.genome,
            params.genome_index,
            params.genome_dict,
            RNASEQ_GATK_RECALIBRATE.out )
      


    format_ch = RNASEQ_CALL_VARIANTS.out.join(RNASEQ_GATK_RECALIBRATE.out)
        .map { it.flatten() }
        .toList()
        .map { it.transpose() }
        .map { [it[0].join(','), it[1], it[2], it[3]] }
        //.view{"format_ch: ${it}"}
        // format_ch = RNASEQ_CALL_VARIANTS.out.join(RNASEQ_GATK_RECALIBRATE.out).toList().map { it.transpose() }.map { it[0].join(','), it[1], it[2], it[3] }.view{"format_ch: ${it}"}
       
    // PART 6: Format the output

    FORMAT_OUTPUT(
        params.gtf,
        format_ch)


    // Part 7: Differential Analysis for base editing
      
    be_ch = FORMAT_OUTPUT.out.view()     
    // Combine channels to form meta and output file
    be_final_ch = deseq_ch
                      .combine(be_ch).view()
    BE_STATS(
         be_final_ch
       )

    // Part 8:Generating the plots
    PLOTS_HYPERTRIBE_SEQ(
        params.diff_threshold_hyper,
        BE_STATS.out,
        be_ch.view().first(),
        "hypertribe")

}