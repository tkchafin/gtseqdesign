/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_acamel-gtseqdesign_pipeline/main.nf'
include { ADMIXPIPE as ADMIXPIPE_PRE } from '../subworkflows/local/admixpipe.nf'
include { ADMIXPIPE as ADMIXPIPE_POST } from '../subworkflows/local/admixpipe.nf'
include { SELECT_CANDIDATES } from '../subworkflows/local/select_candidates.nf'
include { GENERATE_REPORT } from '../subworkflows/local/generate_report.nf'
include { SNPIO_PRE_FILTER as SNPIO_FILTER } from '../modules/local/snpio/pre_filter.nf'
include { LIST_CHROMS } from '../modules/local/list_chroms.nf'
include { GENERATE_CONSENSUS } from '../modules/local/generate_consensus.nf'
include { FILTER_POSITIONS } from '../modules/local/filter_positions.nf'
include { CUSTOMIZE_REPORT } from '../modules/local/report/customize_report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ACAMEL_GTSEQDESIGN {

    take:
    ch_vcf     // [meta, vcf]
    ch_tbi     // [meta, tbi]
    ch_popmap  // [meta, popmap]
    ch_reference

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Process reference catalog if needed
    //
    ch_reference
        .branch {
            loci: it[1].name.endsWith('.loci')
            fasta: true
        }
        .set { ch_ref_to_process }
    GENERATE_CONSENSUS ( ch_ref_to_process.loci )
    ch_ref_ready = ch_ref_to_process.fasta
        | mix ( GENERATE_CONSENSUS.out.fasta )

    //
    // VCF pre-processing
    //
    // This step removes individuals with a large amount of missing data,
    // flanking variation within ${params.primer_length} distance, low variation,
    // and generates SNPio missingness reports
    SNPIO_FILTER(
        ch_vcf,
        ch_tbi,
        ch_popmap
    )
    ch_versions = ch_versions.mix(SNPIO_FILTER.out.versions)
    ch_filtered_vcf = SNPIO_FILTER.out.filtered_vcf.map { meta, file -> tuple(meta + [id: "${meta.id}_filtered"], file) }
    ch_filtered_tbi = SNPIO_FILTER.out.filtered_tbi.map { meta, file -> tuple(meta + [id: "${meta.id}_filtered"], file) }
    ch_snpio_output = SNPIO_FILTER.out.snpio_output.map { meta, dir -> tuple(meta + [id: "${meta.id}_filtered"], dir) }

    //
    // Run admixture pipeline on full (filtered) dataset
    //
    ADMIXPIPE_PRE(
        ch_filtered_vcf,
        ch_popmap
    )
    ch_versions = ch_versions.mix(ADMIXPIPE_PRE.out.versions)

    //
    // Denovo assembly handling
    //
    // Removes SNPs if they are not at least primer_length from locus boundaries
    // number of bases (for denovo assembled loci only)
    // This forces primer+variant to be fully contained ONLY within the sequences locus
    // Otherwise, flanking sequence will be inferred using the provided reference
    if ( params.fully_contained ) {
        // Get list of denovo loci and their lengths
        LIST_CHROMS( ch_ref_ready )

        FILTER_POSITIONS(
            ch_filtered_vcf,
            ch_filtered_tbi,
            LIST_CHROMS.out.chroms
        )
        ch_versions = ch_versions.mix ( FILTER_POSITIONS.out.versions )

        ch_candidates = FILTER_POSITIONS.out.vcf
        ch_candidates_tbi = FILTER_POSITIONS.out.tbi
    } else {
        ch_candidates = ch_filtered_vcf
        ch_candidates_tbi = ch_filtered_tbi
    }

    //
    // Compute locus-wise importance metrics
    //
    SELECT_CANDIDATES(
        ch_candidates,
        ch_candidates_tbi,
        ADMIXPIPE_PRE.out.inds,
        ADMIXPIPE_PRE.out.bestK_clumpp,
        ADMIXPIPE_PRE.out.bestK
    )
    ch_versions = ch_versions.mix(SELECT_CANDIDATES.out.versions)
    ch_selected_vcf = SELECT_CANDIDATES.out.vcf.map { meta, file -> tuple(meta + [id: meta.id.replaceFirst(/_filtered$/, '_selected')], file) }
    ch_selected_tbi = SELECT_CANDIDATES.out.tbi.map { meta, file -> tuple(meta + [id: meta.id.replaceFirst(/_filtered$/, '_selected')], file) }
    ch_selected_snpio_output = SELECT_CANDIDATES.out.snpio_output.map { meta, dir -> tuple(meta + [id: meta.id.replaceFirst(/_filtered$/, '_selected')], dir) }

    //
    // Run admixture pipeline on selected candidates
    //
    ADMIXPIPE_POST(
        ch_selected_vcf,
        ch_popmap
    )
    ch_versions = ch_versions.mix(ADMIXPIPE_POST.out.versions)

    //
    // Generate figures for the report
    //
    GENERATE_REPORT(
        ch_vcf,
        ch_tbi,
        ch_filtered_vcf,
        ch_filtered_tbi,
        ch_selected_vcf,
        ch_selected_tbi,
        ADMIXPIPE_PRE.out.cv_file,
        ch_snpio_output,
        ch_selected_snpio_output,
        ADMIXPIPE_PRE.out.bestK_clumpp,
        ADMIXPIPE_POST.out.bestK_clumpp,
        ADMIXPIPE_POST.out.inds,
        ADMIXPIPE_POST.out.pops,
        SELECT_CANDIDATES.out.metrics,
        SELECT_CANDIDATES.out.top_loci
    )
    ch_versions = ch_versions.mix( GENERATE_REPORT.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( GENERATE_REPORT.out.mqc_files )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    // Customize report
    CUSTOMIZE_REPORT( MULTIQC.out.report )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
