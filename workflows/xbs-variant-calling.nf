/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XBS_VARIANT_CALLING — top-level pipeline workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Builds reference / truth-set / dbsnp value channels from params and calls
    the XBS_VARIANT_CALLING subworkflow with the samplesheet stream.
*/

include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_xbs-variant-calling_pipeline'
include { XBS_VARIANT_CALLING as XBS_SPINE      } from '../subworkflows/local/xbs_variant_calling_wf'


workflow XBS_VARIANT_CALLING {

    take:
    ch_samplesheet // per-row [ meta(study, sample, library), [r1] | [r1, r2] ]
    outdir

    main:

    def ch_versions = channel.empty()

    //
    // Reference bundle — single value channel of [meta, fasta, fai, dict, [bwa index files]]
    //
    def ref_fasta = file("${params.reference_dir}/${params.reference_basename}.fa", checkIfExists: true)
    def ref_fai   = file("${params.reference_dir}/${params.reference_basename}.fa.fai", checkIfExists: true)
    def ref_dict  = file("${params.reference_dir}/${params.reference_basename}.dict", checkIfExists: true)
    def bwa_index = [
        file("${params.reference_dir}/${params.reference_basename}.fa.amb", checkIfExists: true),
        file("${params.reference_dir}/${params.reference_basename}.fa.ann", checkIfExists: true),
        file("${params.reference_dir}/${params.reference_basename}.fa.bwt", checkIfExists: true),
        file("${params.reference_dir}/${params.reference_basename}.fa.pac", checkIfExists: true),
        file("${params.reference_dir}/${params.reference_basename}.fa.sa",  checkIfExists: true),
    ]
    def ch_reference = channel.value([ [id: 'ref'], ref_fasta, ref_fai, ref_dict, bwa_index ])

    //
    // Truth sets — required only when the corresponding filter_mode == 'vqsr'.
    // dbsnp_vcf required only when !skip_bqsr.
    //
    def ch_snp_truth
    if (params.snp_filter_mode == 'vqsr') {
        def snp_truth_vcf = file(params.snp_truth_vcf, checkIfExists: true)
        def snp_truth_tbi = file(params.snp_truth_vcf_tbi ?: "${params.snp_truth_vcf}.tbi", checkIfExists: true)
        ch_snp_truth = channel.value([ [id: 'snp_truth'], snp_truth_vcf, snp_truth_tbi ])
    } else {
        ch_snp_truth = channel.value([ [id: 'snp_truth'], [], [] ])
    }

    def ch_indel_truth
    if (params.indel_filter_mode == 'vqsr') {
        def indel_truth_vcf = file(params.indel_truth_vcf, checkIfExists: true)
        def indel_truth_tbi = file(params.indel_truth_vcf_tbi ?: "${params.indel_truth_vcf}.tbi", checkIfExists: true)
        ch_indel_truth = channel.value([ [id: 'indel_truth'], indel_truth_vcf, indel_truth_tbi ])
    } else {
        ch_indel_truth = channel.value([ [id: 'indel_truth'], [], [] ])
    }

    def ch_dbsnp
    if (!params.skip_bqsr) {
        def dbsnp_vcf = file(params.dbsnp_vcf, checkIfExists: true)
        def dbsnp_tbi = file(params.dbsnp_vcf_tbi ?: "${params.dbsnp_vcf}.tbi", checkIfExists: true)
        ch_dbsnp = channel.value([ [id: 'dbsnp'], dbsnp_vcf, dbsnp_tbi ])
    } else {
        ch_dbsnp = channel.value([ [id: 'dbsnp'], [], [] ])
    }

    //
    // Reshape the samplesheet stream so meta.id is unique per library
    // (study.sample.library) for the per-library mapping stage.
    //
    def ch_reads = ch_samplesheet.map { meta, reads ->
        def new_meta = meta + [ id: "${meta.study}.${meta.sample}.${meta.library}" ]
        [ new_meta, reads ]
    }

    //
    // Run the XBS spine
    //
    XBS_SPINE(ch_reads, ch_reference, ch_snp_truth, ch_indel_truth, ch_dbsnp)

    //
    // Collate software versions (standard nf-core boilerplate, kept verbatim)
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'xbs-variant-calling_software_versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    gvcfs            = XBS_SPINE.out.gvcfs
    sample_bam       = XBS_SPINE.out.sample_bam
    samtools_stats   = XBS_SPINE.out.samtools_stats
    flagstat         = XBS_SPINE.out.flagstat
    wgs_metrics      = XBS_SPINE.out.wgs_metrics
    markdup_metrics  = XBS_SPINE.out.markdup_metrics
    raw_variants     = XBS_SPINE.out.raw_variants
    snp_filtered     = XBS_SPINE.out.snp_filtered
    indel_filtered   = XBS_SPINE.out.indel_filtered
    vqsr_diagnostics = XBS_SPINE.out.vqsr_diagnostics
    versions         = ch_versions
}
