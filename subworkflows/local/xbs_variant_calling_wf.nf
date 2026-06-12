/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XBS_VARIANT_CALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Top-level composer of the XBS spine: glues XBS_PER_SAMPLE
    (stages 1–7) to XBS_COHORT (stages 8–12).

    This is the subworkflow that other pipelines (MAGMA, future
    NF_CORE_TBANALYZER, M. bovis / other bacterial pipelines) import
    via `nf-core subworkflows install` once published.
*/

include { XBS_PER_SAMPLE   } from './xbs_per_sample'
include { XBS_COHORT  } from './xbs_cohort'


workflow XBS_VARIANT_CALLING {

    take:
    ch_reads        // channel: [ meta(study, sample, library, id), [r1, r2] ]
    ch_reference    // value:   [ meta(id:'ref'), fasta, fai, dict, [amb, ann, bwt, pac, sa] ]
    ch_snp_truth    // value:   [ meta(id:'snp_truth'),   vcf, tbi ]   (required)
    ch_indel_truth  // value:   [ meta(id:'indel_truth'), vcf, tbi ]   (required)
    ch_dbsnp        // value:   [ meta(id:'dbsnp'),       vcf, tbi ]   (empty unless !params.skip_bqsr)

    main:

    XBS_PER_SAMPLE(ch_reads, ch_reference, ch_dbsnp)

    XBS_COHORT(
        XBS_PER_SAMPLE.out.gvcf_vcf,
        XBS_PER_SAMPLE.out.gvcf_tbi,
        ch_reference,
        ch_snp_truth,
        ch_indel_truth
    )

    emit:
    // per-sample outputs
    gvcfs              = XBS_PER_SAMPLE.out.gvcf_vcf.join(XBS_PER_SAMPLE.out.gvcf_tbi)   // [ meta, vcf, tbi ]
    sample_bam         = XBS_PER_SAMPLE.out.sample_bam                                   // [ meta, bam, bai ]
    samtools_stats     = XBS_PER_SAMPLE.out.samtools_stats
    flagstat           = XBS_PER_SAMPLE.out.flagstat
    wgs_metrics        = XBS_PER_SAMPLE.out.wgs_metrics
    markdup_metrics    = XBS_PER_SAMPLE.out.markdup_metrics

    // cohort outputs (the XBS primary deliverables)
    raw_variants       = XBS_COHORT.out.raw_variants
    snp_filtered       = XBS_COHORT.out.snp_filtered
    indel_filtered     = XBS_COHORT.out.indel_filtered
    vqsr_diagnostics   = XBS_COHORT.out.vqsr_diagnostics
}
