/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XBS_VARIANT_CALLING — top-level composer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Glues XBS_PER_SAMPLE (stages 1–7) to XBS_COHORT (stages 8–12).

    Consumed by:
      - MAGMA (Plan-3 §11 integration) — drives the Heupink 2021 spine for MTB
      - mtbc-varcaller-nf (Phase-1) — wraps this as the `gatk_vqsr` backend
      - any other haploid-bacterial caller pipeline

    The emit names below are the STABLE INTERFACE CONTRACT. Renames are
    breaking changes for both MAGMA and mtbc-varcaller-nf. See docs/CONTRACT.md.
*/

include { XBS_PER_SAMPLE   } from './xbs_per_sample'
include { XBS_COHORT       } from './xbs_cohort'


workflow XBS_VARIANT_CALLING {

    take:
    ch_reads        // channel: [ meta(study, sample, library, id), [r1, r2] ]
    ch_reference    // value:   [ meta(id:'ref'), fasta, fai, dict, [amb, ann, bwt, pac, sa] ]
    ch_snp_truth    // value:   [ meta(id:'snp_truth'),   vcf, tbi ]   (required iff params.snp_filter_mode == 'vqsr')
    ch_indel_truth  // value:   [ meta(id:'indel_truth'), vcf, tbi ]   (required iff params.indel_filter_mode == 'vqsr')
    ch_dbsnp        // value:   [ meta(id:'dbsnp'),       vcf, tbi ]   (empty unless !params.skip_bqsr)

    main:

    XBS_PER_SAMPLE(ch_reads, ch_reference, ch_dbsnp)

    // Cohort stages are gated by params.skip_cohort. When skipped, downstream
    // emits are empty channels — consumers can still wire to them safely.
    def ch_raw_variants     = channel.empty()
    def ch_snp_filtered     = channel.empty()
    def ch_indel_filtered   = channel.empty()
    def ch_vqsr_diagnostics = channel.empty()

    if (!params.skip_cohort) {
        XBS_COHORT(
            XBS_PER_SAMPLE.out.gvcf_vcf,
            XBS_PER_SAMPLE.out.gvcf_tbi,
            ch_reference,
            ch_snp_truth,
            ch_indel_truth
        )
        ch_raw_variants     = XBS_COHORT.out.raw_variants
        ch_snp_filtered     = XBS_COHORT.out.snp_filtered
        ch_indel_filtered   = XBS_COHORT.out.indel_filtered
        ch_vqsr_diagnostics = XBS_COHORT.out.vqsr_diagnostics
    }

    emit:
    // ============= STABLE INTERFACE CONTRACT (see docs/CONTRACT.md) =============
    // Per-sample outputs (always emitted)
    gvcfs              = XBS_PER_SAMPLE.out.gvcf_vcf.join(XBS_PER_SAMPLE.out.gvcf_tbi)   // [ meta, vcf, tbi ]
    sample_bam         = XBS_PER_SAMPLE.out.sample_bam                                   // [ meta, bam, bai ]
    samtools_stats     = XBS_PER_SAMPLE.out.samtools_stats                               // [ meta, *.stats    ]
    flagstat           = XBS_PER_SAMPLE.out.flagstat                                     // [ meta, *.flagstat ]
    wgs_metrics        = XBS_PER_SAMPLE.out.wgs_metrics                                  // [ meta, *_metrics  ]
    markdup_metrics    = XBS_PER_SAMPLE.out.markdup_metrics                              // [ meta, *.metrics  ]

    // Cohort outputs — empty channels when params.skip_cohort=true
    raw_variants       = ch_raw_variants                                                 // [ meta, vcf, tbi ]
    snp_filtered       = ch_snp_filtered                                                 // [ meta, vcf, tbi ] (vqsr | hard_filters | none)
    indel_filtered     = ch_indel_filtered                                               // [ meta, vcf, tbi ] (vqsr | hard_filters | none)
    vqsr_diagnostics   = ch_vqsr_diagnostics                                             // tranches/recal/.R — only when any *_filter_mode == 'vqsr'
}
