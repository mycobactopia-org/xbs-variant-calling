/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XBS_COHORT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Cohort-level joint calling + VQSR — the core Heupink 2021 step that
    enables low-coverage / contaminated-sample variant calling at scale.

    Stages 8–12 from the canonical XBS-core bash:
      8. GATK CombineGVCFs (cohort-wide)
      9. GATK GenotypeGVCFs (cohort-wide, ploidy=1)
     10. GATK SelectVariants split into raw SNPs + raw INDELs
     11. SNP VQSR (VariantRecalibrator + ApplyVQSR)
     12. INDEL VQSR (VariantRecalibrator + ApplyVQSR)

    Output: filtered SNP VCF + filtered INDEL VCF.
*/

include { GATK4_COMBINEGVCFS                                                } from '../../modules/nf-core/gatk4/combinegvcfs/main'
include { GATK4_GENOTYPEGVCFS                                               } from '../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_SNP                  } from '../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_INDEL                } from '../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_SNP        } from '../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_INDEL      } from '../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNP                            } from '../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDEL                          } from '../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_SNP            } from '../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_INDEL          } from '../../modules/nf-core/gatk4/variantfiltration/main'


workflow XBS_COHORT {

    take:
    ch_gvcf_vcf       // channel: [ per_sample_meta, *.g.vcf.gz ]
    ch_gvcf_tbi       // channel: [ per_sample_meta, *.g.vcf.gz.tbi ]
    ch_reference      // value:   [ meta(id:'ref'), fasta, fai, dict, [bwa_index_files] ]
    ch_snp_truth      // value:   [ meta(id:'snp_truth'),   vcf, tbi ]
    ch_indel_truth    // value:   [ meta(id:'indel_truth'), vcf, tbi ]

    main:

    // ---- reference channel adapters ----
    def ch_fasta_tuple = ch_reference.map { m, f, _fai, _d, _b -> [m, f] }
    def ch_fai_tuple   = ch_reference.map { m, _f, fai, _d, _b -> [m, fai] }
    def ch_dict_tuple  = ch_reference.map { m, _f, _fai, d, _b -> [m, d] }
    def ch_fasta_path  = ch_reference.map { _m, f, _fai, _d, _b -> f }
    def ch_fai_path    = ch_reference.map { _m, _f, fai, _d, _b -> fai }
    def ch_dict_path   = ch_reference.map { _m, _f, _fai, d, _b -> d }

    // ---- truth-set adapters ----
    def ch_snp_truth_vcf   = ch_snp_truth.map   { _m, v, _t -> v }
    def ch_snp_truth_tbi   = ch_snp_truth.map   { _m, _v, t -> t }
    def ch_indel_truth_vcf = ch_indel_truth.map { _m, v, _t -> v }
    def ch_indel_truth_tbi = ch_indel_truth.map { _m, _v, t -> t }

    //
    // Collect all per-sample GVCFs into ONE cohort tuple.
    // ch_combined_input shape: [ cohort_meta, [vcf, vcf, ...], [tbi, tbi, ...] ]
    //
    def cohort_id = params.joint_name ?: 'joint'
    def ch_combined_input = ch_gvcf_vcf
        .join(ch_gvcf_tbi)
        .map { _meta, vcf, tbi -> [vcf, tbi] }
        .collect(flat: false)
        .map { vcf_tbi_pairs ->
            def vcfs = vcf_tbi_pairs.collect { it[0] }
            def tbis = vcf_tbi_pairs.collect { it[1] }
            [ [id: cohort_id], vcfs, tbis ]
        }

    //
    // STAGE 8: CombineGVCFs (cohort)
    //
    GATK4_COMBINEGVCFS(ch_combined_input, ch_fasta_path, ch_fai_path, ch_dict_path)

    //
    // STAGE 9: GenotypeGVCFs (cohort, ploidy=1 via ext.args)
    //
    def ch_gg_input = GATK4_COMBINEGVCFS.out.combined_gvcf
        .join(GATK4_COMBINEGVCFS.out.tbi)
        .map { meta, vcf, tbi -> [meta, vcf, tbi, [], []] }   // empty intervals
    GATK4_GENOTYPEGVCFS(
        ch_gg_input,
        ch_fasta_tuple, ch_fai_tuple, ch_dict_tuple,
        channel.value([[id:'none'], []]),  // dbsnp
        channel.value([[id:'none'], []])   // dbsnp_tbi
    )

    def ch_raw_variants = GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi)

    //
    // STAGE 10: split into SNP / INDEL
    //
    def ch_sv_input = ch_raw_variants.map { meta, vcf, tbi -> [meta, vcf, tbi, []] }  // no intervals
    GATK4_SELECTVARIANTS_SNP(ch_sv_input)
    GATK4_SELECTVARIANTS_INDEL(ch_sv_input)

    //
    // STAGE 11: SNP filtering — three modes (vqsr | hard_filters | none)
    //
    def ch_snp_filtered     = channel.empty()
    def ch_indel_filtered   = channel.empty()
    def ch_vqsr_diagnostics = channel.empty()

    def ch_sv_snp_out = GATK4_SELECTVARIANTS_SNP.out.vcf.join(GATK4_SELECTVARIANTS_SNP.out.tbi)

    if (params.snp_filter_mode == 'vqsr') {
        def ch_vr_snp_in = ch_sv_snp_out
        def labels_snp = ['--resource:5000SNP,known=false,training=true,truth=true,prior=20.0']
        GATK4_VARIANTRECALIBRATOR_SNP(
            ch_vr_snp_in,
            ch_snp_truth_vcf, ch_snp_truth_tbi,
            labels_snp,
            ch_fasta_path, ch_fai_path, ch_dict_path
        )

        def ch_apply_snp_in = ch_sv_snp_out
            .join(GATK4_VARIANTRECALIBRATOR_SNP.out.recal)
            .join(GATK4_VARIANTRECALIBRATOR_SNP.out.idx)
            .join(GATK4_VARIANTRECALIBRATOR_SNP.out.tranches)
        GATK4_APPLYVQSR_SNP(ch_apply_snp_in, ch_fasta_path, ch_fai_path, ch_dict_path)

        ch_snp_filtered = GATK4_APPLYVQSR_SNP.out.vcf.join(GATK4_APPLYVQSR_SNP.out.tbi)
        ch_vqsr_diagnostics = ch_vqsr_diagnostics
            .mix(GATK4_VARIANTRECALIBRATOR_SNP.out.tranches)
            .mix(GATK4_VARIANTRECALIBRATOR_SNP.out.plots.ifEmpty([]))
    } else if (params.snp_filter_mode == 'hard_filters') {
        GATK4_VARIANTFILTRATION_SNP(
            ch_sv_snp_out,
            ch_fasta_tuple, ch_fai_tuple, ch_dict_tuple,
            channel.value([[id:'none'], []])   // gzi unused for bgzipped (non-gzi) VCFs
        )
        ch_snp_filtered = GATK4_VARIANTFILTRATION_SNP.out.vcf.join(GATK4_VARIANTFILTRATION_SNP.out.tbi)
    } else {
        // 'none' — pass-through unfiltered SelectVariants output
        ch_snp_filtered = ch_sv_snp_out
    }

    //
    // STAGE 12: INDEL filtering — three modes (vqsr | hard_filters | none)
    //
    def ch_sv_indel_out = GATK4_SELECTVARIANTS_INDEL.out.vcf.join(GATK4_SELECTVARIANTS_INDEL.out.tbi)

    if (params.indel_filter_mode == 'vqsr') {
        def ch_vr_indel_in = ch_sv_indel_out
        def labels_indel = ['--resource:500INDEL,known=false,training=true,truth=true,prior=20.0']
        GATK4_VARIANTRECALIBRATOR_INDEL(
            ch_vr_indel_in,
            ch_indel_truth_vcf, ch_indel_truth_tbi,
            labels_indel,
            ch_fasta_path, ch_fai_path, ch_dict_path
        )

        def ch_apply_indel_in = ch_sv_indel_out
            .join(GATK4_VARIANTRECALIBRATOR_INDEL.out.recal)
            .join(GATK4_VARIANTRECALIBRATOR_INDEL.out.idx)
            .join(GATK4_VARIANTRECALIBRATOR_INDEL.out.tranches)
        GATK4_APPLYVQSR_INDEL(ch_apply_indel_in, ch_fasta_path, ch_fai_path, ch_dict_path)

        ch_indel_filtered = GATK4_APPLYVQSR_INDEL.out.vcf.join(GATK4_APPLYVQSR_INDEL.out.tbi)
        ch_vqsr_diagnostics = ch_vqsr_diagnostics
            .mix(GATK4_VARIANTRECALIBRATOR_INDEL.out.tranches)
            .mix(GATK4_VARIANTRECALIBRATOR_INDEL.out.plots.ifEmpty([]))
    } else if (params.indel_filter_mode == 'hard_filters') {
        GATK4_VARIANTFILTRATION_INDEL(
            ch_sv_indel_out,
            ch_fasta_tuple, ch_fai_tuple, ch_dict_tuple,
            channel.value([[id:'none'], []])
        )
        ch_indel_filtered = GATK4_VARIANTFILTRATION_INDEL.out.vcf.join(GATK4_VARIANTFILTRATION_INDEL.out.tbi)
    } else {
        // 'none' — pass-through
        ch_indel_filtered = ch_sv_indel_out
    }

    emit:
    raw_variants     = ch_raw_variants     // [ meta, raw_variants.vcf.gz, tbi ]
    snp_filtered     = ch_snp_filtered     // [ meta, vcf, tbi ] — content depends on snp_filter_mode
    indel_filtered   = ch_indel_filtered   // [ meta, vcf, tbi ] — content depends on indel_filter_mode
    vqsr_diagnostics = ch_vqsr_diagnostics // mixed; empty unless any *_filter_mode == 'vqsr'
}
