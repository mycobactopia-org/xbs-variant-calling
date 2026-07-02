/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XBS_PER_SAMPLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per-sample portion of the XBS spine (Heupink 2021). Maps each library,
    merges libraries within a sample, marks duplicates, optionally runs
    BQSR (skipped by default — paper-true), runs HaplotypeCaller in GVCF
    mode, and emits per-sample QC stats.

    Stages 1–7 from the canonical XBS-core bash:
      1. BWA-MEM (sort-on-the-fly via the nf-core module's sort_bam=true)
      2. samtools index
      3. samtools merge across libraries within a sample
      4. GATK MarkDuplicates
      5. (Optional) GATK BaseRecalibrator + ApplyBQSR
      6. samtools stats + flagstat + Picard CollectWgsMetrics
      7. GATK HaplotypeCaller (ploidy=1, -ERC GVCF)
*/

include { BWA_MEM                                              } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_LIB                 } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                                       } from '../../modules/nf-core/samtools/merge/main'
include { GATK4_MARKDUPLICATES                                 } from '../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP             } from '../../modules/nf-core/samtools/index/main'
include { GATK4_BASERECALIBRATOR                               } from '../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR                                      } from '../../modules/nf-core/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RECAL               } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                                       } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT                                    } from '../../modules/nf-core/samtools/flagstat/main'
include { PICARD_COLLECTWGSMETRICS                             } from '../../modules/nf-core/picard/collectwgsmetrics/main'
include { GATK4_HAPLOTYPECALLER                                } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4SPARK_HAPLOTYPECALLER                           } from '../../modules/local/gatk4spark/haplotypecaller/main'


workflow XBS_PER_SAMPLE {

    take:
    ch_reads      // channel: [ meta(study, sample, library, id), [r1, r2] ]
    ch_reference  // value:   [ meta(id:'ref'), fasta, fai, dict, [amb, ann, bwt, pac, sa] ]
    ch_dbsnp      // value:   [ meta, vcf, tbi ]  (empty unless !params.skip_bqsr)

    main:

    // ---- shape adapters for the various reference channel signatures ----
    def ch_fasta_tuple     = ch_reference.map { m, f, _fai, _d, _b -> [m, f] }
    def ch_fai_tuple       = ch_reference.map { m, _f, fai, _d, _b -> [m, fai] }
    def ch_dict_tuple      = ch_reference.map { m, _f, _fai, d, _b -> [m, d] }
    def ch_fasta_fai_tuple = ch_reference.map { m, f, fai, _d, _b -> [m, f, fai] }
    def ch_bwa_index_tuple = ch_reference.map { m, _f, _fai, _d, bwa -> [m, bwa] }
    def ch_fasta_path      = ch_reference.map { _m, f, _fai, _d, _b -> f }
    def ch_fai_path        = ch_reference.map { _m, _f, fai, _d, _b -> fai }
    def ch_dict_path       = ch_reference.map { _m, _f, _fai, d, _b -> d }

    def ch_dbsnp_vcf_tuple = ch_dbsnp.map { m, vcf, _tbi -> [m, vcf] }
    def ch_dbsnp_tbi_tuple = ch_dbsnp.map { m, _vcf, tbi -> [m, tbi] }

    //
    // STAGE 1: BWA-MEM (with sort_bam=true → sorted BAM out-of-the-box)
    //
    BWA_MEM(ch_reads, ch_bwa_index_tuple, ch_fasta_tuple, true)

    //
    // STAGE 2: index per-library BAM
    //
    SAMTOOLS_INDEX_LIB(BWA_MEM.out.bam)

    //
    // STAGE 3: collapse libraries per sample → SAMTOOLS_MERGE
    //          Even single-library samples go through MERGE for shape uniformity.
    //
    def ch_per_sample_bams = BWA_MEM.out.bam
        .join(SAMTOOLS_INDEX_LIB.out.index)
        .map { meta, bam, bai ->
            // promote `sample` → meta.id so groupTuple collapses libraries
            def new_meta = meta + [ id: meta.sample, library: null ]
            [ new_meta, bam, bai ]
        }
        .groupTuple(by: 0)

    def ch_merge_input = ch_per_sample_bams.map { meta, bams, bais -> [meta, bams, bais] }
    SAMTOOLS_MERGE(ch_merge_input, channel.value([ [id: 'none'], [], [], [] ]))

    //
    // STAGE 4: Mark duplicates
    //
    GATK4_MARKDUPLICATES(SAMTOOLS_MERGE.out.bam, [], [])

    //
    // STAGE 5: (optional) BQSR — skipped by default (paper-true)
    //
    def ch_call_bam
    def ch_call_bai
    if (!params.skip_bqsr) {

        SAMTOOLS_INDEX_MARKDUP(GATK4_MARKDUPLICATES.out.bam)

        def ch_brsr_input = GATK4_MARKDUPLICATES.out.bam
            .join(SAMTOOLS_INDEX_MARKDUP.out.index)
            .map { meta, bam, bai -> [meta, bam, bai, []] }     // empty intervals
        GATK4_BASERECALIBRATOR(
            ch_brsr_input,
            ch_fasta_tuple, ch_fai_tuple, ch_dict_tuple,
            ch_dbsnp_vcf_tuple, ch_dbsnp_tbi_tuple
        )

        def ch_applybqsr_input = GATK4_MARKDUPLICATES.out.bam
            .join(SAMTOOLS_INDEX_MARKDUP.out.index)
            .join(GATK4_BASERECALIBRATOR.out.table)
            .map { meta, bam, bai, table -> [meta, bam, bai, table, []] }    // empty intervals
        GATK4_APPLYBQSR(ch_applybqsr_input, ch_fasta_path, ch_fai_path, ch_dict_path)

        SAMTOOLS_INDEX_RECAL(GATK4_APPLYBQSR.out.bam)
        ch_call_bam = GATK4_APPLYBQSR.out.bam
        ch_call_bai = SAMTOOLS_INDEX_RECAL.out.index
    } else {
        SAMTOOLS_INDEX_MARKDUP(GATK4_MARKDUPLICATES.out.bam)
        ch_call_bam = GATK4_MARKDUPLICATES.out.bam
        ch_call_bai = SAMTOOLS_INDEX_MARKDUP.out.index
    }

    def ch_bam_bai = ch_call_bam.join(ch_call_bai)

    //
    // STAGE 6: per-sample QC stats
    //
    SAMTOOLS_STATS(ch_bam_bai, ch_fasta_fai_tuple)
    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    PICARD_COLLECTWGSMETRICS(ch_bam_bai, ch_fasta_tuple, ch_fai_tuple, [])

    //
    // STAGE 7: HaplotypeCaller (per-sample GVCF, ploidy=1 via ext.args)
    //
    // Two interchangeable backends, identical vcf/tbi emit shape:
    //   - GATK4_HAPLOTYPECALLER       (default; non-Spark; nf-core)
    //   - GATK4SPARK_HAPLOTYPECALLER  (opt-in via params.use_spark_haplotypecaller)
    //
    // Versions are collected via nf-core's topic-based `versions` channel — no
    // per-branch versions capture is needed here (nf-core's GATK4_HAPLOTYPECALLER
    // does not emit a `.out.versions` sub-channel; it publishes to `topic: versions`
    // as `versions_gatk4`, which the top-level workflow harvests separately).
    def ch_hc_input = ch_bam_bai.map { meta, bam, bai -> [meta, bam, bai, [], []] }   // empty intervals + dragstr_model
    def ch_hc_vcf
    def ch_hc_tbi
    if (params.use_spark_haplotypecaller) {
        GATK4SPARK_HAPLOTYPECALLER(
            ch_hc_input,
            ch_fasta_tuple, ch_fai_tuple, ch_dict_tuple,
            channel.value([[id:'none'], []]),
            channel.value([[id:'none'], []])
        )
        ch_hc_vcf = GATK4SPARK_HAPLOTYPECALLER.out.vcf
        ch_hc_tbi = GATK4SPARK_HAPLOTYPECALLER.out.tbi
    } else {
        GATK4_HAPLOTYPECALLER(
            ch_hc_input,
            ch_fasta_tuple, ch_fai_tuple, ch_dict_tuple,
            channel.value([[id:'none'], []]),   // dbsnp not used in HC (BQSR-only)
            channel.value([[id:'none'], []])
        )
        ch_hc_vcf = GATK4_HAPLOTYPECALLER.out.vcf
        ch_hc_tbi = GATK4_HAPLOTYPECALLER.out.tbi
    }

    emit:
    gvcf_vcf        = ch_hc_vcf                            // [ meta, *.g.vcf.gz ]
    gvcf_tbi        = ch_hc_tbi                            // [ meta, *.g.vcf.gz.tbi ]
    sample_bam      = ch_bam_bai                          // [ meta, bam, bai ]
    samtools_stats  = SAMTOOLS_STATS.out.stats            // [ meta, *.stats ]
    flagstat        = SAMTOOLS_FLAGSTAT.out.flagstat      // [ meta, *.flagstat ]
    wgs_metrics     = PICARD_COLLECTWGSMETRICS.out.metrics // [ meta, *_metrics ]
    markdup_metrics = GATK4_MARKDUPLICATES.out.metrics    // [ meta, *.metrics ]
}
