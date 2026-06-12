# mycobactopia-org/xbs-variant-calling

[![GitHub Actions CI Status](https://github.com/mycobactopia-org/xbs-variant-calling/actions/workflows/nf-test.yml/badge.svg)](https://github.com/mycobactopia-org/xbs-variant-calling/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/mycobactopia-org/xbs-variant-calling/actions/workflows/linting.yml/badge.svg)](https://github.com/mycobactopia-org/xbs-variant-calling/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**mycobactopia-org/xbs-variant-calling** is a Nextflow re-implementation of the **XBS** (com**pleX** **B**acterial **S**ample) variant-calling pipeline described in [Goig et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743552/), faithful to the canonical bash implementation at [TORCH-Consortium/XBS-variant-calling-core](https://github.com/TORCH-Consortium/XBS-variant-calling-core), packaged as a reusable nf-core-style subworkflow for **any haploid bacterium**.

Given paired-end short reads, a reference genome, and SNP + INDEL truth-set VCFs, it produces:

- per-sample GVCFs (GATK HaplotypeCaller, `-ploidy 1 -ERC GVCF`)
- per-sample QC stats (samtools stats, GATK CollectWgsMetrics, GATK FlagStat)
- cohort-level joint-genotyped VCFs (GATK CombineGVCFs + GenotypeGVCFs)
- filtered SNP and INDEL VCFs (GATK VQSR, separate SNP + INDEL models)
- VQSR diagnostics (tranches, recal tables, R plots)

The pipeline was designed for **low-coverage** and **contaminated** sequencing data — Goig 2022 reports unaffected performance at 5–10× depth and across contamination levels up to >99.99%, where standard hard-filtering pipelines lose sensitivity. The combination of joint calling + VQSR is what enables this.

### Pipeline stages

1. **Map** — BWA-MEM (`-M`, configurable extra args)
2. **Sort + index** — samtools
3. **Library merge** — samtools merge (when a sample has multiple libraries)
4. **Mark duplicates** — GATK MarkDuplicates
5. **(Optional) BQSR** — GATK BaseRecalibrator + ApplyBQSR. **Skipped by default** — Goig 2022 explicitly omits this to avoid contaminant DNA variants being interpreted as systematic errors.
6. **Per-sample QC** — samtools stats, GATK CollectWgsMetrics, GATK FlagStat (emitted only; no gating)
7. **Per-sample HaplotypeCaller** — `-ploidy 1 -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation --read-filter MappingQualityNotZeroReadFilter`
8. **Cohort CombineGVCFs**
9. **Cohort GenotypeGVCFs** — `--sample-ploidy 1`
10. **Split SNP / INDEL** — GATK SelectVariants
11. **SNP VQSR** — VariantRecalibrator + ApplyVQSR (annotations: `AS_MQRankSum, AS_QD, AS_MQ, DP`; SNP truth set required; `--max-gaussians 4`; tranche 99.9; configurable `--target-titv`)
12. **INDEL VQSR** — VariantRecalibrator + ApplyVQSR (annotations: `MQRankSum, QD, DP`; INDEL truth set required; `--max-gaussians 2`; `--lod-score-cutoff 0.0`)

### Genericness contract

The pipeline is *organism-agnostic by configuration*. Only five values change per organism:

| Parameter | M. tuberculosis (default test profile) | What other organisms supply |
|---|---|---|
| `reference` bundle | NC_000962.3 H37Rv | their reference + BWA index |
| `snp_truth_vcf` | Coll2018 | a high-confidence species-level callset |
| `indel_truth_vcf` | curated set | their own |
| `target_titv` | 1.85 | organism-specific, or `null` to disable |
| `dbsnp_vcf` (only if `--skip_bqsr false`) | Benavente2015 | their own known-sites VCF |

Ploidy stays at 1 for all bacteria.

### Use as a subworkflow

The intended consumption path is via the `XBS_VARIANT_CALLING` subworkflow under `subworkflows/local/`, installable into another pipeline via `nf-core subworkflows install` once published to nf-core/subworkflows. The standalone `main.nf` is provided so the pipeline can be tested, validated, and benchmarked independently of any downstream consumer.

Primary downstream consumer: [mycobactopia-org/MAGMA](https://github.com/mycobactopia-org/MAGMA) (M. tuberculosis pipeline) — once Phase 4 integration lands, MAGMA's per-sample HC + cohort joint calling + VQSR will be delegated to this subworkflow.

## Usage

> [!NOTE]
> If you are new to Nextflow, please refer to [the nf-core docs](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set up Nextflow. Test your setup with `-profile test` before running on real data.

Prepare a samplesheet (`assets/samplesheet.csv` is the template):

```csv
study,sample,library,r1,r2
test_study,SAMPLE_A,lib1,/path/to/SAMPLE_A_lib1_R1.fastq.gz,/path/to/SAMPLE_A_lib1_R2.fastq.gz
test_study,SAMPLE_A,lib2,/path/to/SAMPLE_A_lib2_R1.fastq.gz,/path/to/SAMPLE_A_lib2_R2.fastq.gz
test_study,SAMPLE_B,lib1,/path/to/SAMPLE_B_R1.fastq.gz,/path/to/SAMPLE_B_R2.fastq.gz
```

Multiple `library` rows for the same `sample` are merged after mapping.

Run:

```bash
nextflow run mycobactopia-org/xbs-variant-calling \
   -profile <docker|singularity|conda>,test_mtb \
   --input samplesheet.csv \
   --outdir results/
```

> [!WARNING]
> Pass pipeline parameters via the CLI or `-params-file`. Custom config files via `-c` should not override parameters — see [the nf-core docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

mycobactopia-org/xbs-variant-calling was developed by Abhinav Sharma.

The XBS pipeline this implementation is based on was developed by Tim H.H. Coorens and collaborators in the lab of Lennard Epping at the Bavarian Health and Food Safety Authority — see [Goig et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743552/) and the original bash scripts at [TORCH-Consortium/XBS-variant-calling-core](https://github.com/TORCH-Consortium/XBS-variant-calling-core).

## Contributions and Support

If you would like to contribute, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

If you use this pipeline, please cite **Goig et al. 2022** ([PMC8743552](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743552/)) for the XBS methodology, and consult [`CITATIONS.md`](CITATIONS.md) for full citations of every tool the pipeline uses.

This pipeline reuses scaffolding from the [nf-core](https://nf-co.re) community framework under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE):

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> *Nat Biotechnol.* 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
