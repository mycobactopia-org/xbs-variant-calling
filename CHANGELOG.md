# mycobactopia-org/xbs-variant-calling: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.4.0dev - unreleased

GATK optimisations — Pattern 1 from the GATK-optimisations spec ([`abc-universe/specs/active/xbs-variant-calling-gatk-optimizations.md`](../../abc-universe/specs/active/xbs-variant-calling-gatk-optimizations.md)).

### `Added`

- **Optional `HaplotypeCallerSpark` backend** (param `skip_gatk4_haplotypecaller_spark`, default `true`). When `false`, the per-sample stage routes BAMs through a new local module `GATK4SPARK_HAPLOTYPECALLER` that runs `gatk HaplotypeCallerSpark --spark-master local[$task.cpus]` for per-sample parallelism. The Spark and non-Spark modules emit identical `{vcf, tbi, vcf_tbi, versions}` channels so every downstream stage is agnostic. Default `true` keeps the per-sample output byte-equivalent to `v0.3.0` for any current invocation.
- New local module `modules/local/gatk4spark/haplotypecaller/main.nf` (container `quay.io/biocontainers/gatk4:4.6.0.0--py310hdfd78af_0` / singularity equivalent — same GATK4 build that nf-core's non-Spark module uses).

### `Pending`

- Pattern 2 — GenomicsDB joint genotyping (per-chromosome scatter) — separate PR.
- Pattern 3 — gVCF checkpoints with `zstd` / `genozip` codecs — separate PRs.

## v0.3.0 — released 2026-06-30

Stable interface contract for consumption by MAGMA and mtbc-varcaller-nf.

### `Added`

- **`docs/CONTRACT.md`** — explicit interface contract for `XBS_VARIANT_CALLING` subworkflow: `take:` / `emit:` shapes, required-vs-empty rules per channel, filter-mode contract, `skip_cohort` semantics, version semantics for breaking changes, audit checklist.
- **`skip_cohort` parameter** (default `false`) — when true, runs only `XBS_PER_SAMPLE` (per-sample mapping + HaplotypeCaller GVCF). Cohort emits become empty channels but stay safely wired. Supports incremental cohort patterns and consumers that swap in a different cohort-stage backend.
- **`{snp,indel}_filter_mode` parameters** (`'vqsr' | 'hard_filters' | 'none'`, default `'vqsr'`) — tri-state replacement for the historical `skip_{snp,indel}_vqsr` booleans. Adds a proper GATK VariantFiltration path so consumers (notably MAGMA, which always used hard filters for INDELs) get correct filtering rather than a pass-through.
- **Hard-filter expressions** as configurable params (`{snp,indel}_hard_filter_expression`, `{snp,indel}_hard_filter_name`) with GATK best-practice defaults tuned for haploid bacterial calling.
- **`gatk4/variantfiltration`** nf-core module installed; aliased as `GATK4_VARIANTFILTRATION_SNP` and `GATK4_VARIANTFILTRATION_INDEL` for the hard-filter path.
- **`nextflow_schema.json`** sections: `reference_bundle_options`, `variant_filtering_options`, `bqsr_options`, `cohort_options`, `organism_options` — every XBS-specific param now has schema validation and `--help` discoverability.
- **`tests/conformance/compare_command_sh.sh`** extended to cover the two new VariantFiltration aliases — now 15/15 modules pass.

### `Changed`

- **`subworkflows/local/xbs_cohort.nf`** — replaced the per-mode `if (!params.skip_*_vqsr)` branches with a 3-way `switch` on `*_filter_mode` (`vqsr` / `hard_filters` / `none`).
- **`subworkflows/local/xbs_variant_calling_wf.nf`** — wraps `XBS_COHORT` in `if (!params.skip_cohort)` so cohort emits are empty channels (not absent) when skipped. Header comments updated to flag the stable contract.
- **`workflows/xbs-variant-calling.nf`** — truth-set file existence checks now branch on `*_filter_mode == 'vqsr'` (not on `skip_*_vqsr`).
- **`conf/test.config`** — switched to the new flags: `snp_filter_mode='vqsr'` + `indel_filter_mode='hard_filters'` (was `skip_indel_vqsr=true`). MAGMA-parity behaviour now exercised end-to-end by the test profile.

### `Removed`

- `skip_snp_vqsr` and `skip_indel_vqsr` parameters — superseded by `{snp,indel}_filter_mode`. No deprecation shim since the project has had only `v0.x` releases.

## v0.2.0-phase3 - 2026-06-12

Runnable Phase-3 setup with bundled MTB test data + abc-cluster profile. Awaiting abc-cluster end-to-end submission.

### `Added`

- Bundled `resources/genome/NC-000962-3-H37Rv.*` (~9 MB) for the MTB test profile.
- Bundled `resources/truth/Coll2018.UVPapproved.rRNAexcluded.vcf.gz` as the SNP VQSR truth set.
- `assets/samplesheet_test.csv` — 3 EXIT-RIF / PRJNA1026351 samples (matches MAGMA's `test` profile).
- `conf/test.config` end-to-end MTB test profile.
- `conf/abc_cluster.config` — abc-cluster nf-nomad executor with MAGMA-parity container choices.

## v0.1.0dev - 2026-06-12

Initial scaffold from the nf-core template + Phases 0–2 build.

### `Added`

- nf-core v4.0.2 scaffold via `nf-core pipelines create`.
- 16 nf-core registry modules (bwa/mem, samtools/{sort,index,merge,stats,flagstat}, gatk4/{markduplicates,baserecalibrator,applybqsr,haplotypecaller,combinegvcfs,genotypegvcfs,selectvariants,variantrecalibrator,applyvqsr}, picard/collectwgsmetrics).
- Two patched nf-core modules: `gatk4/combinegvcfs` adds `tbi` emit; `gatk4/selectvariants` adds `ext.intervals_mode='exclude'` knob.
- `subworkflows/local/xbs_per_sample.nf` — Heupink 2021 stages 1–7.
- `subworkflows/local/xbs_cohort.nf` — Heupink 2021 stages 8–12.
- `subworkflows/local/xbs_variant_calling_wf.nf` — top-level composer.
- `conf/modules.config` — every `ext.args` reproducing Heupink 2021 flags byte-for-byte.
- `docs/VALIDATION.md` — full MVP test strategy + deferred-benchmark plan.
- `tests/conformance/compare_command_sh.sh` — Layer-0 conformance against the XBS-core bash flags (13/13 pass at this version).
- `.github/workflows/conformance.yml` — CI wiring for the conformance script.
