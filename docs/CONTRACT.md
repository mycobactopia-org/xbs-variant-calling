# XBS_VARIANT_CALLING — Interface Contract

This document is the **stable interface contract** for the `XBS_VARIANT_CALLING` subworkflow. Consumer pipelines bind to this contract; renames or signature changes here are breaking changes for every consumer.

**Consumers (confirmed):**

- **MAGMA** (`mycobactopia-org/MAGMA`, branch `nf-core-port`) — drives the Heupink 2021 spine for *M. tuberculosis*, replacing MAGMA's `CALL_WF` per-sample chain and `MERGE_WF` CombineGVCFs→VQSR chain (Plan-3 §11).
- **mtbc-varcaller-nf** (`mycobactopia-org/mtbc-varcaller-nf`) — wraps XBS as the `gatk_vqsr` backend in a multi-backend MTBC variant-calling substrate. Spec: `abc-universe/specs/active/mtbc-varcaller-nf.md`.

**Stability commitment:** the `take:` and `emit:` shapes below are versioned alongside the subworkflow; any breaking change requires a major version bump in `nextflow.config` `manifest.version` and a CHANGELOG entry.

---

## `take:` channels (subworkflow inputs)

```
workflow XBS_VARIANT_CALLING {
    take:
    ch_reads        // channel  [ meta(study, sample, library, id), [r1, r2] ]
    ch_reference    // value    [ meta(id:'ref'),       fasta, fai, dict, [amb, ann, bwt, pac, sa] ]
    ch_snp_truth    // value    [ meta(id:'snp_truth'),   vcf, tbi ]   — empty list-elements OK when not needed
    ch_indel_truth  // value    [ meta(id:'indel_truth'), vcf, tbi ]   — empty list-elements OK when not needed
    ch_dbsnp        // value    [ meta(id:'dbsnp'),       vcf, tbi ]   — empty list-elements OK when not needed
    ...
}
```

### When channels are required vs may be empty

| Channel | Required when |
|---|---|
| `ch_reads` | always |
| `ch_reference` | always |
| `ch_snp_truth` | `params.snp_filter_mode == 'vqsr'` |
| `ch_indel_truth` | `params.indel_filter_mode == 'vqsr'` |
| `ch_dbsnp` | `params.skip_bqsr == false` |

For non-required cases, consumers pass `channel.value([ [id:'<name>'], [], [] ])` (empty path slots), matching the convention nf-core modules use for optional inputs.

## `emit:` channels (subworkflow outputs)

### Per-sample outputs (always emitted)

| Emit | Shape | File | Notes |
|---|---|---|---|
| `gvcfs` | `[ meta, vcf, tbi ]` | `*.g.vcf.gz` + `.tbi` | GATK HaplotypeCaller GVCF, per row in input samplesheet |
| `sample_bam` | `[ meta, bam, bai ]` | `*.markdup.bam` (or `*.recal.bam` when `skip_bqsr=false`) | Post-MarkDup or post-BQSR BAM. Consumers like MAGMA use this for downstream LoFreq / DELLY |
| `samtools_stats` | `[ meta, *.stats ]` | `*.SamtoolStats` | `samtools stats -F DUP,SUPPLEMENTARY,SECONDARY,UNMAP,QCFAIL` |
| `flagstat` | `[ meta, *.flagstat ]` | `*.FlagStat` | `samtools flagstat` output |
| `wgs_metrics` | `[ meta, *_metrics ]` | `*.WgsMetrics` | Picard CollectWgsMetrics |
| `markdup_metrics` | `[ meta, *.metrics ]` | `*.metrics` | GATK MarkDuplicates metrics |

### Cohort outputs (empty channels when `params.skip_cohort=true`)

| Emit | Shape | File | Notes |
|---|---|---|---|
| `raw_variants` | `[ meta, vcf, tbi ]` | `${joint_name}.raw_variants.vcf.gz` | Post-GenotypeGVCFs, unfiltered |
| `snp_filtered` | `[ meta, vcf, tbi ]` | `${joint_name}.{FilteredSNPs,FilteredSNPs.hardfilter,raw_snps}.vcf.gz` | Contents depend on `snp_filter_mode` (see below) |
| `indel_filtered` | `[ meta, vcf, tbi ]` | `${joint_name}.{FilteredINDELs,FilteredINDELs.hardfilter,raw_indels}.vcf.gz` | Contents depend on `indel_filter_mode` |
| `vqsr_diagnostics` | mixed | `*.tranches`, `*.recal.vcf.gz`, `*plots.R` | Only emitted when any `*_filter_mode == 'vqsr'`; empty otherwise |

### Versions

Versions are emitted via the standard nf-core `versions` topic (channel-via-topic). Consumers should `channel.topic("versions")` to collect them; do not consume a `versions` emit directly from this subworkflow.

---

## Filter-mode contract

`snp_filter_mode` and `indel_filter_mode` are tri-state strings driving the cohort filter stage independently per variant class. The contract holds the **emit channel shape constant** across modes — only the file contents change.

| Mode | What runs | Truth set needed | File suffix |
|---|---|---|---|
| `vqsr` (default) | `GATK4_VARIANTRECALIBRATOR_*` + `GATK4_APPLYVQSR_*` | yes (the `ch_{snp,indel}_truth` channel) | `.FilteredSNPs.vcf.gz` / `.FilteredINDELs.vcf.gz` |
| `hard_filters` | `GATK4_VARIANTFILTRATION_*` with the configurable filter expression | no | `.FilteredSNPs.hardfilter.vcf.gz` / `.FilteredINDELs.hardfilter.vcf.gz` |
| `none` | none — pass-through from `GATK4_SELECTVARIANTS_*` | no | `.raw_snps.vcf.gz` / `.raw_indels.vcf.gz` |

**MAGMA mapping (after Plan-3 §11 integration):**
- `snp_filter_mode = 'vqsr'` (MAGMA always ran VQSR on SNPs)
- `indel_filter_mode = 'hard_filters'` (MAGMA always used hard filters for INDELs)

**Per-organism / per-cohort tuning:**
- `snp_hard_filter_expression` / `snp_hard_filter_name` — defaults are GATK best-practice tuned for haploid bacterial calling
- `indel_hard_filter_expression` / `indel_hard_filter_name` — same
- `vqsr_max_gaussians_{snp,indel}` — knock down for small cohorts (Heupink 2021 used 12+ samples; for N=3 set to 1)
- `target_titv` — organism-specific (1.85 for MTB; leave null otherwise)

---

## `skip_cohort` contract

When `params.skip_cohort = true`:

- `XBS_PER_SAMPLE` runs end-to-end (mapping → MarkDup → optional BQSR → HC GVCF + QC)
- `XBS_COHORT` is skipped entirely
- All cohort emits (`raw_variants`, `snp_filtered`, `indel_filtered`, `vqsr_diagnostics`) are `channel.empty()`
- Consumers can wire to those emits unconditionally — they will simply receive no items

This supports incremental cohort patterns (run per-sample first, defer joint calling) and consumers that want to swap in a different cohort-stage backend (e.g. mtbc-varcaller-nf consensus).

---

## Version semantics

| Change type | Example | Version bump |
|---|---|---|
| Add new `emit:` channel | new per-sample QC metric exposed | minor (no consumer breaks) |
| Add new `take:` parameter with sensible default | new optional truth set | minor |
| Rename an `emit:` channel | `gvcfs` → `per_sample_gvcfs` | **major** (consumers break) |
| Change `emit:` channel shape | `[meta, vcf, tbi]` → `[meta, vcf]` | **major** |
| Change `take:` parameter required-ness | making `ch_dbsnp` always required | **major** |
| Change default of a `params.*` value | `snp_filter_mode = 'vqsr'` → `'hard_filters'` | **major** (changes scientific outputs) |
| Add new behavior behind a flag (default = current behavior) | new optional region-exclusion step | minor |

The subworkflow's pipeline-level version (`nextflow.config` `manifest.version`) is the contract version. Pin against a tagged commit, not against `master`.

---

## Audit checklist for breaking changes

Before merging any PR that touches `subworkflows/local/xbs_variant_calling_wf.nf`'s `take:` or `emit:` blocks:

- [ ] Are channel names unchanged?
- [ ] Are channel shapes unchanged (number + type of tuple elements)?
- [ ] Is the required-vs-empty contract for each `take:` channel unchanged?
- [ ] Does `tests/conformance/compare_command_sh.sh` still pass?
- [ ] If anything above is "no" → bump `manifest.version` major; CHANGELOG; notify MAGMA + mtbc-varcaller-nf maintainers.
