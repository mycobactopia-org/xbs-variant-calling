# XBS Validation & Test Strategy

This document describes how `xbs-variant-calling` is validated against (a) its own runtime correctness, (b) the Heupink 2021 XBS specification, and (c) the broader benchmark landscape including MAGMA and MTBseq-nf. It is deliberately compute-conservative — the **MVP test suite runs on a single laptop in minutes**, while heavyweight publication-grade benchmarks are documented but deferred.

Reference: [Heupink et al. 2021, *Microb Genom* 7(11):000689](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743552/).
Canonical bash implementation: [TORCH-Consortium/XBS-variant-calling-core](https://github.com/TORCH-Consortium/XBS-variant-calling-core).

---

## 1. What XBS emits

The subworkflow is intentionally minimal — every output is "minimal sufficient" for a downstream MTBC analysis layer (MAGMA, TBANALYZER, M. bovis pipelines) to bolt onto.

### Per-sample outputs

| Output | Channel | File pattern | Use |
|---|---|---|---|
| GVCF | `gvcfs` | `*.g.vcf.gz` + `.tbi` | Cohort joint-genotyping input |
| Sample BAM | `sample_bam` | `*.markdup.bam` (or `*.recal.bam` if BQSR) + `.bai` | Consumers' downstream chains (LoFreq, DELLY, etc.) |
| MarkDup metrics | `markdup_metrics` | `*.metrics` | Duplicate rate |
| samtools stats | `samtools_stats` | `*.SamtoolStats` | Read count, mapping rate, error rate |
| samtools flagstat | `flagstat` | `*.FlagStat` | Pair/dup/QC flag counts |
| Picard WgsMetrics | `wgs_metrics` | `*.WgsMetrics` | Coverage depth + breadth |

### Cohort outputs

| Output | Channel | File pattern | Stage |
|---|---|---|---|
| Raw cohort variants | `raw_variants` | `${joint_name}.raw_variants.vcf.gz` + `.tbi` | post-GenotypeGVCFs |
| Filtered SNPs | `snp_filtered` | `${joint_name}.FilteredSNPs.vcf.gz` + `.tbi` | post-ApplyVQSR SNP (or pass-through SelectVariants if `skip_snp_vqsr`) |
| Filtered INDELs | `indel_filtered` | `${joint_name}.FilteredINDELs.vcf.gz` + `.tbi` | post-ApplyVQSR INDEL (or pass-through if `skip_indel_vqsr`) |
| VQSR diagnostics | `vqsr_diagnostics` | `*.tranches`, `*.recal.vcf.gz`, `*plots.R` (per mode) | Inspecting / publishing tranche curves |

### Explicitly not emitted by XBS — left to consumers

LoFreq minor variant VCFs, DELLY structural variant VCFs, TBprofiler resistance JSON/itol summaries, phylogeny trees, cluster picks, SnpEff annotation, lineage assignment, region-excluded VCFs (rRNA / PE-PPE), cohort QC scoring / multi-infection filter. These all sit *on top of* XBS's emits.

---

## 2. MVP test suite (runnable on a laptop, gated by CI)

Optimized for "**catch every class of failure with the fewest CPU-hours**" rather than statistical power.

### Failure classes vs minimum test

| Failure class | Smallest test that catches it | Layer |
|---|---|---|
| Process-level errors (signature mismatch, ext.args typo, channel shape) | `nextflow lint` + `nextflow inspect` | **L0** |
| Per-sample logic (multi-library merge, optional BQSR, HC params) | 3 samples with mixed library counts | **L1** |
| Cohort-level logic (CombineGVCFs + GenotypeGVCFs + VQSR convergence on small N) | 3 samples with `vqsr_max_gaussians_snp=1` | **L1** |
| Heupink 2021 spec conformance | `.command.sh` text-diff against XBS-core bash | **L0** |
| Determinism / non-deterministic ordering bugs | Same input twice → same SHA256 of outputs | **L2** |
| Cross-implementation regression vs MAGMA | XBS-standalone outputs vs MAGMA's equivalent stages on identical samples | **L3** |

### Layer 0 — Static analysis (cost ≈ 0; runs on every push)

| Test | What it catches | Mechanism |
|---|---|---|
| `nextflow lint .` | Syntax, deprecated patterns, identifier rules | GitHub Actions `linting.yml` |
| `nf-core lint` | nf-core conformance, schema validity | GitHub Actions `linting.yml` |
| `nextflow inspect -profile test main.nf` | Workflow composition, missing emits | Add to CI |
| `tests/conformance/compare_command_sh.sh` | Heupink 2021 flag parity at module level — text-diffs our emitted `.command.sh` against the XBS-core bash scripts | Add to CI |

The conformance script is the most under-leveraged test in the suite — it answers "does our `.command.sh` produce the same flags as the canonical XBS bash?" in milliseconds and **without running any analysis**.

### Layer 1 — Small-data end-to-end (cost ≈ 5 min local, ≈ 10 min on abc-cluster)

The current `-profile test` (3 EXIT-RIF samples) exercises:

- BWA-MEM with `sort_bam=true`
- Single-library `SAMTOOLS_MERGE` path
- `GATK4_MARKDUPLICATES`
- BQSR-skip path (`skip_bqsr=true`)
- Post-dedup `SAMTOOLS_INDEX_MARKDUP`
- Per-sample QC (samtools stats + flagstat + Picard wgsmetrics)
- `GATK4_HAPLOTYPECALLER` per sample
- Cohort: `CombineGVCFs` + `GenotypeGVCFs` + `SelectVariants` (SNP+INDEL)
- SNP VQSR with `--max-gaussians 1` (convergence on small N)
- INDEL VQSR skip path

Wired as `tests/default.nf.test` with output-existence assertions and record-count sanity checks. CI runs it on every PR.

### Layer 2 — Reproducibility (cost ≈ 2× L1; runs weekly on CI)

```bash
nextflow run main.nf -profile test --outdir run1
nextflow run main.nf -profile test --outdir run2

# Diff cohort SNP VCF body (ignore header timestamps)
diff <(zcat run1/.../test_cohort.FilteredSNPs.vcf.gz | grep -v ^#) \
     <(zcat run2/.../test_cohort.FilteredSNPs.vcf.gz | grep -v ^#)
```

Catches non-deterministic ordering bugs. Cheap to schedule on a weekly cron rather than per-PR.

### Layer 3 — Cross-implementation regression vs MAGMA (cost ≈ L1 + manual)

The scientific gate before declaring Phase 3 "validated":

```
Run XBS standalone on 3-sample EXIT-RIF set     → outputs A
Pull MAGMA v0.3.0 outputs from SciVer run        → outputs B (already exist)

Compare for the equivalent stages:
   A:test_cohort.raw_variants.vcf.gz   ↔   B:joint.filtered_SNP.RawIndels.vcf.gz
   A:test_cohort.FilteredSNPs.vcf.gz   ↔   B:joint.filtered_SNP_*-rRNA.vcf.gz  (after re-adding rRNA exclusion)
   A:per-sample GVCFs                  ↔   B:per-sample GVCFs
   A:per-sample QC                     ↔   B:per-sample QC
```

If record counts agree within ±5% on a 3-sample run, XBS produces scientifically equivalent output to the MAGMA-embedded GATK chain. Script lives at `tests/regression/compare_vs_magma.sh` (drafted; runs once we have a Layer-1 pass).

---

## 3. Deferred — what NOT to do in the MVP

Documented here so the bigger ambitions don't get forgotten, but explicitly out-of-scope for current validation. Implementation gated on (a) abc-cluster connectivity returning, (b) compute/storage budget becoming available, (c) Phase 4 MAGMA integration landing.

### 3.1 Tiered cohort ladder (mtbseq-nf-style scaling validation)

Run XBS at increasing N and verify statistical convergence + compute scaling holds.

| Tier | N | Approx wall-clock on abc-cluster | What it validates |
|---|---|---|---|
| T1 | 30 | ~1 h | Smoke test, all subworkflows fire |
| T2 | 200 | ~6–8 h | First statistically meaningful F1 CIs; VQSR converges to stable filters |
| T3 | 1000 | 1–3 days | Publication-grade Heupink 2021 replication; first usable resistance & lineage databases |
| T4 | 3000–5000 | 1–2 weeks | Discovery scale: novel resistance markers at MAF ≥ 1%, comprehensive lineage catalog, transmission cluster validation |

Compute cliffs to expect at scale:
- `CombineGVCFs`: O(N) memory; may OOM at N ≥ 500 without interval sharding
- `GenotypeGVCFs`: O(N²) in some operations — hours at N=1000, days at N=5000
- `IQTREE`: O(N² × sites) — days at N=1000 even with `-fast`

Sharding work is a separate Plan-3.5 deliverable if T3 reveals the cliff.

### 3.2 Three-pipeline benchmark (XBS / MAGMA / MTBseq-nf)

Modelled on the MAGMA paper's head-to-head structure with MTBseq, but using **MTBseq-nf** (the Nextflow port) as the comparator since that's what fits our orchestration story.

| Axis | What to compare | Why |
|---|---|---|
| Per-sample SNP set | Jaccard / F1 between pipeline pairs | Concordance, stratified by coverage |
| Per-sample INDEL set | Same | Same |
| Per-sample QC exclusion | Sample-by-sample concordance of "passed cohort QC" sets | MAGMA's claim of fewer false exclusions vs MTBseq |
| DR variant calls | F1 vs phenotypic ground truth per drug | Resistance prediction accuracy |
| Lineage assignment | Confusion matrix vs Coll lineage SNPs | Lineage typing accuracy |
| Cluster membership | Adjusted Rand Index between pipelines' cluster sets | Transmission analysis comparability |
| Wall-clock + cost | From abc-cluster job logs | Practical operational tradeoff |
| Peak memory | From job logs | Scaling cliffs |

Algorithmic distinction worth flagging:

| | XBS / MAGMA | MTBseq / MTBseq-nf |
|---|---|---|
| Variant caller | GATK HaplotypeCaller (local-reassembly) | samtools mpileup + VarScan2 (pileup) |
| Cohort calling | GATK joint genotyping | Per-sample + final merge — no joint genotyping |
| Variant filtering | VQSR (ML) | Hard filters (depth, allele-freq, MAQ) |
| Minor variants | LoFreq (MAGMA add-on) | Same VarScan2 caller at lower AF threshold |
| Lineage typing | TBprofiler | MTBseq's built-in (Coll 2014-based) |
| Clustering | clusterPicker on SNP-distance trees | MTBseq's `classifyTM` |

This is a true three-way **methodological** comparison, not a flag-tweak comparison. Concordance between pipelines isn't the question — the question is which is closer to phenotypic ground truth and which scales further. Suggested home: a separate `mycobactopia-org/magma-benchmark` repo (not this one), so the benchmark machinery lives separately from either pipeline.

### 3.3 Heupink 2021 conformance benchmarks (B1–B3)

Three reproducible benchmarks pulled directly from the paper:

| Benchmark | Setup | Pass criterion |
|---|---|---|
| **B1: Low-coverage F1** | 5 simulated strains × 5× / 10× / 20× / 50× / 100× coverage | F1 ≥ 0.95 at 5× (paper headline) |
| **B2: Contamination resistance** | 1 strain × {0, 10, 25, 50, 90, 99.99}% NTM contamination | F1 ≥ 0.90 at 50% contam |
| **B3: Sputum head-to-head vs MTBseq** | Goig sputum cohort (PRJEB39561) | +13.9% variable sites vs MTBseq baseline (paper claim) |

These belong in `tests/benchmarks/` with `nf-test` snapshot expectations. Compute envelope: each B-test needs ≥ 25 samples, so this is T2-tier work.

### 3.4 MAGMA-modification ablation matrix

MAGMA layers six modifications on top of paper-true XBS — BQSR-on-by-default, LoFreq for minor variants, DELLY for SVs, region exclusion (rRNA, PE/PPE), VQSR grid-search optimization, SnpEff annotation. Each is an additive deviation from Heupink 2021.

The ablation matrix toggles each independently to quantify its contribution:

| Variant | `skip_bqsr` | LoFreq | DELLY | rRNA-excl | VQSR-grid | SnpEff |
|---|---|---|---|---|---|---|
| `xbs_paper_true` | true | off | off | off | off (single-pass) | off |
| `magma_minus_lofreq` | false | **off** | on | on | on | on |
| `magma_minus_delly` | false | on | **off** | on | on | on |
| ... | ... | ... | ... | ... | ... | ... |
| `magma_full` (current default) | false | on | on | on | on | on |

Outputs of each variant on identical T2/T3 samples → per-modification delta tables.

### 3.5 Database-generation pathway

If the long-term goal is **producing databases** (resistance, lineage, phylogeny), the pipeline output needs a stable database-curation handoff.

#### Resistance database extension (TBprofiler-compatible)

- **Input**: per-sample TBprofiler JSON + cohort-wide consolidated DR variant set
- **Curation rule**: variants seen in ≥ 3 phenotypically-resistant samples AND absent in ≥ 10 pan-susceptible samples
- **Output**: drop-in additions to TBprofiler's WHO catalog format
- **Use case**: extending TBprofiler with novel resistance markers discovered in the cohort
- **N required**: ~ 1000 samples (≥ 30 resistant per drug × 15+ drugs)

#### Lineage SNP catalog (Coll-style)

- **Input**: cohort filtered SNP VCF + per-sample lineage assignment
- **Curation rule**: SNPs present in ≥ 95% of one lineage AND ≤ 5% of all others
- **Output**: lineage-defining SNP VCF in Coll 2014 format
- **Use case**: extending Coll 2018 truth set used by VQSR; enabling XBS for other MTBC sub-lineages and other species
- **N required**: ~ 500 samples (≥ 50 per major lineage)

#### Feedback loop

Both curated databases feed back into:
- MAGMA's `resources/truth/Coll2018.UVPapproved.rRNAexcluded.vcf.gz` (next-generation truth set)
- xbs-variant-calling's `snp_truth_vcf` parameter (other-organism generality)
- TBprofiler's WHO catalog (novel resistance markers)

Closing the loop: MAGMA outputs → curation → improved truth sets → better next-generation MAGMA outputs.

---

## 4. Cross-references

- **Phase 3 (standalone XBS validation)**: this document is the validation plan for that phase.
- **Phase 4 (MAGMA integration)**: `mycobactopia-org/MAGMA` `docs/plans/PLAN_3_XBS_EXTRACTION.md §11` — has its own SciVer impact table, file-change list, and resumption checklist.
- **Plan 2 (M. bovis integration)**: `mycobactopia-org/MAGMA` `docs/plans/PLAN_2_MBOVIS_INTEGRATION.md` — Path B routes through XBS once Phase 4 lands.
- **mtbseq-nf benchmark methodology**: future work, separate repo `mycobactopia-org/magma-benchmark`.

---

## 5. Storage budget

The MVP suite is storage-cheap by design — single test run persists < 1 GB to S3 outdir.

| Artifact | Size | Lifecycle |
|---|---|---|
| 3 input FASTQs | ~ 600 MB total (downloaded from EBI SRA on each run; never cached locally) | scratch |
| Per-sample BAMs | ~ 200 MB × 3 | publish to outdir |
| Cohort VCFs | ~ 5 MB total | publish to outdir |
| Logs + trace + `.nextflow/` | ~ 50 MB | scratch |

Layer 0 + Layer 1 fit on a laptop; Layer 2 doubles that; Layer 3 needs MAGMA's existing SciVer outputs as comparator (already on S3, no extra cost).
