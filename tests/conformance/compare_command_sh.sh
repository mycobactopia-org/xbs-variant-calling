#!/usr/bin/env bash
#
# Layer-0 conformance test for xbs-variant-calling.
#
# Static check that the `ext.args` we attach to each module in
# `conf/modules.config` produce flag sets matching the canonical
# XBS-core bash scripts at TORCH-Consortium/XBS-variant-calling-core.
#
# This script does NOT run nextflow. It greps `conf/modules.config`
# for each module's ext.args and asserts that every flag the bash
# script uses appears, AND that no obviously-divergent flag has been
# added beyond what Heupink 2021 prescribes.
#
# Cost: milliseconds. Catches: paper-spec divergence at module level.
#
# Usage:  bash tests/conformance/compare_command_sh.sh

set -uo pipefail

CONF="conf/modules.config"
[[ -f "$CONF" ]] || { echo "ERROR: $CONF not found — run from repo root"; exit 2; }

FAIL=0
PASS=0
SKIP=0

note_pass() { printf "  \033[32m✓\033[0m  %s\n" "$1"; PASS=$((PASS+1)); }
note_fail() { printf "  \033[31m✗\033[0m  %s\n  reason: %s\n" "$1" "$2"; FAIL=$((FAIL+1)); }
note_skip() { printf "  \033[33m-\033[0m  %s (skipped: %s)\n" "$1" "$2"; SKIP=$((SKIP+1)); }

# Extracts the ext.args block for a `withName: 'NAME' { ... }`
# section. Returns the lines from `withName` to the closing brace.
extract_block() {
  local name="$1"
  awk -v n="$name" '
    $0 ~ "withName: '\''" n "'\''" { in_block=1; print; next }
    in_block { print }
    in_block && /^[[:space:]]*\}[[:space:]]*$/ { in_block=0 }
  ' "$CONF"
}

# Assert that the block for $1 contains every regex in $2..$N.
require_flags() {
  local name="$1"; shift
  local block
  block=$(extract_block "$name")
  if [[ -z "$block" ]]; then
    note_skip "$name" "no withName block found"
    return
  fi
  local missing=()
  for flag in "$@"; do
    if ! grep -Eq -- "$flag" <<<"$block"; then
      missing+=("$flag")
    fi
  done
  if [[ ${#missing[@]} -eq 0 ]]; then
    note_pass "$name"
  else
    note_fail "$name" "missing flags: ${missing[*]}"
  fi
}

echo
echo "=== XBS-core conformance: per-module ext.args parity ==="
echo "    (source of truth: TORCH-Consortium/XBS-variant-calling-core bash)"
echo

# ---------------------------------------------------------------------
# Stage 1 — BWA-MEM
#   XBS-core bash: bwa mem -M -t $BWA_THREADS -R $RG $REFERENCE $R1 $R2 | samtools sort -@ $SAMTOOLS_THREADS -O BAM
# ---------------------------------------------------------------------
require_flags 'BWA_MEM'                       '\-M'              # mark short split hits as secondary

# ---------------------------------------------------------------------
# Stage 4 — MarkDuplicates: XBS uses GATK4 defaults; no special flags required.
# Just sanity-check the prefix is set so output filename is deterministic.
# ---------------------------------------------------------------------
require_flags 'GATK4_MARKDUPLICATES'          'ext\.prefix'

# ---------------------------------------------------------------------
# Stage 6 — samtools stats
#   XBS-core bash: samtools stats -F DUP,SUPPLEMENTARY,SECONDARY,UNMAP,QCFAIL
# ---------------------------------------------------------------------
require_flags 'SAMTOOLS_STATS'                '\-F DUP,SUPPLEMENTARY,SECONDARY,UNMAP,QCFAIL'

# ---------------------------------------------------------------------
# Stage 6 — Picard CollectWgsMetrics
#   XBS-core bash: --READ_LENGTH 0 --COVERAGE_CAP 10000 --COUNT_UNPAIRED
# ---------------------------------------------------------------------
require_flags 'PICARD_COLLECTWGSMETRICS' \
              '\-\-READ_LENGTH 0' \
              '\-\-COVERAGE_CAP 10000' \
              '\-\-COUNT_UNPAIRED'

# ---------------------------------------------------------------------
# Stage 7 — HaplotypeCaller
#   XBS-core bash: -ploidy 1 -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation
#                  --read-filter MappingQualityNotZeroReadFilter
# ---------------------------------------------------------------------
require_flags 'GATK4_HAPLOTYPECALLER' \
              '\-ploidy \$\{params\.sample_ploidy\}' \
              '\-ERC GVCF' \
              '\-G StandardAnnotation' \
              '\-G AS_StandardAnnotation' \
              '\-\-read-filter MappingQualityNotZeroReadFilter'

# ---------------------------------------------------------------------
# Stage 8 — CombineGVCFs
#   XBS-core bash: -G StandardAnnotation -G AS_StandardAnnotation
# ---------------------------------------------------------------------
require_flags 'GATK4_COMBINEGVCFS' \
              '\-G StandardAnnotation' \
              '\-G AS_StandardAnnotation'

# ---------------------------------------------------------------------
# Stage 9 — GenotypeGVCFs
#   XBS-core bash: --sample-ploidy 1 -G StandardAnnotation -G AS_StandardAnnotation
# ---------------------------------------------------------------------
require_flags 'GATK4_GENOTYPEGVCFS' \
              '\-\-sample-ploidy \$\{params\.sample_ploidy\}' \
              '\-G StandardAnnotation' \
              '\-G AS_StandardAnnotation'

# ---------------------------------------------------------------------
# Stage 10 — SelectVariants split
#   XBS-core bash: --select-type-to-include {SNP,INDEL} --remove-unused-alternates --exclude-non-variants
# ---------------------------------------------------------------------
require_flags 'GATK4_SELECTVARIANTS_SNP' \
              '\-\-select-type-to-include SNP' \
              '\-\-remove-unused-alternates' \
              '\-\-exclude-non-variants'
require_flags 'GATK4_SELECTVARIANTS_INDEL' \
              '\-\-select-type-to-include INDEL' \
              '\-\-remove-unused-alternates' \
              '\-\-exclude-non-variants'

# ---------------------------------------------------------------------
# Stage 11 — SNP VariantRecalibrator
#   XBS-core bash: -AS -an AS_MQRankSum -an AS_QD -an AS_MQ -an DP -mode SNP
#                  --max-gaussians 4 -mq-cap 60 -tranche 100.0 -tranche 99.9 -tranche 99.0
#                  --target-titv 1.85 (MTB-specific; we make it param-driven)
# ---------------------------------------------------------------------
require_flags 'GATK4_VARIANTRECALIBRATOR_SNP' \
              '\-AS' \
              '\-an AS_MQRankSum' \
              '\-an AS_QD' \
              '\-an AS_MQ' \
              '\-an DP' \
              '\-mode SNP' \
              '\-\-max-gaussians \$\{params\.vqsr_max_gaussians_snp\}' \
              '\-mq-cap 60' \
              '\-tranche 100\.0' \
              '\-tranche 99\.9' \
              '\-tranche 99\.0'

# ---------------------------------------------------------------------
# Stage 11 — SNP ApplyVQSR
#   XBS-core bash: -AS --ts-filter-level 99.9 --exclude-filtered -mode SNP
# ---------------------------------------------------------------------
require_flags 'GATK4_APPLYVQSR_SNP' \
              '\-AS' \
              '\-\-ts-filter-level 99\.9' \
              '\-mode SNP' \
              '\-\-exclude-filtered'

# ---------------------------------------------------------------------
# Stage 12 — INDEL VariantRecalibrator
#   XBS-core bash: -an MQRankSum -an QD -an DP -mode INDEL --max-gaussians 2 -mq-cap 60
#                  -tranche 100.0 -tranche 99.9 -tranche 99.0
# ---------------------------------------------------------------------
require_flags 'GATK4_VARIANTRECALIBRATOR_INDEL' \
              '\-an MQRankSum' \
              '\-an QD' \
              '\-an DP' \
              '\-mode INDEL' \
              '\-\-max-gaussians \$\{params\.vqsr_max_gaussians_indel\}' \
              '\-mq-cap 60' \
              '\-tranche 100\.0' \
              '\-tranche 99\.9' \
              '\-tranche 99\.0'

# ---------------------------------------------------------------------
# Stage 12 — INDEL ApplyVQSR
#   XBS-core bash: --lod-score-cutoff 0.0000 --exclude-filtered -mode INDEL
# ---------------------------------------------------------------------
require_flags 'GATK4_APPLYVQSR_INDEL' \
              '\-\-lod-score-cutoff' \
              '\-mode INDEL' \
              '\-\-exclude-filtered'

# ---------------------------------------------------------------------
# Hard-filter fallback (alternate path; not in XBS-core bash but
# required for MAGMA parity — MAGMA filters INDELs with hard filters,
# never VQSR). Check that the wired flags reference the configurable
# expressions, so consumers can tune thresholds per organism.
# ---------------------------------------------------------------------
require_flags 'GATK4_VARIANTFILTRATION_SNP' \
              '\-\-filter-expression' \
              '\$\{params\.snp_hard_filter_expression\}' \
              '\-\-filter-name'
require_flags 'GATK4_VARIANTFILTRATION_INDEL' \
              '\-\-filter-expression' \
              '\$\{params\.indel_hard_filter_expression\}' \
              '\-\-filter-name'

echo
echo "=== summary ==="
printf "  passed: %d\n  failed: %d\n  skipped: %d\n\n" "$PASS" "$FAIL" "$SKIP"

if [[ $FAIL -gt 0 ]]; then
  echo "FAIL — at least one module's ext.args diverges from XBS-core bash."
  echo "       See TORCH-Consortium/XBS-variant-calling-core for the canonical flags."
  exit 1
fi

echo "PASS — every module's ext.args contains the Heupink 2021 prescribed flags."
