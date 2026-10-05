#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 04 - nf-core/tumourevo: clonal/subclonal deconvolution from sarek output
#      (Mutect2 VCF + ASCAT copy number). Needs 02 to be finished.
#
#   ./04_run_tumourevo.sh
#   TEVO_TOOLS=tinc,mobster,viber,pyclone-vi,sparsesignatures ./04_run_tumourevo.sh
#   TEVO_FILTER=true ./04_run_tumourevo.sh    # clonality on CNAqc-PASS segments only
#
# tumourevo does NOT accept oncoanalyser/PURPLE output (CNA callers: ASCAT,
# sequenza, Battenberg, facets), so sarek is its only input here.
#
# Signature tools are OFF by default. They need a cohort, not one sample, and
# upstream does not test them either: the pipeline's own nf-test overrides tools
# to "tinc,mobster,pyclone-vi". On sparse input SparseSignatures returns NA for
# every cross-validation MSE and then dies picking K (`if (K < 2)`). Add
# sparsesignatures / sigprofiler through TEVO_TOOLS once the rest has run.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_env
make_dirs

TEVO_TOOLS="${TEVO_TOOLS:-tinc,mobster,viber,pyclone-vi}"
# --filter: true = subclonal / signature deconvolution only on segments that pass CNAqc
# (purity + copy number fit the mutations' VAF peaks); false (tumourevo's default) = all.
# A filtered run gets its own outdir, so the two can be compared; with -resume only the
# steps after CNAqc run again (while work/wgs/tumourevo_RJALS still exists).
TEVO_FILTER="${TEVO_FILTER:-false}"
[[ "$TEVO_FILTER" == true || "$TEVO_FILTER" == false ]] || die "TEVO_FILTER must be true or false"
SAREK_OUT="$RESULTS_BASE/sarek/$DATASET"
OUT="$RESULTS_BASE/tumourevo/$DATASET"
[[ "$TEVO_FILTER" == true ]] && OUT="$RESULTS_BASE/tumourevo/${DATASET}_cnaqcPASS"
PAIR="${TUMOUR_ID}_vs_${NORMAL_ID}"
FASTA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

[[ -n "$CANCER_TYPE" ]] || die "set CANCER_TYPE (IntOGen code of the tumour) in 00_config.sh"
[[ -d "$VEP_CACHE/homo_sapiens/${VEP_CACHE_VERSION}_GRCh38" ]] || die "VEP cache missing - run ../wgs_test/04_download_references.sh vep"
[[ -f "$FASTA" ]] || die "missing $FASTA - run ../wgs_test/04_download_references.sh sarek"

# ---- Locate sarek outputs ---------------------------------------------------
pick() { find "$SAREK_OUT/variant_calling" -name "$1" 2>/dev/null | head -1; }
VCF=$(pick "${PAIR}.mutect2.filtered.vcf.gz")
SEG=$(pick "${PAIR}.segments.txt")
PP=$(pick "${PAIR}.purityploidy.txt")
for v in VCF SEG PP; do
    [[ -n "${!v}" ]] || die "sarek output not found ($v) under $SAREK_OUT/variant_calling - did 02_run_sarek.sh finish with mutect2,ascat?"
    ok "$v = ${!v}"
done
[[ -f "$VCF.tbi" ]] || die "missing $VCF.tbi"

# VCF sample names must match tumour_sample / normal_sample exactly
T_SM="${PATIENT}_${TUMOUR_ID}"; N_SM="${PATIENT}_${NORMAL_ID}"
activate_tools
samples=$(bcftools query -l "$VCF" | tr '\n' ' ')
[[ " $samples " == *" $T_SM "* && " $samples " == *" $N_SM "* ]] \
    || die "VCF samples are '$samples', expected $T_SM and $N_SM"
ok "VCF samples: $samples"

# ---- Samplesheet ----------------------------------------------------------------
SHEET="$SAMPLESHEET_DIR/tumourevo_${DATASET}.csv"
cat > "$SHEET" <<EOF
dataset,patient,tumour_sample,normal_sample,vcf,tbi,cna_segments,cna_extra,cna_caller,cancer_type
$DATASET,$PATIENT,$T_SM,$N_SM,$VCF,$VCF.tbi,$SEG,$PP,ASCAT,$CANCER_TYPE
EOF
ok "$SHEET"

# tumourevo's per-task limits are tight: 2 h by default, 6 / 8 / 10 h for its
# process_low / medium / high labels. Nextflow kills a local task that passes its
# time (SIGTERM -> exit 143), which is what killed sarek's markdup twice. Inside a
# PBS job the walltime is the only limit that should apply. Times are not part of
# a task hash, so this does not affect resuming.
TUNING="$NXF_WORK_BASE/tumourevo_${DATASET}/tumourevo_tuning.config"
mkdir -p "$(dirname "$TUNING")"
{
    echo "// Written by 04_run_tumourevo.sh at $(date '+%F %T')."
    echo "process {"
    echo "    time = 240.h"
    for _lab in process_single process_low process_medium process_high process_long process_high_memory; do
        echo "    withLabel: $_lab { time = 240.h }"
    done
    echo "}"
} > "$TUNING"
log "task time: 240 h (tumourevo's default is 2 h)"
log "tools    : $TEVO_TOOLS"
log "filter   : $TEVO_FILTER (true = CNAqc-PASS segments only)"
log "outdir   : $OUT"

# tumourevo dev still uses pre-strict syntax; Nextflow 26.04 parses strict by default.
export NXF_SYNTAX_PARSER="${NXF_SYNTAX_PARSER:-v1}"

nf_run "tumourevo_${DATASET}" "$NXF_PROFILE" nf-core/tumourevo -r "$TUMOUREVO_REV" \
    --input "$SHEET" \
    --outdir "$OUT" \
    --genome GRCh38 \
    --fasta "$FASTA" \
    --tools "$TEVO_TOOLS" \
    --filter "$TEVO_FILTER" \
    --download_cache_vep false \
    --vep_cache "$VEP_CACHE" \
    --vep_cache_version "$VEP_CACHE_VERSION" \
    --vep_genome GRCh38 \
    --vep_species homo_sapiens \
    -c "$TUNING"

log "tumourevo done: $OUT"
