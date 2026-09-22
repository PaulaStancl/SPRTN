#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 12 - nf-core/tumourevo: clonal/subclonal deconvolution + signatures from
#      sarek output (Mutect2 VCF + ASCAT copy number). Needs 10 to be finished.
#
#   ./12_run_tumourevo.sh
#   TEVO_TOOLS=tinc,mobster,viber,pyclone-vi,sparsesignatures ./12_run_tumourevo.sh
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
SAREK_OUT="$RESULTS_BASE/sarek/$DATASET"
OUT="$RESULTS_BASE/tumourevo/$DATASET"
PAIR="${TUMOUR_ID}_vs_${NORMAL_ID}"
FASTA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

[[ -d "$VEP_CACHE/homo_sapiens/${VEP_CACHE_VERSION}_GRCh38" ]] || die "VEP cache missing - run ./04_download_references.sh vep"
[[ -f "$FASTA" ]] || die "missing $FASTA - run ./04_download_references.sh sarek"

# ---- Locate sarek outputs ---------------------------------------------------
pick() { find "$SAREK_OUT/variant_calling" -name "$1" 2>/dev/null | head -1; }
VCF=$(pick "${PAIR}.mutect2.filtered.vcf.gz")
SEG=$(pick "${PAIR}.segments.txt")
PP=$(pick "${PAIR}.purityploidy.txt")
for v in VCF SEG PP; do
    [[ -n "${!v}" ]] || die "sarek output not found ($v) under $SAREK_OUT/variant_calling - did 10_run_sarek.sh finish with mutect2,ascat?"
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

log "tools    : $TEVO_TOOLS"
log "outdir   : $OUT"

# tumourevo dev still uses pre-strict syntax; Nextflow 26.04 parses strict by default.
export NXF_SYNTAX_PARSER="${NXF_SYNTAX_PARSER:-v1}"

nf_run "tumourevo_${DATASET}" "$NXF_PROFILE" nf-core/tumourevo -r "$TUMOUREVO_REV" \
    --input "$SHEET" \
    --outdir "$OUT" \
    --genome GRCh38 \
    --fasta "$FASTA" \
    --tools "$TEVO_TOOLS" \
    --download_cache_vep false \
    --vep_cache "$VEP_CACHE" \
    --vep_cache_version "$VEP_CACHE_VERSION" \
    --vep_genome GRCh38 \
    --vep_species homo_sapiens

log "tumourevo done: $OUT"
