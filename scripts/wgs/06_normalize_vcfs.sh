#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 06 - bcftools norm on every somatic SNV/indel VCF the analysis compares, so the
#      callers' indels are written the same way. Needs 02 (and 03 for SAGE/PURPLE).
#
#   ./06_normalize_vcfs.sh          # a few minutes; fine on the login node
#
# Same arguments as sarek's own --normalize_vcfs (off by default, not used in our runs):
# left-align and trim indels against the GATK GRCh38 FASTA sarek aligned to, split
# multiallelic sites (-m -both), drop exact duplicates. MNVs are NOT split (no --atomize) -
# the analysis scripts do that themselves. Inputs, whichever exist:
#   sarek/RJALS and sarek/RJALS_vc (extra callers, SAREK_STEP=variant_calling):
#     Mutect2 filtered, Strelka2 somatic snvs + indels, MuSE, FreeBayes, LoFreq
#   oncoanalyser: PURPLE's final somatic VCF (SAGE + PAVE + PURPLE) - Hartwig's GRCh38 has the
#     same coordinates but some masked bases, so a REF mismatch there only warns (-c w)
# Output: $RESULTS_BASE/sarek/RJALS/normalized_bcftools/<pair>/*.norm.vcf.gz(.tbi) and
#         $RESULTS_BASE/oncoanalyser/RJALS/normalized_bcftools/*.norm.vcf.gz(.tbi)
# - analysis/01_sarek.R and 03_compare_callers.R use these instead of the raw VCFs.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_tools
umask 077

PAIR="${TUMOUR_ID}_vs_${NORMAL_ID}"
VC="$RESULTS_BASE/sarek/$DATASET/variant_calling"
OUT="$RESULTS_BASE/sarek/$DATASET/normalized_bcftools/$PAIR"
FASTA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

command -v bcftools >/dev/null || die "bcftools not in $TOOLS_PREFIX/bin"
[[ -f "$FASTA" && -f "$FASTA.fai" ]] || die "missing $FASTA (+ .fai) - run ../wgs_test/04_download_references.sh sarek"
mkdir -p "$OUT"
SORT_TMP=$(mktemp -d "${TMPDIR:-/tmp}/bcfsort.XXXXXX"); trap 'rm -rf "$SORT_TMP"' EXIT

# norm_one <input vcf> <output dir> [check-ref mode]
norm_one() {
    local in="$1" od="$2" cref="${3:-e}" out
    mkdir -p "$od"
    out="$od/$(basename "${in%.vcf.gz}").norm.vcf.gz"
    log "$(basename "$in")"
    # left-aligning can move a record past the next one, hence the sort
    bcftools norm -f "$FASTA" -c "$cref" --multiallelics -both --rm-dup exact -Ou "$in" \
        | bcftools sort -T "$SORT_TMP" -Oz -o "$out"
    bcftools index -t -f "$out"
    ok "$(bcftools index -n "$in" 2>/dev/null || echo '?') records in -> $(bcftools index -n "$out") out: $out"
}

# sarek callers, in the main run and in the extra-callers run (RJALS_vc)
n=0
for vc in "$VC" "$RESULTS_BASE/sarek/${DATASET}_vc/variant_calling"; do
    [[ -d "$vc" ]] || continue
    while IFS= read -r in; do
        norm_one "$in" "$OUT"; n=$((n + 1))
    done < <(find "$vc" -path "*/$PAIR/*" \( -name "${PAIR}.mutect2.filtered.vcf.gz" \
                  -o -name "${PAIR}.strelka.somatic_snvs.vcf.gz" -o -name "${PAIR}.strelka.somatic_indels.vcf.gz" \
                  -o -name "${PAIR}.muse.vcf.gz" -o -name "${PAIR}.freebayes.vcf.gz" -o -name "${PAIR}*lofreq*.vcf.gz" \) | sort)
done
(( n )) || warn "no sarek SNV/indel VCFs found under $VC"

# oncoanalyser: PURPLE's final somatic VCF
PURPLE_VCF="$RESULTS_BASE/oncoanalyser/$DATASET/$DATASET/purple/${TUMOUR_ID}.purple.somatic.vcf.gz"
if [[ -f "$PURPLE_VCF" ]]; then
    norm_one "$PURPLE_VCF" "$RESULTS_BASE/oncoanalyser/$DATASET/normalized_bcftools" w
else
    warn "not found: $PURPLE_VCF - skipped"
fi
log "done: $OUT   (bcftools norm's own summary - lines total/split/realigned/skipped - is above)"
