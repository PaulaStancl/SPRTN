#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 06 - bcftools norm on sarek's somatic SNV/indel VCFs (Mutect2, Strelka2 snvs
#      + indels), so the callers' indels can be compared. Needs 02 finished.
#
#   ./06_normalize_vcfs.sh          # a few minutes; fine on the login node
#
# Same arguments as sarek's own --normalize_vcfs (off by default, not used in
# our run): left-align and trim indels against the GATK GRCh38 FASTA sarek
# aligned to, split multiallelic sites (-m -both), drop exact duplicates. MNVs
# are NOT split (no --atomize) - analysis/01_sarek.R does that itself.
# Output: $RESULTS_BASE/sarek/RJALS/normalized_bcftools/<pair>/*.norm.vcf.gz(.tbi)
# - analysis/01_sarek.R uses these instead of the raw VCFs once they exist.
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

pick() { find "$VC" -path "*/$PAIR/*" -name "$1" 2>/dev/null | head -1; }
for name in "${PAIR}.mutect2.filtered.vcf.gz" \
            "${PAIR}.strelka.somatic_snvs.vcf.gz" \
            "${PAIR}.strelka.somatic_indels.vcf.gz"; do
    in=$(pick "$name")
    [[ -n "$in" ]] || { warn "not found under $VC: $name - skipped"; continue; }
    out="$OUT/${name%.vcf.gz}.norm.vcf.gz"
    log "$name"
    # left-aligning can move a record past the next one, hence the sort
    bcftools norm -f "$FASTA" --multiallelics -both --rm-dup exact -Ou "$in" \
        | bcftools sort -T "$SORT_TMP" -Oz -o "$out"
    bcftools index -t -f "$out"
    ok "$(bcftools index -n "$in" 2>/dev/null || echo '?') records in -> $(bcftools index -n "$out") out: $out"
done
log "done: $OUT   (bcftools norm's own summary - lines total/split/realigned/skipped - is above)"
