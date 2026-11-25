#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Author: Paula Breyer
# 15 - sample identity check: are RJALS_Tm and RJALS_N the same individual? (somalier, from the BAM/CRAMs)
#
#   ./15_somalier_identity.sh
#
# somalier (Pedersen et al., Genome Medicine 2020) genotypes each sample at ~17.8k common human SNPs
# (its own GRCh38 site list) and compares them pairwise. Same individual: relatedness ~1 and IBS0 ~0
# (IBS0 = opposite homozygotes); different people: relatedness ~0 and IBS0 well above 0. A tumour keeps the
# person's germline genotypes, so a matched pair should look identical; LOH in the tumour lowers the
# relatedness a little (hets become homs) but does not create IBS0.
#
# Inputs: sarek's recalibrated CRAMs (falls back to md.cram / bam) and the GATK GRCh38 FASTA sarek used.
# somalier binary and site list: /common/WORK/pstancl/tools/somalier/ (installed by the osteosarcoma
# 00_somalier.sh; the site list is somalier's own sites.hg38.vcf.gz, chr-prefixed like GATK GRCh38).
# relate gets the tumour/normal pair as one expected group (-g), so the report flags it if they disagree.
# Output: $RESULTS_BASE/analysis/somalier/ (somalier.pairs.tsv, .samples.tsv, .html, RUN_LOG.md)
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ./00_config.sh
umask 077

SOM_DIR=/common/WORK/pstancl/tools/somalier
SOMALIER=$SOM_DIR/somalier
SITES=$SOM_DIR/sites.hg38.vcf.gz
SITES_URL=https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
FA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
OUT_DIR="$RESULTS_BASE/analysis/somalier"
PRE="$RESULTS_BASE/sarek/$DATASET/preprocessing"

[[ -x $SOMALIER ]] || die "somalier not found: $SOMALIER (see scripts/mouse_wgs/mouse_wgs_somalier.R header to install)"
if [[ ! -s $SITES ]]; then
    echo "downloading $SITES_URL"
    curl -sSfL -o "$SITES" "$SITES_URL" || die "no internet here? download sites.hg38.vcf.gz on the Mac and scp it to $SITES"
fi
[[ -r $FA && -r $FA.fai ]] || die "reference FASTA or .fai missing: $FA"

find_aln() {   # <sample>: first of recal.cram, md.cram, recal.bam, md.bam under sarek's preprocessing/
    local s="$1" p hits
    for p in "$s.recal.cram" "$s.md.cram" "$s.recal.bam" "$s.md.bam"; do
        mapfile -t hits < <(find "$PRE" -type f -name "$p" 2>/dev/null)
        if (( ${#hits[@]} == 1 )); then echo "${hits[0]}"; return; fi
        (( ${#hits[@]} > 1 )) && die "several $p under $PRE: ${hits[*]}"
    done
    die "no alignment for $s under $PRE"
}
T_ALN=$(find_aln "$TUMOUR_ID"); N_ALN=$(find_aln "$NORMAL_ID")

rm -rf "$OUT_DIR/extracted"; mkdir -p "$OUT_DIR/extracted"
LOG="$OUT_DIR/15_somalier_identity.log"
exec > >(tee "$LOG") 2>&1
echo "=== somalier identity $(date '+%F %T') on $(hostname) ==="
echo "tumour: $T_ALN"; echo "normal: $N_ALN"; echo "sites:  $SITES"; echo "fasta:  $FA"

# extract: one .somalier file per sample (named after the read-group SM tag)
for aln in "$N_ALN" "$T_ALN"; do
    "$SOMALIER" extract -d "$OUT_DIR/extracted" --sites "$SITES" -f "$FA" "$aln"
done
ls -l "$OUT_DIR"/extracted/*.somalier

# relate: tumour + normal declared as one expected individual (group)
smp=$(cd "$OUT_DIR/extracted" && ls *.somalier | sed 's/\.somalier$//' | paste -sd, -)
echo "$smp" > "$OUT_DIR/groups.txt"
"$SOMALIER" relate -g "$OUT_DIR/groups.txt" -o "$OUT_DIR/somalier" "$OUT_DIR"/extracted/*.somalier

echo; echo "=== pairs (expected_relatedness 1 = declared same individual) ==="
column -t "$OUT_DIR/somalier.pairs.tsv"
echo; echo "=== samples ==="
cut -f2,8-17 "$OUT_DIR/somalier.samples.tsv" | column -t
echo
echo "Same individual: relatedness near 1 and ibs0 near 0 (ibs0/n well below 0.01)."
echo "Different people: relatedness near 0 and ibs0 well above 0."
echo "HTML: $OUT_DIR/somalier.html"

# provenance
{
    echo "## $(date '+%F %T') - 15_somalier_identity.sh"
    echo "- Run by:   $(whoami)@$(hostname); $("$SOMALIER" 2>&1 | grep -m1 'somalier version')"
    echo "- Script:   scripts/wgs/15_somalier_identity.sh  git: $(git -C "$PROJECT_DIR" rev-parse --short HEAD 2>/dev/null || echo unknown)"
    echo "- Inputs:   $T_ALN, $N_ALN; sites $SITES; FASTA $FA"
    echo "- Command:  somalier extract --sites ... -f ...; somalier relate -g groups.txt (defaults: -d 7, --min-ab 0.3)"
    echo "- Outputs:  $OUT_DIR (somalier.pairs.tsv, somalier.samples.tsv, somalier.html, $LOG)"
    echo
    [[ -f "$OUT_DIR/RUN_LOG.md" ]] && cat "$OUT_DIR/RUN_LOG.md"
} > "$OUT_DIR/RUN_LOG.md.new" && mv "$OUT_DIR/RUN_LOG.md.new" "$OUT_DIR/RUN_LOG.md"
