#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 08 - every variant record in a region, from every VCF under results/wgs:
#      somatic (raw / unfiltered, filtered, normalised), germline, SV - sarek,
#      sarek extra callers (RJALS_vc), oncoanalyser (SAGE, PURPLE), tumourevo.
#
#   ./08_region_variants.sh                                # SPRTN (default)
#   ./08_region_variants.sh chr1:231338243-231338654 SPRTN_exon1
#
# Per file: records in the region (all / PASS), then the records (CHROM POS REF ALT
# FILTER, plus PAVE's IMPACT where PURPLE annotated it). gVCF reference blocks
# (ALT ".") are left out. A file without an index is read whole (slower).
# Output: printed and saved to $RESULTS_BASE/region_variants/<name>.txt
# SPRTN: Ensembl ENSG00000010072, GRCh38 chr1:231,337,104-231,375,416 (igv_regions.tsv)
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ./00_config.sh
activate_tools
umask 077

REGION="${1:-chr1:231337104-231375416}"
NAME="${2:-SPRTN}"
OUT_DIR="$RESULTS_BASE/region_variants"; mkdir -p "$OUT_DIR"
REPORT="$OUT_DIR/${NAME}.txt"
command -v bcftools >/dev/null || die "bcftools not in $TOOLS_PREFIX/bin"

{
echo "# $NAME  $REGION   $(date '+%F %T')"
echo "# per file: records in region (all / PASS), then CHROM POS REF ALT FILTER [IMPACT]"
find "$RESULTS_BASE" -path "$OUT_DIR" -prune -o -name '*.vcf.gz' -print | grep -v '/work/' | sort | while read -r f; do
    if [[ -f "$f.tbi" || -f "$f.csi" ]]; then sel=(-r "$REGION"); else sel=(-t "$REGION"); fi
    recs=$(bcftools view -H "${sel[@]}" -e 'ALT="."' "$f" 2>/dev/null | wc -l) || recs="?"
    pass=$(bcftools view -H "${sel[@]}" -f PASS -e 'ALT="."' "$f" 2>/dev/null | wc -l) || pass="?"
    echo "== ${f#$RESULTS_BASE/}   all: $recs   PASS: $pass"
    (( recs > 0 )) 2>/dev/null || continue
    if bcftools view -h "$f" 2>/dev/null | grep -q '^##INFO=<ID=IMPACT,'; then
        bcftools view "${sel[@]}" -e 'ALT="."' "$f" | bcftools query -f '  %CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/IMPACT\n'
    else
        bcftools view "${sel[@]}" -e 'ALT="."' "$f" | bcftools query -f '  %CHROM\t%POS\t%REF\t%ALT\t%FILTER\n'
    fi
done
} | tee "$REPORT"
log "saved: $REPORT"
