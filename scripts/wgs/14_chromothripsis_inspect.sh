#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 14 - chromothripsis step 1: find the inputs for run_chromothripsis.R in the oncoanalyser and
#      sarek results, and print what the R script needs to know about them - metadata only
#      (paths, sizes, record counts, VCF meta lines, table column names). No data rows.
#
#   ./14_chromothripsis_inspect.sh [oncoanalyser_dir] [sarek_dir]
#
# Defaults: $RESULTS_BASE/oncoanalyser/RJALS/RJALS and $RESULTS_BASE/sarek/RJALS
# Output: stdout and $RESULTS_BASE/analysis/chromothripsis/inspect.txt
# ---------------------------------------------------------------------------
set -eu   # no pipefail: the "| head" pipelines below end early by design (SIGPIPE upstream)
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ./00_config.sh
activate_tools
umask 077

ONCO="${1:-$RESULTS_BASE/oncoanalyser/$PATIENT/$PATIENT}"
SAREK="${2:-$RESULTS_BASE/sarek/$PATIENT}"
OUT_DIR="$RESULTS_BASE/analysis/chromothripsis"; mkdir -p "$OUT_DIR"
exec > >(tee "$OUT_DIR/inspect.txt") 2>&1

section() { printf '\n==== %s ====\n' "$*"; }
size()    { du -h "$1" | cut -f1; }

# every file matching a name pattern under a folder, with size; MISSING if none
find_all() {   # find_all <label> <dir> <name glob> [extra find args]
    local label=$1 dir=$2 glob=$3; shift 3
    local hits
    hits=$(find "$dir" -name "$glob" "$@" 2>/dev/null | sort || true)
    if [[ -z "$hits" ]]; then printf '%-28s MISSING  (%s under %s)\n' "$label" "$glob" "$dir"; return 0; fi
    while read -r f; do printf '%-28s %6s  %s\n' "$label" "$(size "$f")" "$f"; done <<< "$hits"
}

# VCF: meta lines that matter (no records) - reference, source, contig naming, FILTER and the
# INFO/FORMAT keys used to pair breakends and read strands, sample columns, record counts
vcf_meta() {
    local f=$1
    section "VCF $f"
    bcftools view -h "$f" | grep -E '^##(fileformat|reference|source|assembly)' | cut -c1-200 || true
    bcftools view -h "$f" | grep -E '^##(gridss|esvee|cmdline|.*[Cc]ommand)' | cut -c1-200 | head -5 || true
    local contigs; contigs=$(bcftools view -h "$f" | grep -c '^##contig=' || true)
    echo "contig lines: $contigs; first 3 / chrX / chrY / chrM entries:"
    bcftools view -h "$f" | grep '^##contig=' | head -3 | cut -c1-120
    bcftools view -h "$f" | grep -E '^##contig=<ID=(chr)?(X|Y|M|MT),' | cut -c1-120 || true
    echo "FILTER IDs: $(bcftools view -h "$f" | grep '^##FILTER=' | sed 's/.*ID=\([^,]*\),.*/\1/' | tr '\n' ' ')"
    echo "INFO IDs:   $(bcftools view -h "$f" | grep '^##INFO=' | sed 's/.*ID=\([^,]*\),.*/\1/' | tr '\n' ' ')"
    echo "FORMAT IDs: $(bcftools view -h "$f" | grep '^##FORMAT=' | sed 's/.*ID=\([^,]*\),.*/\1/' | tr '\n' ' ')"
    echo "samples:    $(bcftools query -l "$f" | tr '\n' ' ')"
    echo "records: $(bcftools view -H "$f" | wc -l)   PASS: $(bcftools view -H -f PASS,. "$f" | wc -l)"
    echo "records per contig naming (counts only):"
    bcftools view -H "$f" | cut -f1 | sed -E 's/^chr//' | sort | uniq -c | sort -k2,2V | awk '{printf "  %s:%s", $2, $1} END {print ""}'
}

# table: column names and row count only
tsv_meta() {
    local f=$1
    section "TABLE $f"
    local cat=cat; [[ "$f" == *.gz ]] && cat="zcat"
    echo "columns: $($cat "$f" | grep -v '^##' | head -1 | tr '\t' ' ')"
    echo "rows (excl. header): $(( $($cat "$f" | grep -vc '^#') - 1 ))"
}

section "inputs"
echo "oncoanalyser dir: $ONCO"; [[ -d "$ONCO" ]] || echo "  MISSING"
echo "sarek dir:        $SAREK"; [[ -d "$SAREK" ]] || echo "  MISSING"
echo "sample ids (00_config.sh): tumour $TUMOUR_ID  normal $NORMAL_ID"

section "oncoanalyser files"
find_all purple.sv.vcf         "$ONCO" '*.purple.sv.vcf.gz'
find_all esvee.somatic.vcf     "$ONCO" '*.esvee.somatic.vcf.gz'
find_all gridss.vcf            "$ONCO" '*gridss*.vcf.gz'
find_all purple.cnv.somatic    "$ONCO" '*.purple.cnv.somatic.tsv'
find_all purple.purity         "$ONCO" '*.purple.purity.tsv'
find_all purple.segment        "$ONCO" '*.purple.segment.tsv'
find_all amber.baf             "$ONCO" '*.amber.baf.tsv.gz'
find_all linx.clusters         "$ONCO" '*.linx.clusters.tsv'
find_all linx.svs              "$ONCO" '*.linx.svs.tsv'
find_all linx.breakend         "$ONCO" '*.linx.breakend.tsv'
find_all linx.vis_sv_data      "$ONCO" '*.linx.vis_sv_data.tsv'
find_all linx.vis_segments     "$ONCO" '*.linx.vis_segments.tsv'
find_all linx.vis_copy_number  "$ONCO" '*.linx.vis_copy_number.tsv'
echo "LINX plots: $(find "$ONCO" -path '*linx*' -name '*.png' 2>/dev/null | wc -l) png"
find "$ONCO" -path '*linx*' -name '*.png' 2>/dev/null | sed 's|.*/||' | sort | head -40 | sed 's/^/  /'
find_all linxreport.html       "$ONCO" '*linx*.html'
find_all orange.pdf            "$ONCO" '*.orange.pdf'
find_all orange.json           "$ONCO" '*.orange.json'

section "sarek files"
find_all manta.somatic_sv       "$SAREK" '*.somatic_sv.vcf.gz' -path '*manta*'
find_all ascat.segments        "$SAREK" '*.segments.txt'    -path '*ascat*'
find_all ascat.purityploidy    "$SAREK" '*purityploidy*'    -path '*ascat*'
find_all ascat.metrics         "$SAREK" '*metrics*'         -path '*ascat*'
find_all ascat.cnvs            "$SAREK" '*.cnvs.txt'        -path '*ascat*'

section "versions (LINX jar / container for SvVisualiser)"
find "$ONCO/.." "$ONCO/../.." -maxdepth 3 -path '*pipeline_info*' \( -name '*versions*.yml' -o -name 'params_*.json' \) 2>/dev/null | sort -u | head -5
for f in $(find "$ONCO/.." "$ONCO/../.." -maxdepth 3 -path '*pipeline_info*' -name '*versions*.yml' 2>/dev/null | sort -u); do
    grep -i -A2 'linx' "$f" | head -12 || true
done
sif=$(find "${NXF_SINGULARITY_CACHEDIR:-/nonexistent}" -maxdepth 1 -iname '*linx*' 2>/dev/null | head -3 || true)
echo "LINX images in NXF_SINGULARITY_CACHEDIR (${NXF_SINGULARITY_CACHEDIR:-unset}): ${sif:-none found}"

section "VCF headers"
for f in $(find "$ONCO" -name '*.purple.sv.vcf.gz' 2>/dev/null) \
         $(find "$SAREK" -path '*manta*' -name '*.somatic_sv.vcf.gz' 2>/dev/null); do
    vcf_meta "$f"
done

section "table headers"
for f in $(find "$ONCO" \( -name '*.purple.cnv.somatic.tsv' -o -name '*.purple.purity.tsv' \
                          -o -name '*.linx.clusters.tsv' -o -name '*.linx.svs.tsv' \
                          -o -name '*.linx.vis_sv_data.tsv' -o -name '*.linx.breakend.tsv' \) 2>/dev/null) \
         $(find "$SAREK" -path '*ascat*' \( -name '*.segments.txt' -o -name '*purityploidy*' -o -name '*metrics*' \) 2>/dev/null); do
    tsv_meta "$f"
done

# purity / ploidy are sample-level summaries (already in the README), printed in full
section "purity / ploidy"
for f in $(find "$ONCO" -name '*.purple.purity.tsv' 2>/dev/null); do
    echo "$f"; cut -f1-6 "$f" | column -t
done
for f in $(find "$SAREK" -path '*ascat*' -name '*purityploidy*' 2>/dev/null); do
    echo "$f"; column -t "$f"
done

echo; echo "written: $OUT_DIR/inspect.txt"
