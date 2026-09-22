#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01 - Write the sarek and oncoanalyser samplesheets from the Novogene delivery.
#      (tumourevo's samplesheet is built from sarek output by 04_run_tumourevo.sh)
#
#   ./01_make_samplesheets.sh
#
# Novogene names each lane-pair  <sample>_<library>_<flowcell>_L<n>_{1,2}.fq.gz
# and each sample here spans 3 lanes on 2 flowcells. Every lane-pair gets its
# own row with lane = <flowcell>_L<n>, so read groups (and sarek's BQSR, which
# models error per read group) stay separate per flowcell and lane.
#
# The sarek sheet is only written once SEX is set in 00_config.sh; the
# oncoanalyser sheet does not need it.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
mkdir -p "$SAMPLESHEET_DIR"
check_data
shopt -s nullglob

SAREK_SHEET="$SAMPLESHEET_DIR/sarek_${DATASET}.csv"
ONCO_SHEET="$SAMPLESHEET_DIR/oncoanalyser_${DATASET}.csv"
sarek_rows=""; onco_rows=""

# <sample> <sarek status: 0 normal / 1 tumour> <oncoanalyser sample_type>
add_sample() {
    local sample="$1" status="$2" stype="$3" r1 r2 base
    local r1s=("$FASTQ_DIR/$sample/${sample}_"*_1.fq.gz)
    (( ${#r1s[@]} > 0 )) || die "no ${sample}_*_1.fq.gz in $FASTQ_DIR/$sample"
    for r1 in "${r1s[@]}"; do
        r2="${r1%_1.fq.gz}_2.fq.gz"
        [[ -s "$r2" ]] || die "missing mate of $(basename "$r1"): $r2"
        base=$(basename "$r1" _1.fq.gz)
        [[ "${base#"${sample}_"}" =~ ^([^_]+)_([^_]+)_(L[0-9]+)$ ]] \
            || die "unexpected FASTQ name (want <sample>_<library>_<flowcell>_L<n>): $base"
        local lib="${BASH_REMATCH[1]}" lane="${BASH_REMATCH[2]}_${BASH_REMATCH[3]}"
        sarek_rows+="$PATIENT,$SEX,$status,$sample,$lane,$r1,$r2"$'\n'
        onco_rows+="$PATIENT,$PATIENT,$sample,$stype,dna,fastq,library_id:$lib;lane:$lane,$r1;$r2"$'\n'
        ok "$sample  $lane  ($lib)"
    done
}

log "FASTQ    : $FASTQ_DIR"
add_sample "$NORMAL_ID" 0 normal
add_sample "$TUMOUR_ID" 1 tumor

# ---- oncoanalyser: one row per lane, R1;R2 in filepath ------------------------
{ echo "group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath"
  printf '%s' "$onco_rows"; } > "$ONCO_SHEET"
ok "$ONCO_SHEET"

# ---- sarek: status 0 = normal, 1 = tumour ---------------------------------------
case "$SEX" in
    XX|XY|NA)
        { echo "patient,sex,status,sample,lane,fastq_1,fastq_2"
          printf '%s' "$sarek_rows"; } > "$SAREK_SHEET"
        ok "$SAREK_SHEET"
        [[ "$SEX" == NA ]] && warn "SEX=NA - ASCAT uses the patient's sex; set XX/XY if it is known"
        ;;
    "") rm -f "$SAREK_SHEET"
        warn "SEX not set in 00_config.sh - sarek sheet NOT written (oncoanalyser can run meanwhile)" ;;
    *)  die "SEX must be XX, XY or NA, got '$SEX'" ;;
esac

echo
for f in "$ONCO_SHEET" "$SAREK_SHEET"; do
    [[ -f "$f" ]] && { column -s, -t < "$f" | cut -c1-150; echo; }
done
log "Next:  ./02_run_sarek.sh  and/or  ./03_run_oncoanalyser.sh   (one after the other)"
