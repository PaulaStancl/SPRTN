#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 11 - nf-core/oncoanalyser (Hartwig WiGiTS): FASTQ -> bwa-mem2 -> REDUX ->
#      SAGE / ESVEE / AMBER / COBALT / PURPLE / LINX / LILAC / CHORD / CUPPA /
#      virus / ORANGE. Tumour/normal DNA only, GRCh38_hmf.
#      Run inside tmux/screen (1-2 days). Independent of sarek.
#
#   ./11_run_oncoanalyser.sh
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_env
make_dirs

SHEET="$SAMPLESHEET_DIR/oncoanalyser_${DATASET}.csv"
OUT="$RESULTS_BASE/oncoanalyser/$DATASET"

[[ -f "$SHEET" ]] || die "missing $SHEET - run ./07_make_samplesheets.sh"
[[ -f "$ONCO_REFDATA_CONFIG" ]] || die "missing $ONCO_REFDATA_CONFIG - run ./04_download_references.sh oncoanalyser"

log "outdir   : $OUT"

nf_run "oncoanalyser_${DATASET}" "$NXF_PROFILE" nf-core/oncoanalyser -r "$ONCOANALYSER_REV" \
    -c "$ONCO_REFDATA_CONFIG" \
    --input "$SHEET" \
    --outdir "$OUT" \
    --mode wgts \
    --genome GRCh38_hmf

log "oncoanalyser done: $OUT   (summary report: $OUT/$PATIENT/orange/)"
log "then: rm -rf $NXF_WORK_BASE/oncoanalyser_${DATASET}"
