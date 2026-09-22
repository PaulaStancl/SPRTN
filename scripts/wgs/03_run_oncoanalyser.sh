#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 03 - nf-core/oncoanalyser (Hartwig WiGiTS): FASTQ -> bwa-mem2 -> REDUX ->
#      SAGE / ESVEE / AMBER / COBALT / PURPLE / LINX / LILAC / CHORD / CUPPA /
#      virus / ORANGE. RJALS tumour/normal DNA, GRCh38_hmf.
#      Run inside screen/tmux - several days at this depth. Independent of sarek.
#
#   ./03_run_oncoanalyser.sh
#
# Does not need SEX, so it can run before the patient's sex is confirmed.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_env
make_dirs
check_data

SHEET="$SAMPLESHEET_DIR/oncoanalyser_${DATASET}.csv"
OUT="$RESULTS_BASE/oncoanalyser/$DATASET"

# oncoanalyser only uses its prebuilt SIFs under the singularity profile; under
# apptainer the OCI conversion of hmftools-esvee fails (see 00_config.sh).
[[ "$NXF_PROFILE" == singularity ]] || die "oncoanalyser needs the singularity profile, got '$NXF_PROFILE'"
[[ -f "$SHEET" ]] || die "missing $SHEET - run ./01_make_samplesheets.sh"
[[ -f "$ONCO_REFDATA_CONFIG" ]] \
    || die "missing $ONCO_REFDATA_CONFIG - copy ../wgs_test/conf/oncoanalyser_refdata.config here"

log "outdir   : $OUT"

nf_run "oncoanalyser_${DATASET}" "$NXF_PROFILE" nf-core/oncoanalyser -r "$ONCOANALYSER_REV" \
    -c "$ONCO_REFDATA_CONFIG" \
    --input "$SHEET" \
    --outdir "$OUT" \
    --mode wgts \
    --genome GRCh38_hmf

log "oncoanalyser done: $OUT   (summary report: $OUT/$PATIENT/orange/)"
log "then: rm -rf $NXF_WORK_BASE/oncoanalyser_${DATASET}"
