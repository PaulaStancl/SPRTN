#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 07 - Write the sarek and oncoanalyser samplesheets for the test pair.
#      (tumourevo's samplesheet is built from sarek output by 12_run_tumourevo.sh)
#
#   ./07_make_samplesheets.sh          # subsampled FASTQ (default)
#   ./07_make_samplesheets.sh --full   # full-depth raw FASTQ
#
# Sample names are the same in both pipelines. sarek writes VCF sample names as
# <patient>_<sample>, i.e. HCC1395_HCC1395T - tumourevo needs exactly that.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
mkdir -p "$SAMPLESHEET_DIR"

FQ_DIR="$SUB_DIR"
[[ "${1:-}" == "--full" ]] && FQ_DIR="$RAW_DIR"

for s in "$TUMOUR_ID" "$NORMAL_ID"; do
    for r in R1 R2; do
        [[ -s "$FQ_DIR/${s}_${r}.fastq.gz" ]] || die "missing $FQ_DIR/${s}_${r}.fastq.gz"
    done
done

T1="$FQ_DIR/${TUMOUR_ID}_R1.fastq.gz"; T2="$FQ_DIR/${TUMOUR_ID}_R2.fastq.gz"
N1="$FQ_DIR/${NORMAL_ID}_R1.fastq.gz"; N2="$FQ_DIR/${NORMAL_ID}_R2.fastq.gz"

# ---- sarek: status 0 = normal, 1 = tumour -------------------------------------
SAREK_SHEET="$SAMPLESHEET_DIR/sarek_${DATASET}.csv"
cat > "$SAREK_SHEET" <<EOF
patient,sex,status,sample,lane,fastq_1,fastq_2
$PATIENT,$SEX,0,$NORMAL_ID,L001,$N1,$N2
$PATIENT,$SEX,1,$TUMOUR_ID,L001,$T1,$T2
EOF
ok "$SAREK_SHEET"

# ---- oncoanalyser: one row per lane, R1;R2 in filepath -----------------------
ONCO_SHEET="$SAMPLESHEET_DIR/oncoanalyser_${DATASET}.csv"
cat > "$ONCO_SHEET" <<EOF
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
$PATIENT,$PATIENT,$TUMOUR_ID,tumor,dna,fastq,library_id:${TUMOUR_ID}_lib1;lane:001,$T1;$T2
$PATIENT,$PATIENT,$NORMAL_ID,normal,dna,fastq,library_id:${NORMAL_ID}_lib1;lane:001,$N1;$N2
EOF
ok "$ONCO_SHEET"

echo; column -s, -t < "$SAREK_SHEET" | cut -c1-150
echo; column -s, -t < "$ONCO_SHEET"  | cut -c1-150
echo
log "Next:  ./08_test_pipelines.sh all     (small nf-core test profiles)"
