#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02 - nf-core/sarek: FASTQ -> bwa-mem2 -> markdup -> BQSR -> somatic calling
#      RJALS tumour/normal WGS on GATK.GRCh38. Run inside screen/tmux - on 8
#      cores at this depth expect about a week (see README).
#
#   ./02_run_sarek.sh
#   SAREK_TOOLS=strelka,manta,ascat ./02_run_sarek.sh    # faster, no Mutect2
#
# Tools: mutect2 + ascat are what tumourevo (04) consumes; strelka + manta are
# the second SNV/indel caller and the SV caller. Mutect2 is the slowest step.
# No VEP here - tumourevo annotates with its own VEP.
# Re-running resumes (-resume) from the last finished task.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_env
make_dirs
check_data

SAREK_TOOLS="${SAREK_TOOLS:-mutect2,strelka,manta,ascat}"
SHEET="$SAMPLESHEET_DIR/sarek_${DATASET}.csv"
OUT="$RESULTS_BASE/sarek/$DATASET"

[[ -f "$SHEET" ]] || die "missing $SHEET - set SEX in 00_config.sh, then run ./01_make_samplesheets.sh"
# A sheet written before SEX was changed would start a run with the wrong sex -
# and fixing that later restarts every task.
sheet_sex=$(awk -F, 'NR > 1 { print $2 }' "$SHEET" | sort -u | tr '\n' ' ')
[[ "$sheet_sex" == "$SEX " ]] || die "$SHEET has sex '$sheet_sex' but 00_config.sh says '$SEX' - re-run ./01_make_samplesheets.sh"
[[ -d "$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem2Index" ]] \
    || die "iGenomes not staged - run ../wgs_test/04_download_references.sh sarek"

# Per-launch tuning, written next to the run:
#
# 1. Time. sarek's default is 8 h per task (16 h / 32 h for its process_medium /
#    process_high labels), and GATK4_MARKDUPLICATES sets only cpus and memory, so
#    it inherits 8 h. On this 90x pair it needs more, and Nextflow kills a local
#    task that passes its time (SIGTERM -> exit 143); that killed the runs of
#    2026-09-22 and 2026-09-23 at exactly 8 h. Inside a PBS job the walltime is
#    the only limit that should apply, so lift every task to it.
# 2. bwa-mem2 asks 24 cpus per chunk, so in a job of 32+ cpus only one chunk fits
#    and the rest idle through the whole alignment. Give it half the job instead.
#
# Neither changes a task's command line, so a resume still reuses finished tasks.
TUNING="$NXF_WORK_BASE/sarek_${DATASET}/sarek_tuning.config"
mkdir -p "$(dirname "$TUNING")"
{
    echo "// Written by 02_run_sarek.sh at $(date '+%F %T')."
    echo "process {"
    echo "    time = 240.h"
    for _lab in process_single process_low process_medium process_high process_long; do
        echo "    withLabel: $_lab { time = 240.h }"
    done
    if (( ${NCPUS:-0} >= 32 )); then
        echo "    withName: 'BWAMEM2_MEM' { cpus = $(( NCPUS / 2 )) }"
    fi
    echo "}"
} > "$TUNING"
log "task time: 240 h (sarek's default is 8 h - it killed markdup twice)"
(( ${NCPUS:-0} >= 32 )) && log "bwa-mem2 : $(( NCPUS / 2 )) cpus per chunk, two chunks at a time"
log "sex      : $SEX"
log "tools    : $SAREK_TOOLS"
log "outdir   : $OUT"

nf_run "sarek_${DATASET}" "$NXF_PROFILE" nf-core/sarek -r "$SAREK_REV" \
    --input "$SHEET" \
    --outdir "$OUT" \
    --genome GATK.GRCh38 \
    --igenomes_base "$IGENOMES_BASE" \
    --aligner bwa-mem2 \
    --tools "$SAREK_TOOLS" \
    -c "$TUNING"

log "sarek done: $OUT"
log "Next:  ./04_run_tumourevo.sh   (then: rm -rf $NXF_WORK_BASE/sarek_${DATASET})"
