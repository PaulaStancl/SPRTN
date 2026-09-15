#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 10 - nf-core/sarek: FASTQ -> bwa-mem2 -> markdup -> BQSR -> somatic calling
#      Tumour/normal WGS on GATK.GRCh38. Run inside tmux/screen (1-3 days).
#
#   ./10_run_sarek.sh
#   SAREK_TOOLS=strelka,manta,ascat ./10_run_sarek.sh    # faster, no Mutect2
#
# Tools: mutect2 + ascat are what tumourevo (12) consumes; strelka + manta are
# the second SNV/indel caller and the SV caller. Mutect2 is the slowest step.
# No VEP here - tumourevo annotates with its own VEP.
# Re-running resumes (-resume) from the last finished task.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_env
make_dirs

SAREK_TOOLS="${SAREK_TOOLS:-mutect2,strelka,manta,ascat}"
SHEET="$SAMPLESHEET_DIR/sarek_${DATASET}.csv"
OUT="$RESULTS_BASE/sarek/$DATASET"

[[ -f "$SHEET" ]] || die "missing $SHEET - run ./07_make_samplesheets.sh"
[[ -d "$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem2Index" ]] \
    || die "iGenomes not staged - run ./04_download_references.sh sarek"

log "tools    : $SAREK_TOOLS"
log "outdir   : $OUT"

nf_run "sarek_${DATASET}" "$NXF_PROFILE" nf-core/sarek -r "$SAREK_REV" \
    --input "$SHEET" \
    --outdir "$OUT" \
    --genome GATK.GRCh38 \
    --igenomes_base "$IGENOMES_BASE" \
    --aligner bwa-mem2 \
    --tools "$SAREK_TOOLS"

log "sarek done: $OUT"
log "Next:  ./12_run_tumourevo.sh   (then: rm -rf $NXF_WORK_BASE/sarek_${DATASET})"
