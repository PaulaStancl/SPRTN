#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 06 - Subsample to ~30x tumour / ~20x normal so a full run fits in days, not a
#      week, on 8 cores. Same seed + same fraction on R1 and R2 keeps pairs.
#
#   ./06_subsample_fastq.sh
#   TUMOUR_FRACTION=1 ./06_subsample_fastq.sh     # skip subsampling for one side
#
# Do not go below ~20x tumour: AMBER/PURPLE (oncoanalyser) and ASCAT (sarek)
# lose het SNP sites and the purity/ploidy fit gets unstable.
# To run at FULL depth instead, skip this script and use:
#   ./07_make_samplesheets.sh --full
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_tools
mkdir -p "$SUB_DIR" "$LOG_DIR"

SEED=100
for s in "$TUMOUR_ID" "$NORMAL_ID"; do
    for r in R1 R2; do
        [[ -f "$RAW_DIR/${s}_${r}.fastq.gz.md5ok" ]] || die "$RAW_DIR/${s}_${r}.fastq.gz not verified - run ./05_download_test_data.sh"
    done
done

pids=()
for s in "$TUMOUR_ID" "$NORMAL_ID"; do
    frac=$NORMAL_FRACTION; [[ "$s" == "$TUMOUR_ID" ]] && frac=$TUMOUR_FRACTION
    for r in R1 R2; do
        in="$RAW_DIR/${s}_${r}.fastq.gz"
        out="$SUB_DIR/${s}_${r}.fastq.gz"
        [[ -f "$out.done" ]] && { ok "$(basename "$out") already done"; continue; }
        log "  $s $r  fraction $frac  -> $out"
        ( seqtk sample -s"$SEED" "$in" "$frac" | pigz -p 1 > "$out.tmp" \
          && mv "$out.tmp" "$out" && touch "$out.done" ) \
          > "$LOG_DIR/subsample_${s}_${r}.log" 2>&1 &
        pids+=($!)
    done
done

status=0
for p in "${pids[@]}"; do wait "$p" || status=1; done
[[ "$status" -eq 0 ]] || die "subsampling failed - see $LOG_DIR/subsample_*.log"

# Pairing check on the first 10k reads (cheap): read names must match R1 vs R2
for s in "$TUMOUR_ID" "$NORMAL_ID"; do
    if cmp -s <(pigz -dc "$SUB_DIR/${s}_R1.fastq.gz" | awk 'NR%4==1{sub(/[ \/].*/,""); print}' | head -10000) \
              <(pigz -dc "$SUB_DIR/${s}_R2.fastq.gz" | awk 'NR%4==1{sub(/[ \/].*/,""); print}' | head -10000); then
        ok "$s R1/R2 read names in sync"
    else
        die "$s R1/R2 out of sync - delete $SUB_DIR/${s}_* and re-run"
    fi
done

ls -lh "$SUB_DIR"/*.fastq.gz
log "Raw FASTQ no longer needed for the subsampled test:  rm $RAW_DIR/*.fastq.gz  (191 GB)"
log "Next:  ./07_make_samplesheets.sh"
