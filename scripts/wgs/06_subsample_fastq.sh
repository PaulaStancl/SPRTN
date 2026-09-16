#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 06 - Subsample to ~30x tumour / ~20x normal so a full run fits in days, not a
#      week, on 8 cores. Same seed + same fraction on R1 and R2 keeps pairs.
#
#   ./06_subsample_fastq.sh
#   TUMOUR_FRACTION=1 ./06_subsample_fastq.sh     # keep that side at full depth
#   SKIP_PAIR_CHECK=1 ./06_subsample_fastq.sh     # skip the full R1/R2 name check
#
# A fraction of 1 (or more) is NOT passed to seqtk. Its usage line documents the
# argument as `<frac>|<number>` but never says where the switch happens; seqtk
# has no manual beyond that. The threshold is only in the source, seqtk.c:
#   if (frac >= 1.0) num = (uint64_t)(frac + .499), frac = 0.;
# so `seqtk sample in.fq 1` returns ONE read, not every read. Full depth is
# handled here by symlinking the raw FASTQ instead.
#
# Coverage: ~30x tumour / ~20x normal is a compromise for a test run on 8 cores,
# not a validated threshold. Lower coverage mainly costs heterozygous SNP sites
# and somatic SNV recall, so purity/ploidy fits (AMBER/COBALT/PURPLE, ASCAT) get
# noisier. AMBER keeps sites with tumour depth >= 8 and reference depth within
# 50-150% of the genome-wide median, so as depth drops more sites fall out.
# ASCAT documents a >50x preset for its tumour-only WGS mode and no minimum for
# matched mode. Hartwig's own WGS cohort is ~100x tumour / ~30x normal.
#   seqtk sample       https://github.com/lh3/seqtk (README: same -s seed keeps pairs)
#   AMBER              https://github.com/hartwigmedical/hmftools/blob/master/amber/README.md
#   PURPLE             https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md
#   ASCAT              https://github.com/VanLoo-lab/ascat (WGS_hg38_50X preset)
#   SEQC2 depth study  https://doi.org/10.1038/s41587-021-00994-5
#
# To run everything at FULL depth, skip this script and use:
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
        # Trust .done only if the output is actually there and non-empty
        [[ -s "$out" && -f "$out.done" ]] && { ok "$(basename "$out") already done"; continue; }
        rm -f "$out" "$out.tmp" "$out.done"

        # seqtk reads any fraction >= 1 as a read COUNT - keep full depth by symlink
        if awk -v f="$frac" 'BEGIN{exit !(f+0 >= 1)}'; then
            ln -s "$in" "$out" && touch "$out.done"
            ok "$s $r  full depth (symlink to raw)"
            continue
        fi

        log "  $s $r  fraction $frac  -> $out"
        ( seqtk sample -s"$SEED" "$in" "$frac" | pigz -p 1 > "$out.tmp" \
          && mv "$out.tmp" "$out" && touch "$out.done" ) \
          > "$LOG_DIR/subsample_${s}_${r}.log" 2>&1 &
        pids+=($!)
    done
done

status=0
for p in ${pids[@]+"${pids[@]}"}; do wait "$p" || status=1; done
[[ "$status" -eq 0 ]] || die "subsampling failed - see $LOG_DIR/subsample_*.log"

# Pairing check over the WHOLE file, not just the head: sampling R1 and R2
# separately is only correct if the raw files hold the same records in the same
# order, and a mismatch anywhere shifts every name after it. Minutes per sample.
if [[ -n "${SKIP_PAIR_CHECK:-}" ]]; then
    warn "SKIP_PAIR_CHECK set - R1/R2 pairing not verified"
else
    names() { pigz -dc "$1" | awk 'NR%4==1{sub(/[ \/].*/,""); print}'; }
    for s in "$TUMOUR_ID" "$NORMAL_ID"; do
        log "  checking R1/R2 read names: $s"
        cmp -s <(names "$SUB_DIR/${s}_R1.fastq.gz") <(names "$SUB_DIR/${s}_R2.fastq.gz") \
            && ok "$s R1/R2 in sync (all reads)" \
            || die "$s R1/R2 out of sync - delete $SUB_DIR/${s}_* and re-run"
    done
fi

ls -lh "$SUB_DIR"/*.fastq.gz
log "Raw FASTQ no longer needed for the subsampled test:  rm $RAW_DIR/*.fastq.gz  (191 GB)"
log "Next:  ./07_make_samplesheets.sh"
