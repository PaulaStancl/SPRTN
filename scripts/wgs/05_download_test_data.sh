#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 05 - Download the test tumour/normal WGS pair + somatic truth set.
#
# SEQC2 HCC1395 (breast cancer cell line) / HCC1395BL (matched B-lymphocytes)
# BioProject PRJNA489865, Fudan site (WGS_FD), HiSeq X Ten, 2x150, one run each
#   tumour  SRR7890829  WGS_FD_T_1  550 M read pairs  ~53x  93.7 GB
#   normal  SRR7890826  WGS_FD_N_1  570 M read pairs  ~55x  97.4 GB
# Truth set: SEQC2 v1.2.1 high-confidence somatic SNV/indel, GRCh38 (chr names)
#
# Resumable (wget -c) and md5-verified. Run inside tmux/screen - several hours.
# ENA FASTQ headers carry no flowcell/lane, so the samplesheets use lane 001.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_tools
mkdir -p "$RAW_DIR" "$DATA_DIR/truth_set"

ENA=https://ftp.sra.ebi.ac.uk/vol1/fastq
#          local name              url                                              md5
FILES=(
  "${TUMOUR_ID}_R1.fastq.gz  $ENA/SRR789/009/SRR7890829/SRR7890829_1.fastq.gz  e182a4ae1f85db81e2d4c5ba98ad85d0"
  "${TUMOUR_ID}_R2.fastq.gz  $ENA/SRR789/009/SRR7890829/SRR7890829_2.fastq.gz  a45365530f316e9e87dd2a2388705b84"
  "${NORMAL_ID}_R1.fastq.gz  $ENA/SRR789/006/SRR7890826/SRR7890826_1.fastq.gz  2feac99549d4a0c63f571d7cedc7a7ee"
  "${NORMAL_ID}_R2.fastq.gz  $ENA/SRR789/006/SRR7890826/SRR7890826_2.fastq.gz  5d5b94c0242b899c09ae24820daaa937"
)

# ---- Truth set (small) ------------------------------------------------------
TRUTH=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest
log "Truth set -> $DATA_DIR/truth_set"
for f in high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz \
         high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz \
         High-Confidence_Regions_v1.2.bed; do
    wget -c -q -O "$DATA_DIR/truth_set/$f" "$TRUTH/$f" && ok "$f"
done

# ---- FASTQ: 4 files in parallel ------------------------------------------------
log "FASTQ (191 GB) -> $RAW_DIR"
cd "$RAW_DIR"
pids=()
for entry in "${FILES[@]}"; do
    read -r name url md5 <<<"$entry"
    if [[ -f "$name.md5ok" ]]; then ok "$name already verified"; continue; fi
    (
        wget -c -q -O "$name" "$url" || { echo "download failed: $name" >&2; exit 1; }
        echo "$md5  $name" | md5sum -c --quiet - || { echo "md5 MISMATCH: $name (delete it and re-run)" >&2; exit 1; }
        touch "$name.md5ok"
    ) > "$LOG_DIR/download_$name.log" 2>&1 &
    pids+=($!)
    log "  started $name  (log: $LOG_DIR/download_$name.log)"
done

status=0
for p in "${pids[@]}"; do wait "$p" || status=1; done
[[ "$status" -eq 0 ]] || die "a download failed - see $LOG_DIR/download_*.log, then re-run (resumes)"

ls -lh "$RAW_DIR"/*.fastq.gz
log "All FASTQ verified. Next:  ./06_subsample_fastq.sh"
