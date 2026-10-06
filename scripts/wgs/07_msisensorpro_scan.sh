#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 07 - microsatellite list for sarek's MSIsensor-pro, built with the msisensor-pro
#      in sarek's own container from the GATK GRCh38 FASTA sarek aligned to.
#
#   qsub -I -q q2 -l select=1:ncpus=2:mem=16gb -l walltime=04:00:00   # interactive job, then:
#   cd /common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/wgs && ./07_msisensorpro_scan.sh
# single-threaded; roughly an hour for the whole genome.
#
# Why: iGenomes' Homo_sapiens_assembly38.msisensor_scan.list was made with an older
# msisensor-pro. sarek 3.10's msisensor-pro loads 0 sites from it and aborts with "Same
# reference genome file should be used in both scan and msi/pro/baseline steps" (the
# RJALS_vc run of 2026-10-06). 02_run_sarek.sh uses the list written here automatically.
# Output: $REF_BASE/msisensorpro/Homo_sapiens_assembly38.msisensorpro_scan.list
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"     # the script's own folder, wherever it is started from
source ./00_config.sh

FASTA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
OUT_DIR="$REF_BASE/msisensorpro"
LIST="$OUT_DIR/Homo_sapiens_assembly38.msisensorpro_scan.list"
[[ -f "$FASTA" ]] || die "missing $FASTA"
command -v singularity >/dev/null || die "singularity not in PATH"

# the msisensor-pro image sarek pulled for the RJALS_vc run (same version as its MSI step)
IMG=$(find "$CONTAINER_DIR" -maxdepth 1 -iname '*msisensor*pro*' \( -name '*.img' -o -name '*.sif' \) | sort | tail -1)
[[ -n "$IMG" ]] || die "no msisensor-pro image in $CONTAINER_DIR - run sarek with msisensorpro once (qsub_sarek_vc.sh) so it is pulled"
log "image : $IMG"
log "$(singularity exec "$IMG" msisensor-pro 2>&1 | grep -i -m1 version || echo 'version: ?')"
log "fasta : $FASTA"

mkdir -p "$OUT_DIR"
singularity exec -B "$(dirname "$FASTA")" -B "$OUT_DIR" "$IMG" \
    msisensor-pro scan -d "$FASTA" -o "$LIST.tmp"
mv "$LIST.tmp" "$LIST"
ok "$(($(wc -l < "$LIST") - 1)) microsatellite sites -> $LIST"
log "now: cd scripts/wgs && qsub qsub_sarek_vc.sh   (02 picks the list up automatically)"
