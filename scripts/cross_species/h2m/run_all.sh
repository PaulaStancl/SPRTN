#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_all.sh - the human -> mouse prediction scripts in order (00 - 07b), on the server.
#
#   micromamba activate /common/WORK/pstancl/envs/cross-species-h2m
#   ./run_all.sh                         # config/config.yaml, Ensembl 110
#   ./run_all.sh --release latest        # same chain on the current Ensembl -> ens<N>/
#   CONFIG=/path/other.yaml ./run_all.sh
#   FROM=03 ./run_all.sh                 # start at a later step (earlier outputs must exist)
#   FREEZE=0 ./run_all.sh                # stop after the report, do not freeze
#
# Any extra arguments are passed to every script (e.g. --release latest).
# ---------------------------------------------------------------------------
set -euo pipefail
umask 077                                 # patient-derived outputs: readable only by you
cd "$(dirname "${BASH_SOURCE[0]}")"
CONFIG="${CONFIG:-config/config.yaml}"
FROM="${FROM:-00}"
FREEZE="${FREEZE:-1}"

for tool in Rscript python bcftools liftOver; do
    command -v "$tool" >/dev/null || { echo "missing $tool - activate the cross-species-h2m env" >&2; exit 1; }
done

steps=(00_resources 01_snv_coding_h2m 02_snv_noncoding_liftover 03_cn_to_mouse 04_sv_to_mouse
       05_signatures_expected 05b_pathways 06_build_matrix 06b_build_event_table 07_report)
[[ "$FREEZE" == 1 ]] && steps+=(07b_freeze)

for s in "${steps[@]}"; do
    [[ "${s:0:2}" < "${FROM:0:2}" ]] && continue
    echo "=== $s $(date '+%F %T')"
    Rscript "scripts/$s.R" --config "$CONFIG" "$@"
done
echo "=== done $(date '+%F %T')"
