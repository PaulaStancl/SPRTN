#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_synteny.sh - compare our gene-order synteny blocks with public human-mouse synteny datasets.
# Needs the cross-species-h2m env and h2m step 03 (synteny_blocks.GRCm39.tsv) to have run.
#
#   micromamba activate /common/WORK/pstancl/envs/cross-species-h2m
#   ./run_synteny.sh                     # config.yaml next to this script
#   CONFIG=/path/other.yaml ./run_synteny.sh
# ---------------------------------------------------------------------------
set -euo pipefail
umask 077
cd "$(dirname "${BASH_SOURCE[0]}")"
for tool in Rscript liftOver; do
    command -v "$tool" >/dev/null || { echo "missing $tool - activate the cross-species-h2m env" >&2; exit 1; }
done
Rscript compare_synteny.R --config "${CONFIG:-config.yaml}"
