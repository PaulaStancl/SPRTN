#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 03 - Pull the three pipelines at their pinned revisions into NXF_HOME.
#      Containers are pulled on first use into the shared cache
#      ($CONTAINER_DIR); 08_test_pipelines.sh triggers most of them.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_env
make_dirs

log "nf-core/sarek        -r $SAREK_REV"
nextflow pull nf-core/sarek -r "$SAREK_REV"

log "nf-core/oncoanalyser -r $ONCOANALYSER_REV"
nextflow pull nf-core/oncoanalyser -r "$ONCOANALYSER_REV"

log "nf-core/tumourevo    -r $TUMOUREVO_REV  (dev, no release yet)"
nextflow pull nf-core/tumourevo -r "$TUMOUREVO_REV"

echo
nextflow list
log "Next:  ./04_download_references.sh all   (run in tmux/screen)"
