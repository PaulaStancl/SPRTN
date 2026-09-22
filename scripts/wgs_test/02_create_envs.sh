#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02 - micromamba envs.
#   nextflow-<version> : reused if it already exists (made for CHLOCK spatial)
#   wgs-tools          : seqtk, pigz, awscli, samtools, bcftools
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
command -v micromamba >/dev/null 2>&1 || die "micromamba not found"
make_dirs

# ---- Nextflow ---------------------------------------------------------------
if [[ -x "$ENV_PREFIX/bin/nextflow" ]]; then
    log "Nextflow env already exists - reusing: $ENV_PREFIX"
else
    log "Creating $ENV_PREFIX"
    micromamba create -y -p "$ENV_PREFIX" -f env_nextflow.yml
fi

# ---- WGS helper tools -------------------------------------------------------
if [[ -x "$TOOLS_PREFIX/bin/seqtk" ]]; then
    log "Tools env already exists - reusing: $TOOLS_PREFIX"
else
    log "Creating $TOOLS_PREFIX"
    micromamba create -y -p "$TOOLS_PREFIX" -f env_wgs_tools.yml
fi

# ---- Verify -----------------------------------------------------------------
activate_env
log "JAVA_HOME = $JAVA_HOME"
java -version 2>&1 | head -1
nf_ver=$(nextflow -v 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
[[ "$(printf '%s\n25.10.4\n' "$nf_ver" | sort -V | head -1)" == "25.10.4" ]] \
    && ok "Nextflow $nf_ver (sarek/oncoanalyser need >=25.10.4)" \
    || die "Nextflow $nf_ver is too old"

activate_tools
for t in seqtk pigz aws samtools bcftools wget; do
    command -v "$t" >/dev/null && ok "$t" || die "$t missing in $TOOLS_PREFIX"
done

log "Next:  ./03_pull_pipelines.sh"
