#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Shared configuration for the WGS nf-core pipelines on lobsang.
#
#   source 00_config.sh
#
# Every other script sources this one - change paths and versions HERE.
# Anything already exported in your shell wins, e.g.
#   SAREK_TOOLS=strelka,manta,ascat ./10_run_sarek.sh
# ---------------------------------------------------------------------------

_CFG_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export WGS_SCRIPTS="$_CFG_DIR"

# ---- Shared server layout (same as every other project) --------------------
export WORK_BASE="${WORK_BASE:-/common/WORK/pstancl}"
export ENV_ROOT="${ENV_ROOT:-$WORK_BASE/envs}"
export PROGRAMS_DIR="${PROGRAMS_DIR:-$WORK_BASE/PROGRAMI}"
export CONTAINER_DIR="${CONTAINER_DIR:-$WORK_BASE/singularity_cache}"   # ONE cache for all pipelines
export NXF_HOME="${NXF_HOME:-$PROGRAMS_DIR/nextflow}"                   # pulled pipelines live here
export TMPDIR="${TMPDIR:-$WORK_BASE/tmp}"                               # never $HOME, never /tmp

# ---- Shared references (outside the project, reused by other projects) -----
export REF_BASE="${REF_BASE:-$WORK_BASE/references}"
export IGENOMES_BASE="${IGENOMES_BASE:-$REF_BASE/igenomes}"             # sarek: GATK.GRCh38
export HMF_REF_DIR="${HMF_REF_DIR:-$REF_BASE/hmf/oncoanalyser}"         # oncoanalyser: GRCh38_hmf
export VEP_CACHE="${VEP_CACHE:-$REF_BASE/vep_cache}"                    # tumourevo: VEP cache

# ---- Pipeline versions ------------------------------------------------------
export SAREK_REV="${SAREK_REV:-3.10.0}"                  # released 2026-08-12
export ONCOANALYSER_REV="${ONCOANALYSER_REV:-3.0.0}"     # released 2026-09-01
# tumourevo has no release yet - pinned to the dev HEAD of 2026-09-08
export TUMOUREVO_REV="${TUMOUREVO_REV:-738cb052fd51f47563ff0f96eaa4b2eb9d4b44f3}"
export VEP_CACHE_VERSION="${VEP_CACHE_VERSION:-115}"     # tumourevo ships ensembl-vep 115.2

# ---- Nextflow env: version read from env_nextflow.yml (single source) -----
export NEXTFLOW_VERSION="${NEXTFLOW_VERSION:-$(sed -n 's/^[[:space:]]*-[[:space:]]*nextflow=\([0-9][0-9.]*\).*/\1/p' "$_CFG_DIR/env_nextflow.yml" | head -1)}"
export ENV_PREFIX="${ENV_PREFIX:-$ENV_ROOT/nextflow-$NEXTFLOW_VERSION}"
export TOOLS_PREFIX="${TOOLS_PREFIX:-$ENV_ROOT/wgs-tools}"              # seqtk, pigz, awscli, samtools

# ---- This project -----------------------------------------------------------
# The project dir is the folder that contains scripts/, data/, results/.
export PROJECT_DIR="${PROJECT_DIR:-$( cd "$_CFG_DIR/../.." && pwd )}"
export LOG_DIR="${LOG_DIR:-$PROJECT_DIR/logs/wgs}"
export NXF_WORK_BASE="${NXF_WORK_BASE:-$PROJECT_DIR/work/wgs}"          # huge, delete after a run
export RESULTS_BASE="${RESULTS_BASE:-$PROJECT_DIR/results/wgs}"
export SAMPLESHEET_DIR="$_CFG_DIR/samplesheets"
export SITE_CONFIG="$_CFG_DIR/conf/lobsang.config"
export ONCO_REFDATA_CONFIG="$_CFG_DIR/conf/oncoanalyser_refdata.config"  # written by 04

# ---- Test dataset: SEQC2 HCC1395 (breast cancer) tumour / HCC1395BL normal --
export DATASET="${DATASET:-HCC1395}"
export PATIENT="HCC1395"
export TUMOUR_ID="HCC1395T"
export NORMAL_ID="HCC1395BL"
export SEX="XX"
export CANCER_TYPE="BRCA"                                   # IntOGen code used by tumourevo
export DATA_DIR="${DATA_DIR:-$PROJECT_DIR/data/wgs_test/$DATASET}"
export RAW_DIR="$DATA_DIR/fastq_raw"                        # ~53x / ~55x, 191 GB
export SUB_DIR="$DATA_DIR/fastq_subsampled"                 # ~30x / ~20x
export TUMOUR_FRACTION="${TUMOUR_FRACTION:-0.56}"           # 53x -> ~30x
export NORMAL_FRACTION="${NORMAL_FRACTION:-0.36}"           # 55x -> ~20x

# ---- Container engine ------------------------------------------------------
if [[ -z "${CONTAINER_ENGINE:-}" ]]; then
    if   command -v apptainer   >/dev/null 2>&1; then CONTAINER_ENGINE=apptainer
    elif command -v singularity >/dev/null 2>&1; then CONTAINER_ENGINE=singularity
    elif command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then CONTAINER_ENGINE=docker
    else CONTAINER_ENGINE=none; fi
fi
export CONTAINER_ENGINE
export NXF_PROFILE="${NXF_PROFILE:-$CONTAINER_ENGINE}"

# ---- Nextflow runtime -------------------------------------------------------
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx8g}"
export NXF_SINGULARITY_CACHEDIR="$CONTAINER_DIR"
export NXF_APPTAINER_CACHEDIR="$CONTAINER_DIR"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$TMPDIR/apptainer_cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$TMPDIR/apptainer_tmp}"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
export NXF_HTTP_TIMEOUT="${NXF_HTTP_TIMEOUT:-30m}"

# ---- Helpers ---------------------------------------------------------------
log()  { printf '\033[1;34m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }
ok()   { printf '  \033[1;32mOK\033[0m    %s\n' "$*"; }
warn() { printf '  \033[1;33mWARN\033[0m  %s\n' "$*"; }
fail() { printf '  \033[1;31mFAIL\033[0m  %s\n' "$*"; }
die()  { printf '\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }

make_dirs() {
    mkdir -p "$CONTAINER_DIR" "$NXF_HOME" "$TMPDIR" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" \
             "$REF_BASE" "$LOG_DIR" "$NXF_WORK_BASE" "$RESULTS_BASE" "$SAMPLESHEET_DIR"
}

# Nextflow prefers JAVA_HOME over PATH; the system one on the server is Java 11.
set_java_home() {
    local prefix="${1:-$ENV_PREFIX}"
    if   [[ -x "$prefix/bin/java" ]];         then export JAVA_HOME="$prefix"
    elif [[ -x "$prefix/lib/jvm/bin/java" ]]; then export JAVA_HOME="$prefix/lib/jvm"
    else return 1; fi
    export JAVA_CMD="$JAVA_HOME/bin/java"
    export PATH="$JAVA_HOME/bin:$PATH"
}

activate_env() {
    command -v micromamba >/dev/null 2>&1 || die "micromamba not found in PATH"
    eval "$(micromamba shell hook --shell bash)"
    [[ -d "$ENV_PREFIX" ]] || die "env not found: $ENV_PREFIX  (run ./02_create_envs.sh)"
    micromamba activate "$ENV_PREFIX"
    set_java_home "$ENV_PREFIX" || die "no JDK inside $ENV_PREFIX"
}

activate_tools() {
    [[ -d "$TOOLS_PREFIX/bin" ]] || die "tools env not found: $TOOLS_PREFIX  (run ./02_create_envs.sh)"
    export PATH="$TOOLS_PREFIX/bin:$PATH"
}

# nextflow run wrapper: site config, per-run work dir, timestamped log, -resume
#   nf_run <run_name> <profile> <pipeline> [nextflow/pipeline args...]
nf_run() {
    local name="$1" profile="$2" pipeline="$3"; shift 3
    local logf="$LOG_DIR/${name}.$(date +%Y%m%d_%H%M%S).log"
    mkdir -p "$LOG_DIR" "$NXF_WORK_BASE/$name"
    log "pipeline : $pipeline"
    log "profile  : $profile"
    log "work dir : $NXF_WORK_BASE/$name"
    log "log      : $logf"
    echo
    nextflow -log "$logf" run "$pipeline" \
        -profile "$profile" \
        -c "$SITE_CONFIG" \
        -work-dir "$NXF_WORK_BASE/$name" \
        -resume \
        "$@"
}
