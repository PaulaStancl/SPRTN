#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Configuration for the WGS runs on the real SPRTN data: Novogene delivery
# X208SC25056159-Z01-F001, patient RJALS, tumour RJALS_Tm / normal RJALS_N.
#
#   source 00_config.sh
#
# Every other script in this folder sources this one - change paths HERE.
# Run settings (tool lists, versions, SEX, CANCER_TYPE) can be overridden from
# the shell, e.g.  SAREK_TOOLS=strelka,manta,ascat ./02_run_sarek.sh
# Project paths, the dataset name and the container engine can NOT - see below.
#
# One-time setup (envs, pipeline pulls, references) lives in ../wgs_test/
# (01-04) and is already done on the server; these scripts only reuse it.
#
# Why project paths and the container engine ignore the shell: other projects'
# configs (CHLOCK's env_setup/00_config.sh) export PROJECT_DIR, CONTAINER_ENGINE
# and NXF_PROFILE too. In a shell that had sourced one of them, this pipeline
# wrote into CHLOCK's folders under the apptainer profile.
# ---------------------------------------------------------------------------

_CFG_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export WGS_SCRIPTS="$_CFG_DIR"

# Everything these runs write derives from patient data - keep it private.
umask 077

# ---- Shared server layout (same as every other project) --------------------
export WORK_BASE="${WORK_BASE:-/common/WORK/pstancl}"
export ENV_ROOT="${ENV_ROOT:-$WORK_BASE/envs}"
export PROGRAMS_DIR="${PROGRAMS_DIR:-$WORK_BASE/PROGRAMI}"
export CONTAINER_DIR="${CONTAINER_DIR:-$WORK_BASE/singularity_cache}"   # ONE cache for all pipelines
export NXF_HOME="${NXF_HOME:-$PROGRAMS_DIR/nextflow}"                   # pulled pipelines live here
# Never $HOME or /tmp - and not the node-local TMPDIR a PBS job gets either.
export TMPDIR="$WORK_BASE/tmp"

# ---- Shared references (staged by ../wgs_test/04_download_references.sh) ---
export REF_BASE="${REF_BASE:-$WORK_BASE/references}"
export IGENOMES_BASE="${IGENOMES_BASE:-$REF_BASE/igenomes}"             # sarek: GATK.GRCh38
export HMF_REF_DIR="${HMF_REF_DIR:-$REF_BASE/hmf/oncoanalyser}"         # oncoanalyser: GRCh38_hmf
export VEP_CACHE="${VEP_CACHE:-$REF_BASE/vep_cache}"                    # tumourevo: VEP cache

# ---- Pipeline versions (same pins the test run passed with) ----------------
export SAREK_REV="${SAREK_REV:-3.10.0}"
export ONCOANALYSER_REV="${ONCOANALYSER_REV:-3.0.0}"
export TUMOUREVO_REV="${TUMOUREVO_REV:-738cb052fd51f47563ff0f96eaa4b2eb9d4b44f3}"   # dev, no release yet
export VEP_CACHE_VERSION="${VEP_CACHE_VERSION:-115}"     # tumourevo ships ensembl-vep 115.2

# ---- Envs (created by ../wgs_test/02_create_envs.sh) ------------------------
export NEXTFLOW_VERSION="${NEXTFLOW_VERSION:-26.04.6}"
export ENV_PREFIX="${ENV_PREFIX:-$ENV_ROOT/nextflow-$NEXTFLOW_VERSION}"
export TOOLS_PREFIX="${TOOLS_PREFIX:-$ENV_ROOT/wgs-tools}"              # bcftools for 04

# ---- This project -----------------------------------------------------------
# The project dir is the folder that contains scripts/, data/, results/.
# logs/, work/ and results/ are shared with ../wgs_test; DATASET keeps them apart.
export PROJECT_DIR="$( cd "$_CFG_DIR/../.." && pwd )"
export LOG_DIR="$PROJECT_DIR/logs/wgs"
export NXF_WORK_BASE="$PROJECT_DIR/work/wgs"                            # huge, delete after a run
export RESULTS_BASE="$PROJECT_DIR/results/wgs"
export SAMPLESHEET_DIR="$_CFG_DIR/samplesheets"
export SITE_CONFIG="$_CFG_DIR/conf/lobsang.config"
export ONCO_REFDATA_CONFIG="$_CFG_DIR/conf/oncoanalyser_refdata.config"  # copy of the one 04 wrote

# ---- Dataset: Novogene delivery, archived on /common/RAW (see md5check*.log) -
export DATASET="RJALS"
export PATIENT="RJALS"
export TUMOUR_ID="RJALS_Tm"
export NORMAL_ID="RJALS_N"
export DELIVERY_DIR="${DELIVERY_DIR:-/common/RAW/pstancl/MariaBoskovic/SPRTN/wgs/X208SC25056159-Z01-F001}"
export FASTQ_DIR="$DELIVERY_DIR/01.RawData"   # <sample>/<sample>_<library>_<flowcell>_L<n>_{1,2}.fq.gz

# Sex chromosomes of the patient: XY = male (confirmed by Paula 2026-09-22).
# sarek's ASCAT uses it, and because it is part of every task's inputs,
# changing it after sarek has started restarts sarek from scratch.
# oncoanalyser does not need it (AMBER infers sex, PURPLE reports it).
export SEX="${SEX:-XY}"

# IntOGen cancer-type code of the tumour, used by tumourevo for driver
# annotation (04 only): HCC = hepatocellular carcinoma, 204 driver genes in
# tumourevo's Compendium_Cancer_Genes.tsv. A code that is NOT in that table does
# not fail - tumourevo silently falls back to PANCANCER drivers.
export CANCER_TYPE="${CANCER_TYPE:-HCC}"

# ---- Container engine ------------------------------------------------------
# 'singularity' is preferred over 'apptainer' on purpose, even when the binary is
# just apptainer's compat symlink (lobsang: apptainer 1.4.2 provides both).
# oncoanalyser's local modules pick the prebuilt Galaxy SIF only when
# workflow.containerEngine == 'singularity'; under 'apptainer' they fall back to
# quay.io Docker images, and apptainer's OCI->SIF conversion fails on some of them
# (hmftools-esvee 2.0.1: FATAL "no descriptor found for reference ..."), which
# clearing the cache does not fix. The SIF route is a plain https download.
#
# Detected every time, never taken from the shell: CHLOCK's config exports
# CONTAINER_ENGINE=apptainer / NXF_PROFILE=apptainer under the same names.
#
# No system-wide singularity on this host? Use the micromamba-installed one
# (CHLOCK's env_setup/02b_install_apptainer.sh -> $ENV_ROOT/apptainer). It goes
# on PATH directly instead of being activated: activate_env switches to the
# nextflow env and would drop it again, leaving a PBS job with no engine.
if ! command -v singularity >/dev/null 2>&1 && [[ -x "$ENV_ROOT/apptainer/bin/singularity" ]]; then
    export PATH="$ENV_ROOT/apptainer/bin:$PATH"
fi
if   command -v singularity >/dev/null 2>&1; then CONTAINER_ENGINE=singularity
elif command -v apptainer   >/dev/null 2>&1; then CONTAINER_ENGINE=apptainer
elif command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then CONTAINER_ENGINE=docker
else CONTAINER_ENGINE=none; fi
export CONTAINER_ENGINE
export NXF_PROFILE="$CONTAINER_ENGINE"

# ---- Nextflow runtime -------------------------------------------------------
export NXF_OPTS="-Xms1g -Xmx8g"
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
    mkdir -p "$TMPDIR" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" \
             "$LOG_DIR" "$NXF_WORK_BASE" "$RESULTS_BASE" "$SAMPLESHEET_DIR"
}

# The data sits on /common/RAW, which the site config binds read-only into the
# containers. Fail now, not days into a run, if this host cannot read it.
check_data() {
    local s
    for s in "$TUMOUR_ID" "$NORMAL_ID"; do
        [[ -r "$FASTQ_DIR/$s" && -x "$FASTQ_DIR/$s" ]] \
            || die "cannot read $FASTQ_DIR/$s on $(hostname) - is /common/RAW mounted here?"
    done
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
    [[ -d "$ENV_PREFIX" ]] || die "env not found: $ENV_PREFIX  (run ../wgs_test/02_create_envs.sh)"
    micromamba activate "$ENV_PREFIX"
    set_java_home "$ENV_PREFIX" || die "no JDK inside $ENV_PREFIX"
}

activate_tools() {
    [[ -d "$TOOLS_PREFIX/bin" ]] || die "tools env not found: $TOOLS_PREFIX  (run ../wgs_test/02_create_envs.sh)"
    export PATH="$TOOLS_PREFIX/bin:$PATH"
}

# nextflow run wrapper: site config, resource caps, per-run launch + work dir,
# timestamped log
#   nf_run <run_name> <profile> <pipeline> [nextflow/pipeline args...]
#
# Resource caps are written for every launch to <run dir>/resources.config:
# inside a PBS job (qsub_*.sh) the job's allocation - PBS sets NCPUS, the qsub
# script sets JOB_MEMORY_GB because PBS exports no memory variable - and
# otherwise lobsang's 8 cpus / 400 GB. Nextflow must never schedule more than
# the job holds, or PBS kills the whole run. 12 GB are left for the Nextflow
# JVM itself (-Xmx8g) and the OS.
#
# Each run is launched from its own directory ($NXF_WORK_BASE/<run_name>), which
# holds that run's .nextflow/ history and its work/ dir. A bare -resume resumes
# the LAST run started in the launch directory, so with one shared launch dir,
# running oncoanalyser between two sarek attempts would make the second sarek
# attempt start from scratch. Every path passed in must therefore be absolute.
nf_run() {
    local name="$1" profile="$2" pipeline="$3"; shift 3
    local run_dir="$NXF_WORK_BASE/$name"
    local logf="$LOG_DIR/${name}.$(date +%Y%m%d_%H%M%S).log"
    local cpus=8 mem_gb=400
    if [[ -n "${PBS_JOBID:-}" ]]; then
        [[ -n "${NCPUS:-}" ]]         || die "PBS job without NCPUS - request cpus with -l select=1:ncpus=N"
        [[ -n "${JOB_MEMORY_GB:-}" ]] || die "PBS job without JOB_MEMORY_GB - set it in the qsub script, equal to mem="
        cpus="$NCPUS"; mem_gb="$JOB_MEMORY_GB"
    fi
    local task_gb=$(( mem_gb - 12 ))
    (( task_gb >= 16 )) || die "only $mem_gb GB for this run - request at least 32 GB"
    mkdir -p "$LOG_DIR" "$run_dir"
    cat > "$run_dir/resources.config" <<EOF
// Written by nf_run for the launch at $(date '+%F %T')${PBS_JOBID:+ in PBS job $PBS_JOBID}.
process {
    resourceLimits = [
        cpus  : $cpus,
        memory: ${task_gb}.GB,
        time  : 240.h
    ]
}
executor {
    cpus   = $cpus
    memory = ${task_gb}.GB
}
EOF
    log "pipeline : $pipeline"
    log "profile  : $profile"
    log "resources: $cpus cpus, $task_gb GB for tasks${PBS_JOBID:+   (PBS job $PBS_JOBID)}"
    log "run dir  : $run_dir   (history + work/)"
    log "log      : $logf"
    echo
    ( cd "$run_dir" && nextflow -log "$logf" run "$pipeline" \
        -profile "$profile" \
        -c "$SITE_CONFIG" \
        -c "$run_dir/resources.config" \
        -work-dir "$run_dir/work" \
        -resume \
        "$@" )
}
