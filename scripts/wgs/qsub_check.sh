#!/bin/bash
# ---------------------------------------------------------------------------
# PBS job (a few minutes): does a q2 node have what the pipeline jobs need?
# Run it once before the first real submission.
#
#   cd scripts/wgs && qsub qsub_check.sh        # then read check_q2.o<jobid>
#
# Changes nothing. Every line should read OK; a MISSING one means the real jobs
# would fail on that node.
# ---------------------------------------------------------------------------
#PBS -q q2
# q2 rejects jobs with fewer than 2 cpus. (No comments after a #PBS directive:
# PBS reads them as part of the value.)
#PBS -l select=1:ncpus=2:mem=2gb
#PBS -l walltime=00:15:00
#PBS -N check_q2
#PBS -j oe

cd "$PBS_O_WORKDIR" && [[ -f 00_config.sh ]] \
    || { echo "ERROR: submit from scripts/wgs/" >&2; exit 1; }
PBS_TMPDIR="${TMPDIR:-unset}"
source ./00_config.sh

row() { printf '  %-8s %-24s %s\n' "$1" "$2" "$3"; }
yes_no() { if eval "$2"; then row OK "$1" "$3"; else row MISSING "$1" "$4"; fi; }

echo "node $(hostname): $(nproc) cores, $(free -g | awk '/^Mem:/{print $2}') GB RAM   (this job: NCPUS=${NCPUS:-unset}, PBS TMPDIR=$PBS_TMPDIR)"
yes_no micromamba        'command -v micromamba >/dev/null'  "$(command -v micromamba)"  "not in PATH in a batch job - activate_env will stop"
yes_no singularity       'command -v singularity >/dev/null' "$(command -v singularity)" "oncoanalyser refuses to start without it"
yes_no "nextflow env"    '[[ -d $ENV_PREFIX ]]'              "$ENV_PREFIX"               "$ENV_PREFIX"
yes_no "tools env"       '[[ -d $TOOLS_PREFIX/bin ]]'        "$TOOLS_PREFIX"             "tumourevo (bcftools) needs it"
yes_no "RAW FASTQ"       '[[ -r $FASTQ_DIR/$TUMOUR_ID && -r $FASTQ_DIR/$NORMAL_ID ]]' "$FASTQ_DIR" "/common/RAW not mounted here"
yes_no "WORK writable"   'touch "$TMPDIR/.check_q2.$$" 2>/dev/null && rm -f "$TMPDIR/.check_q2.$$"' "$TMPDIR" "cannot write $TMPDIR"
yes_no "container cache" '[[ -d $CONTAINER_DIR ]]'           "$CONTAINER_DIR"            "$CONTAINER_DIR"
# Actually start a cached image: a micromamba-installed singularity only works
# where the kernel allows unprivileged user namespaces, which can differ by node.
img=$(find "$CONTAINER_DIR" -maxdepth 1 \( -name '*.img' -o -name '*.sif' \) 2>/dev/null | head -1)
yes_no "run container"   '[[ -n $img ]] && singularity exec "$img" true 2>/dev/null' "$(basename "${img:-none}")" "singularity cannot start a container on this node"
yes_no internet          'curl -sfI -m 15 https://github.com >/dev/null' "github.com reachable" "no internet - add  export NXF_OFFLINE=true  to the qsub scripts"
echo "  engine   : $CONTAINER_ENGINE (profile $NXF_PROFILE)"
