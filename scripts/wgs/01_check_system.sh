#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01 - Pre-flight check on lobsang. Read-only: installs nothing.
# ---------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh

problems=0

echo; log "1. micromamba"
if command -v micromamba >/dev/null 2>&1; then ok "micromamba $(micromamba --version)"
else fail "micromamba not on PATH"; problems=$((problems+1)); fi

echo; log "2. Container engine (all three pipelines run their tools in containers)"
case "$CONTAINER_ENGINE" in
    apptainer)   ok "apptainer $(apptainer version 2>/dev/null)" ;;
    singularity) ok "singularity $(singularity --version 2>/dev/null)" ;;
    docker)      ok "docker $(docker --version 2>/dev/null)" ;;
    none)        fail "no apptainer / singularity / docker"; problems=$((problems+1)) ;;
esac

echo; log "3. Nextflow env"
if [[ -x "$ENV_PREFIX/bin/nextflow" ]]; then ok "exists: $ENV_PREFIX"
else warn "missing: $ENV_PREFIX  -> ./02_create_envs.sh"; fi
[[ -d "$TOOLS_PREFIX/bin" ]] && ok "exists: $TOOLS_PREFIX" || warn "missing: $TOOLS_PREFIX  -> ./02_create_envs.sh"

echo; log "4. CPU / RAM  (compare with conf/lobsang.config)"
ncpu=$(nproc)
memgb=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) / 1024 / 1024 ))
cfg_cpu=$(grep -oE 'cpus *: *[0-9]+' "$SITE_CONFIG" | grep -oE '[0-9]+' | head -1)
cfg_mem=$(grep -oE 'memory *: *[0-9]+' "$SITE_CONFIG" | grep -oE '[0-9]+' | head -1)
ok "machine: $ncpu cpus, $memgb GB RAM"
[[ "$cfg_cpu" -le "$ncpu" ]]  && ok "config cpus $cfg_cpu <= $ncpu" || { fail "config cpus $cfg_cpu > $ncpu - edit $SITE_CONFIG"; problems=$((problems+1)); }
[[ "$cfg_mem" -le "$memgb" ]] && ok "config memory ${cfg_mem} GB <= $memgb GB" || { fail "config memory ${cfg_mem} GB > $memgb GB - edit $SITE_CONFIG"; problems=$((problems+1)); }
[[ "$memgb" -ge 96 ]] || warn "oncoanalyser recommends 72-96 GB for alignment / SAGE / ESVEE"

echo; log "5. Disk  (budget ~1.5 TB for this test)"
cat <<TXT
        references : sarek iGenomes ~47 GB, oncoanalyser ~31 GB, VEP cache ~28 GB
        FASTQ      : raw 191 GB + subsampled ~70 GB (both kept)
        work dirs  : sarek ~500 GB, oncoanalyser ~400 GB (delete after each run)
TXT
for d in "$WORK_BASE" "$PROJECT_DIR"; do
    avail=$(df -BG --output=avail "$d" 2>/dev/null | tail -1 | tr -dc '0-9')
    [[ "${avail:-0}" -ge 1500 ]] && ok "$d : ${avail} GB free" || warn "$d : ${avail:-?} GB free"
done

echo; log "6. Network"
urls=(
    "https://github.com/nf-core/sarek"
    "https://ngi-igenomes.s3.amazonaws.com/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai"
    "https://data.oncoanalyser.com/r2/reference/dist/v1/genomes/GRCh38_hmf/26.1/samtools_index-1.16/GRCh38_masked_exclusions_alts_hlas.fasta.fai"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR789/009/SRR7890829/"
)
for u in "${urls[@]}"; do
    code=$(curl -sIL -o /dev/null -w '%{http_code}' --max-time 20 "$u" 2>/dev/null)
    host=${u#https://}; host=${host%%/*}
    [[ "$code" =~ ^[23] ]] && ok "$host (HTTP $code)" || { fail "$host NOT reachable (HTTP ${code:-none})"; problems=$((problems+1)); }
done

echo; log "7. Layout"
cat <<TXT
  project      : $PROJECT_DIR
  data         : $DATA_DIR
  results      : $RESULTS_BASE
  work dirs    : $NXF_WORK_BASE
  nextflow env : $ENV_PREFIX
  containers   : $CONTAINER_DIR
  iGenomes     : $IGENOMES_BASE
  HMF refs     : $HMF_REF_DIR
  VEP cache    : $VEP_CACHE
TXT

echo
[[ "$problems" -eq 0 ]] && log "All good. Next:  ./02_create_envs.sh" \
                        || log "$problems problem(s) above - fix before continuing"
