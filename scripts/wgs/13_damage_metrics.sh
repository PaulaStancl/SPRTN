#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 13 - DNA damage metrics from the reads, independent of variant calling, tumour and normal
#      (for analysis/07_roq_damage.R: are Mutect2's low-ROQ calls deamination / oxidation?).
#
#   qsub qsub_damage.sh                         # PBS job (a few hours)
#   ./13_damage_metrics.sh                      # or directly
#   DAMAGE_REGIONS=chr2,chr20 ./13_damage_metrics.sh   # more of the genome (default chr20)
#   STEPS=reads ./13_damage_metrics.sh          # only the read-position step (after 07)
#   GATK_IMG=/path/to/gatk4.img ./13_damage_metrics.sh   # if the GATK4 image is not found
#
# Per sample, on a slice of sarek's recalibrated CRAM (DAMAGE_REGIONS, written to $TMPDIR as BAM
# and deleted afterwards - artefact rates need millions of bases, not the whole genome):
#   1. Picard CollectSequencingArtifactMetrics (GATK4 image from sarek's container cache):
#      pre_adapter_summary - TOTAL_QSCORE per substitution, ARTIFACT_NAME Deamination (C>T) and
#      OxoG (G>T); lower Q = more artefact. bait_bias is written too but means little for WGS.
#      Defaults: duplicates excluded, BQ >= 20, MAPQ >= 30.
#   2. damage_giv.py: GIV score (Chen et al. 2017), R1 vs R2 mismatch rate per type; > 1.5 = damaged.
# Then, if 07 has written read_level/sites.tsv (low-ROQ, high-ROQ and orientation-filtered
# C>T/G>A and G>T/C>A sites):
#   3. damage_read_positions.py on the tumour CRAM: allele, cycle, distance to fragment end,
#      BQ, MQ and mate of every read at those sites.
# Output: $RESULTS_BASE/analysis/qc_damage/sample_metrics/<sample>.* and read_level/reads_<sample>.tsv.gz
# Needs the qc-damage env: micromamba create -y -p $ENV_ROOT/qc-damage -f env_qc_damage.yml
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ./00_config.sh
activate_tools
umask 077

DAMAGE_REGIONS="${DAMAGE_REGIONS:-chr20}"
STEPS="${STEPS:-metrics,reads}"
PY="$ENV_ROOT/qc-damage/bin/python"
FASTA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
QC="$RESULTS_BASE/analysis/qc_damage"
SM="$QC/sample_metrics"; RL="$QC/read_level"
mkdir -p "$SM" "$RL" "$TMPDIR"
THREADS="${NCPUS:-4}"

[[ -x "$PY" ]] || die "no $PY - micromamba create -y -p $ENV_ROOT/qc-damage -f env_qc_damage.yml"
[[ -f "$FASTA.fai" ]] || die "missing $FASTA(.fai)"
command -v samtools >/dev/null || die "samtools not in $TOOLS_PREFIX/bin"

find_cram() {   # <sample>: sarek's recalibrated CRAM (else markduplicates)
    local s="$1" p hits
    for p in "$s.recal.cram" "$s.md.cram"; do
        mapfile -t hits < <(find "$RESULTS_BASE/sarek/$DATASET/preprocessing" -type f -name "$p" 2>/dev/null)
        (( ${#hits[@]} == 1 )) && { echo "${hits[0]}"; return 0; }
    done
    die "no $s.recal.cram / $s.md.cram under $RESULTS_BASE/sarek/$DATASET/preprocessing"
}

if [[ ",$STEPS," == *",metrics,"* ]]; then
    # sarek 3.10's GATK4 image is a Seqera community image, cached under its URL with / -> -
    # (community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-..-<hash>-data.img), so the
    # name has no "gatk4" in it: take the hash from the pulled sarek's own module. GATK_IMG overrides.
    # The module's location depends on the Nextflow version (assets/nf-core/sarek or
    # assets/.repos/...), hence find; sarek 3.10.0's hash is the fallback. Every lookup ends in
    # "|| true": under set -e + pipefail a grep with no match would end the script silently.
    IMG="${GATK_IMG:-}"; H=""
    if [[ -z "$IMG" ]]; then
        MOD=$(find "$NXF_HOME/assets" -path '*sarek*' -path '*gatk4/filtermutectcalls/main.nf' 2>/dev/null | head -1 || true)
        [[ -n "$MOD" ]] && H=$(grep -o 'sha256/[0-9a-f]*/[0-9a-f]\{64\}' "$MOD" 2>/dev/null | head -1 | sed 's#.*/##' || true)
        H="${H:-ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0}"
        IMG=$(find "$CONTAINER_DIR" -maxdepth 1 -name "*${H}*" \( -name '*.img' -o -name '*.sif' \) 2>/dev/null | head -1 || true)
    fi
    [[ -n "$IMG" ]] || IMG=$(find "$CONTAINER_DIR" -maxdepth 1 -iname '*gatk4*' ! -iname '*spark*' \( -name '*.img' -o -name '*.sif' \) | sort | tail -1)
    [[ -n "$IMG" && -f "$IMG" ]] || die "no GATK4 image in $CONTAINER_DIR (looked for sarek's ${H:-?} and *gatk4*) - set GATK_IMG=<path>"
    singularity exec "$IMG" gatk --version >/dev/null 2>&1 || die "$IMG has no working gatk"
    log "GATK4 image: $(basename "$IMG")"
    log "regions    : $DAMAGE_REGIONS"
    for s in "$TUMOUR_ID" "$NORMAL_ID"; do
        CRAM=$(find_cram "$s"); ok "$s: $CRAM"
        SLICE="$TMPDIR/qc_damage_${s}.bam"
        # shellcheck disable=SC2086
        samtools view -@ "$THREADS" -b -T "$FASTA" -o "$SLICE" "$CRAM" ${DAMAGE_REGIONS//,/ }
        samtools index "$SLICE"
        ok "slice: $(samtools view -c -@ "$THREADS" "$SLICE") reads"

        singularity exec -B "$(dirname "$FASTA")" -B "$TMPDIR" -B "$SM" "$IMG" \
            gatk --java-options "-Xmx8g -Djava.io.tmpdir=$TMPDIR" CollectSequencingArtifactMetrics \
                -I "$SLICE" -R "$FASTA" -O "$SM/$s" --TMP_DIR "$TMPDIR"
        ok "Picard: $SM/$s.pre_adapter_summary_metrics"

        "$PY" damage_giv.py --bam "$SLICE" --fasta "$FASTA" --regions "$DAMAGE_REGIONS" \
            --sample "$s" --out "$SM/$s.giv.tsv"
        ok "GIV: $SM/$s.giv.tsv"
        rm -f "$SLICE" "$SLICE.bai"
    done
    printf 'regions\t%s\ndate\t%s\n' "$DAMAGE_REGIONS" "$(date '+%F %T')" > "$SM/run_info.tsv"
fi

if [[ ",$STEPS," == *",reads,"* ]]; then
    if [[ -f "$RL/sites.tsv" ]]; then
        CRAM=$(find_cram "$TUMOUR_ID")
        "$PY" damage_read_positions.py --bam "$CRAM" --fasta "$FASTA" --sites "$RL/sites.tsv" \
            --sample "$TUMOUR_ID" --out "$RL/reads_${TUMOUR_ID}.tsv.gz"
        ok "read positions: $RL/reads_${TUMOUR_ID}.tsv.gz"
    else
        warn "no $RL/sites.tsv - run analysis/07_roq_damage.R first, then STEPS=reads ./13_damage_metrics.sh"
    fi
fi
log "done: $QC  ->  Rscript ../analysis/07_roq_damage.R"
