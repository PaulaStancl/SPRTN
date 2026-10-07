#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02 - nf-core/sarek: FASTQ -> bwa-mem2 -> markdup -> BQSR -> somatic calling
#      RJALS tumour/normal WGS on GATK.GRCh38. Run inside screen/tmux - on 8
#      cores at this depth expect about a week (see README).
#
#   ./02_run_sarek.sh
#   SAREK_TOOLS=strelka,manta,ascat ./02_run_sarek.sh    # faster, no Mutect2
#   SAREK_STEP=variant_calling SAREK_TOOLS=muse,msisensorpro ./02_run_sarek.sh
#   SAREK_STEP=annotate ./02_run_sarek.sh                # VEP on every caller's VCF
#
# SAREK_STEP=variant_calling adds callers to a finished run without realigning:
# it starts from that run's recalibrated CRAMs (csv/recalibrated.csv) and writes
# to its own outdir (<dataset>_vc) and launch dir, so the main run's MultiQC,
# csv/ and resume history stay untouched.
#
# SAREK_STEP=annotate runs only sarek's annotation (VEP $SAREK_VEP_CACHE_VERSION, the version
# sarek 3.10 ships) on the somatic VCFs of the finished runs: the bcftools-normalised SNV/indel
# VCFs from 06 (Mutect2, Strelka2 snvs + indels, MuSE), so they match the analysis, plus Manta's
# somatic SVs. Output: <dataset>_annotate/annotation/<caller>/<pair>/*_VEP.ann.vcf.gz.
# Needs the VEP cache: VEP_CACHE_VERSION=116 ../wgs_test/04_download_references.sh vep
# VEP runs with sarek's default arguments minus --filter_common, which would drop somatic
# calls that sit on a common germline SNP position.
#
# Tools: mutect2 + ascat are what tumourevo (04) consumes; strelka + manta are
# the second SNV/indel caller and the SV caller. Mutect2 is the slowest step.
# No VEP in the mapping run - tumourevo annotates with its own VEP 115; SAREK_STEP=annotate
# adds VEP 116 for all sarek callers afterwards.
# Re-running resumes (-resume) from the last finished task.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_env
make_dirs
check_data

SAREK_TOOLS="${SAREK_TOOLS:-mutect2,strelka,manta,ascat}"
SAREK_STEP="${SAREK_STEP:-mapping}"
case "$SAREK_STEP" in
    mapping)
        RUN="sarek_${DATASET}"
        SHEET="$SAMPLESHEET_DIR/sarek_${DATASET}.csv"
        OUT="$RESULTS_BASE/sarek/$DATASET"
        [[ -f "$SHEET" ]] || die "missing $SHEET - set SEX in 00_config.sh, then run ./01_make_samplesheets.sh" ;;
    variant_calling)
        RUN="sarek_${DATASET}_vc"
        SHEET="$RESULTS_BASE/sarek/$DATASET/csv/recalibrated.csv"
        OUT="$RESULTS_BASE/sarek/${DATASET}_vc"
        [[ -f "$SHEET" ]] || die "missing $SHEET - the main sarek run (SAREK_STEP=mapping) has to finish first" ;;
    annotate)
        RUN="sarek_${DATASET}_annotate"
        SHEET="$SAMPLESHEET_DIR/sarek_${DATASET}_annotate.csv"
        OUT="$RESULTS_BASE/sarek/${DATASET}_annotate"
        SAREK_TOOLS=vep
        PAIR="${TUMOUR_ID}_vs_${NORMAL_ID}"
        NORM="$RESULTS_BASE/sarek/$DATASET/normalized_bcftools/$PAIR"
        [[ -d "$NORM" ]] || die "missing $NORM - run ./06_normalize_vcfs.sh first"
        [[ -d "$VEP_CACHE/homo_sapiens/${SAREK_VEP_CACHE_VERSION}_GRCh38" ]] \
            || die "VEP cache ${SAREK_VEP_CACHE_VERSION} missing - VEP_CACHE_VERSION=${SAREK_VEP_CACHE_VERSION} ../wgs_test/04_download_references.sh vep"
        # patient,sample,variantcaller,vcf - one row per VCF; sarek names outputs by caller
        {
            echo "patient,sample,variantcaller,vcf"
            for f in "$NORM"/*.norm.vcf.gz; do
                case "$(basename "$f")" in
                    *.mutect2.*) c=mutect2 ;; *.strelka.*) c=strelka ;; *.muse.*) c=muse ;;
                    *.freebayes.*) c=freebayes ;; *lofreq*) c=lofreq ;; *) continue ;;
                esac
                echo "$PATIENT,$PAIR,$c,$f"
            done
            sv=$(find "$RESULTS_BASE/sarek/$DATASET/variant_calling/manta/$PAIR" -name "${PAIR}.manta.somatic_sv.vcf.gz" 2>/dev/null | head -1)
            [[ -n "$sv" ]] && echo "$PATIENT,$PAIR,manta,$sv"
        } > "$SHEET"
        (( $(wc -l < "$SHEET") > 1 )) || die "no VCFs found for $SHEET"
        log "annotate sheet: $SHEET"; sed 's/^/    /' "$SHEET" ;;
    *)  die "SAREK_STEP must be mapping, variant_calling or annotate, not '$SAREK_STEP'" ;;
esac
# A sheet written before SEX was changed would start a run with the wrong sex -
# and fixing that later restarts every task.
if [[ "$SAREK_STEP" != annotate ]]; then
    sheet_sex=$(awk -F, 'NR > 1 { print $2 }' "$SHEET" | sort -u | tr '\n' ' ')
    [[ "$sheet_sex" == "$SEX " ]] || die "$SHEET has sex '$sheet_sex' but 00_config.sh says '$SEX' - re-run ./01_make_samplesheets.sh"
fi
[[ -d "$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem2Index" ]] \
    || die "iGenomes not staged - run ../wgs_test/04_download_references.sh sarek"

# Per-launch tuning, written next to the run:
#
# 1. Time. sarek's default is 8 h per task (16 h / 32 h for its process_medium /
#    process_high labels), and GATK4_MARKDUPLICATES sets only cpus and memory, so
#    it inherits 8 h. On this 90x pair it needs more, and Nextflow kills a local
#    task that passes its time (SIGTERM -> exit 143); that killed the runs of
#    2026-09-22 and 2026-09-23 at exactly 8 h. Inside a PBS job the walltime is
#    the only limit that should apply, so lift every task to it.
# 2. bwa-mem2 asks 24 cpus per chunk, so in a job of 32+ cpus only one chunk fits
#    and the rest idle through the whole alignment. Give it half the job instead.
#
# Neither changes a task's command line, so a resume still reuses finished tasks.
TUNING="$NXF_WORK_BASE/$RUN/sarek_tuning.config"
mkdir -p "$(dirname "$TUNING")"
{
    echo "// Written by 02_run_sarek.sh at $(date '+%F %T')."
    echo "process {"
    echo "    time = 240.h"
    for _lab in process_single process_low process_medium process_high process_long; do
        echo "    withLabel: $_lab { time = 240.h }"
    done
    if (( ${NCPUS:-0} >= 32 )); then
        echo "    withName: 'BWAMEM2_MEM' { cpus = $(( NCPUS / 2 )) }"
    fi
    echo "}"
} > "$TUNING"
log "task time: 240 h (sarek's default is 8 h - it killed markdup twice)"
(( ${NCPUS:-0} >= 32 )) && log "bwa-mem2 : $(( NCPUS / 2 )) cpus per chunk, two chunks at a time"
log "sex      : $SEX"
log "step     : $SAREK_STEP   (sheet: $SHEET)"
log "tools    : $SAREK_TOOLS"

# MSIsensor-pro: the microsatellite list from iGenomes was made with an older msisensor-pro,
# whose format sarek's version no longer reads (it loads 0 sites, then aborts: "Same reference
# genome file should be used in both scan and msi steps"). 07_msisensorpro_scan.sh builds the
# list with sarek's own container; it is used automatically once it exists.
EXTRA=()
MSISENSORPRO_SCAN="${MSISENSORPRO_SCAN:-$REF_BASE/msisensorpro/Homo_sapiens_assembly38.msisensorpro_scan.list}"
if [[ ",$SAREK_TOOLS," == *",msisensorpro,"* ]]; then
    [[ -s "$MSISENSORPRO_SCAN" ]] || die "missing $MSISENSORPRO_SCAN - run ./07_msisensorpro_scan.sh first (iGenomes' list does not work with sarek's msisensor-pro)"
    EXTRA+=(--msisensorpro_scan "$MSISENSORPRO_SCAN")
    log "msisensor: $MSISENSORPRO_SCAN"
fi
if [[ "$SAREK_STEP" == annotate ]]; then
    # In a params file, not on the command line: there Nextflow 26 passes "--x true" as the
    # string "true" (sarek's schema rejects it), and a value starting with "--" (the VEP
    # arguments) can be taken for options. download_cache is false by default.
    VEP_PARAMS="$NXF_WORK_BASE/$RUN/vep_params.yml"
    cat > "$VEP_PARAMS" <<EOF
vep_cache: "$VEP_CACHE"
vep_cache_version: "$SAREK_VEP_CACHE_VERSION"
vep_include_fasta: true
vep_custom_args: "--everything --per_gene --total_length --offline --format vcf"
EOF
    EXTRA+=(-params-file "$VEP_PARAMS")
    log "VEP      : ${SAREK_VEP_CACHE_VERSION}_GRCh38 from $VEP_CACHE (no --filter_common)"
fi
log "outdir   : $OUT"

nf_run "$RUN" "$NXF_PROFILE" nf-core/sarek -r "$SAREK_REV" \
    --input "$SHEET" \
    --step "$SAREK_STEP" \
    --outdir "$OUT" \
    --genome GATK.GRCh38 \
    --igenomes_base "$IGENOMES_BASE" \
    --aligner bwa-mem2 \
    --tools "$SAREK_TOOLS" \
    ${EXTRA[@]+"${EXTRA[@]}"} \
    -c "$TUNING"

log "sarek done: $OUT"
if [[ "$SAREK_STEP" == mapping ]]; then
    log "Next:  ./04_run_tumourevo.sh   (then: rm -rf $NXF_WORK_BASE/$RUN)"
elif [[ "$SAREK_STEP" == annotate ]]; then
    log "Annotated VCFs: $OUT/annotation/<caller>/${TUMOUR_ID}_vs_${NORMAL_ID}/   VEP summaries: $OUT/reports/EnsemblVEP/"
    log "Work dir no longer needed: rm -rf $NXF_WORK_BASE/$RUN"
else
    log "Work dir no longer needed: rm -rf $NXF_WORK_BASE/$RUN"
fi
