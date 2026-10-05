#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# IGV snapshots (PNG) of tumour + normal reads at chosen regions, once for the
# sarek CRAMs and once for the oncoanalyser BAMs.
#
#   ./05_igv_snapshots.sh                          # every region in igv_regions.tsv
#   ./05_igv_snapshots.sh chr1:231351254-231351621 # or loci given here instead
#   PIPELINES=sarek ./05_igv_snapshots.sh          # only one of the two
#
# Run it in an interactive job (qsub -I ... mem=16gb), not on the login node.
# Output: $RESULTS_BASE/igv/<pipeline>/<region>.png, plus the IGV batch files.
#
# IGV has no headless mode: it is run under xvfb-run (a virtual display). An IGV
# already in one of the envs is used; otherwise it is downloaded once.
# Each pipeline is shown on its own reference: a CRAM can only be decoded with
# the FASTA it was written against, and the HMF FASTA is masked, unlike GATK's.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/00_config.sh"

PIPELINES="${PIPELINES:-sarek oncoanalyser}"
REGIONS="${REGIONS:-$WGS_SCRIPTS/igv_regions.tsv}"
TRACKS="${TRACKS:-$WGS_SCRIPTS/igv_SPRTN_exons.bed}"        # extra feature tracks, space-separated
OUT_IGV="$RESULTS_BASE/igv"
IGV_VERSION="${IGV_VERSION:-2.19.8}"

SAREK_FASTA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

# ---- IGV and a display -------------------------------------------------------
# First IGV found: $IGV (set it to force one), any env in $ENV_ROOT (bioconda's
# `igv`), an unzipped IGV in $PROGRAMS_DIR; only if none, download IGV once.
find_first() { local f; for f in "$@"; do [[ -x "$f" ]] && { echo "$f"; return; }; done; }
IGV="${IGV:-$(find_first "$ENV_ROOT"/*/bin/igv "$ENV_ROOT"/*/bin/igv.sh "$PROGRAMS_DIR"/IGV_Linux_*/igv.sh)}"
if [[ -z "$IGV" ]]; then
    log "no IGV in $ENV_ROOT or $PROGRAMS_DIR - downloading $IGV_VERSION to $PROGRAMS_DIR"
    zip="$TMPDIR/IGV_Linux_${IGV_VERSION}_WithJava.zip"
    wget -q -O "$zip" "https://data.broadinstitute.org/igv/projects/downloads/${IGV_VERSION%.*}/IGV_Linux_${IGV_VERSION}_WithJava.zip" \
        || die "download failed (no internet on this node? run once on the login node)"
    unzip -q -o "$zip" -d "$PROGRAMS_DIR" && rm -f "$zip"
    IGV="$PROGRAMS_DIR/IGV_Linux_${IGV_VERSION}/igv.sh"
    [[ -x "$IGV" ]] || die "unexpected zip layout - no $IGV"
fi
ok "IGV: $IGV"

# xvfb-run from PATH or from any env (conda-forge ships it with Xvfb).
XVFB_RUN="$(command -v xvfb-run || find_first "$ENV_ROOT"/*/bin/xvfb-run || true)"
if   [[ -n "$XVFB_RUN" ]]; then
    export PATH="$(dirname "$XVFB_RUN"):$PATH"                    # its Xvfb sits next to it
    RUN=("$XVFB_RUN" --auto-servernum --server-args="-screen 0 1600x1200x24"); ok "display: $XVFB_RUN"
elif [[ -n "${DISPLAY:-}" ]]; then RUN=(); warn "no xvfb-run - using your display $DISPLAY (X forwarding)"
else die "neither xvfb-run nor a DISPLAY on $(hostname) - ssh -X and rerun, or install xvfb-run"; fi

# ---- regions -----------------------------------------------------------------
# name<TAB>locus lines; loci given as arguments are named after the locus.
regions=()
if (( $# )); then
    for l in "$@"; do regions+=("$(tr ':-' '__' <<<"${l//,/}")	$l"); done
else
    [[ -r "$REGIONS" ]] || die "no regions file $REGIONS"
    while IFS= read -r line; do
        [[ -z "$line" || "$line" == \#* ]] || regions+=("$line")
    done < "$REGIONS"
fi
(( ${#regions[@]} )) || die "no regions"
log "${#regions[@]} region(s)"

# ---- inputs per pipeline -------------------------------------------------------
# Exactly one file per sample; the first pattern that matches wins.
find_aln() {   # <dir> <sample> <pattern>...
    local dir="$1" s="$2" p f hits; shift 2
    for p in "$@"; do
        hits=()
        while IFS= read -r f; do hits+=("$f"); done < <(find "$dir" -type f -name "${p//@/$s}" 2>/dev/null)
        if   (( ${#hits[@]} == 1 )); then echo "${hits[0]}"; return; fi
        if   (( ${#hits[@]} > 1 ));  then die "several ${p//@/$s} under $dir: ${hits[*]}"; fi
    done
    die "no alignment for $s under $dir (looked for: $*)"
}

has_index() { [[ -e "$1.crai" || -e "$1.bai" || -e "${1%.*}.bai" || -e "${1%.*}.crai" ]]; }

inputs() {     # sets FASTA, TUMOUR_ALN, NORMAL_ALN for pipeline $1
    case "$1" in
        sarek)
            local d="$RESULTS_BASE/sarek/$DATASET/preprocessing"
            FASTA="$SAREK_FASTA"
            TUMOUR_ALN=$(find_aln "$d" "$TUMOUR_ID" "@.recal.cram" "@.md.cram" "@.recal.bam" "@.md.bam")
            NORMAL_ALN=$(find_aln "$d" "$NORMAL_ID" "@.recal.cram" "@.md.cram" "@.recal.bam" "@.md.bam") ;;
        oncoanalyser)
            local d="$RESULTS_BASE/oncoanalyser/$DATASET"
            FASTA=$(sed -n 's/^ *fasta *= *"\(.*\)".*/\1/p' "$ONCO_REFDATA_CONFIG" | head -1)
            # alignments/dna/<sample>.redux.bam; never the RNA ones (wgts mode).
            TUMOUR_ALN=$(find_aln "$d" "$TUMOUR_ID" "@.redux.bam" "@.bam")
            NORMAL_ALN=$(find_aln "$d" "$NORMAL_ID" "@.redux.bam" "@.bam") ;;
        *) die "unknown pipeline $1 (sarek or oncoanalyser)" ;;
    esac
    [[ -r "$FASTA" && -r "$FASTA.fai" ]] || die "$1: reference FASTA or its .fai missing: $FASTA"
    local f; for f in "$TUMOUR_ALN" "$NORMAL_ALN"; do has_index "$f" || die "no index next to $f"; done
}

# ---- one IGV session per pipeline ------------------------------------------------
for pl in $PIPELINES; do
    inputs "$pl"
    od="$OUT_IGV/$pl"; mkdir -p "$od"
    log "$pl"; ok "reference: $FASTA"; ok "tumour   : $TUMOUR_ALN"; ok "normal   : $NORMAL_ALN"

    batch="$od/igv_batch.txt"
    {
        echo "new"
        echo "genome $FASTA"
        echo "load $TUMOUR_ALN name=${TUMOUR_ID}_tumour_$pl"
        echo "load $NORMAL_ALN name=${NORMAL_ID}_normal_$pl"
        for t in $TRACKS; do [[ -r "$t" ]] && echo "load $t"; done
        echo "snapshotDirectory $od"
        echo "maxPanelHeight 1000"
        echo "preference SAM.SHOW_SOFT_CLIPPED false"
        for r in "${regions[@]}"; do
            name="${r%%	*}"; locus="${r#*	}"
            echo "goto $locus"
            echo "sort base"
            echo "collapse"
            echo "snapshot ${name}.png"
        done
        echo "exit"
    } > "$batch"

    log "running IGV ($pl) - batch: $batch"
    ${RUN[@]+"${RUN[@]}"} "$IGV" --batch "$batch" > "$od/igv.log" 2>&1 \
        || { tail -20 "$od/igv.log"; die "IGV failed for $pl - see $od/igv.log"; }
    n=$(find "$od" -maxdepth 1 -name '*.png' -newer "$batch" | wc -l | tr -d " ")
    (( n == ${#regions[@]} )) && ok "$n snapshot(s) in $od" \
        || { warn "$n of ${#regions[@]} snapshots written for $pl - see $od/igv.log"; }
done
