#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 09 - small BAM slices around chosen sites, to open in IGV on a laptop - no IGV
#      or display needed on the server (the alternative to 05_igv_snapshots.sh).
#
#   ./09_igv_slices.sh                                   # SPRTN Y117C, chr1:231347825
#   ./09_igv_slices.sh chr1:231347825 chr1:231351400     # any positions or ranges
#   PAD=1000 NAME=SPRTN ./09_igv_slices.sh chr1:231337104-231375416
#
# For each pipeline (sarek CRAMs, oncoanalyser BAMs) and sample (tumour, normal):
# reads within PAD bp of the sites -> <sample>_<pipeline>.slice.bam + .bai. sarek's
# CRAMs are written out as BAM, so IGV needs no reference FASTA to read them.
# Output: $RESULTS_BASE/igv_slices/<NAME>/
#   *.slice.bam(.bai)   4 slices (tumour/normal x sarek/oncoanalyser)
#   sites.tsv           the sites, IGV loci
#   regions.bed         the sliced regions (0-based, padded)
#   SPRTN_exons.bed     exon track
#   igv_batch.txt       IGV batch with RELATIVE paths: run it from this folder
#                       (IGV > Tools > Run Batch Script); PNGs -> igv_snapshots/
# The slices contain patient reads: copy them to a local, non-synced folder only.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ./00_config.sh
activate_tools
umask 077

PAD="${PAD:-500}"
NAME="${NAME:-SPRTN_Y117C}"
(( $# )) || set -- chr1:231347825                     # SPRTN c.350A>G, p.Tyr117Cys (rs527236213)
OUT="$RESULTS_BASE/igv_slices/$NAME"; mkdir -p "$OUT"
SAREK_FASTA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
command -v samtools >/dev/null || die "samtools not in $TOOLS_PREFIX/bin"
[[ -f "$SAREK_FASTA.fai" ]] || die "missing $SAREK_FASTA(.fai) - needed to decode sarek's CRAMs"

# ---- sites -> sites.tsv + regions.bed ------------------------------------------
: > "$OUT/sites.tsv"; : > "$OUT/regions.bed"
printf 'n\tlocus\tchrom\tstart\tend\n' > "$OUT/sites.tsv"
i=0
for l in "$@"; do
    l="${l//,/}"; chrom="${l%%:*}"; rng="${l#*:}"
    [[ "$chrom" != "$l" && -n "$rng" ]] || die "not a locus: $l (expected chr:pos or chr:start-end)"
    start="${rng%%-*}"; end="${rng##*-}"
    [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]] || die "not a locus: $l"
    i=$((i + 1))
    printf '%d\t%s\t%s\t%s\t%s\n' "$i" "$l" "$chrom" "$start" "$end" >> "$OUT/sites.tsv"
    lo=$(( start - 1 - PAD )); (( lo < 0 )) && lo=0
    printf '%s\t%d\t%d\tsite%d\n' "$chrom" "$lo" $(( end + PAD )) "$i" >> "$OUT/regions.bed"
done
log "$i site(s), +/- $PAD bp -> $OUT"

# ---- alignments: exactly one file per sample and pipeline ------------------------
find_aln() {   # <dir> <pattern...>  (first pattern with exactly one hit wins)
    local dir="$1" p hits; shift
    for p in "$@"; do
        mapfile -t hits < <(find "$dir" -type f -name "$p" 2>/dev/null)
        (( ${#hits[@]} == 1 )) && { echo "${hits[0]}"; return 0; }
        (( ${#hits[@]} > 1 )) && die "several $p under $dir: ${hits[*]}"
    done
    die "no alignment under $dir (looked for: $*)"
}

slice() {      # <pipeline> <sample> <alignment>
    local pl="$1" s="$2" aln="$3" out="$OUT/${2}_${1}.slice.bam" bed="$OUT/regions.bed" first
    ok "$pl $s: $aln"
    # chr1 or 1? Match the BED to the file's own contig names. awk reading to the end (no
    # head -1, no awk exit): stopping early, samtools would get SIGPIPE and pipefail + set -e end the script.
    first=$(samtools view -H "$aln" | awk '$1 == "@SQ" && !done { sub(/^SN:/, "", $2); print $2; done = 1 }')
    if [[ "$first" != chr* ]]; then sed 's/^chr//' "$OUT/regions.bed" > "$OUT/regions.nochr.bed"; bed="$OUT/regions.nochr.bed"; fi
    # -T: the reference to decode a CRAM (ignored for BAM); -M: regions may overlap
    samtools view -b -h -M -T "$SAREK_FASTA" -L "$bed" -o "$out" "$aln"
    samtools index "$out"
    ok "  -> $(basename "$out")  $(samtools view -c "$out") reads"
}

SP="$RESULTS_BASE/sarek/$DATASET/preprocessing"
OP="$RESULTS_BASE/oncoanalyser/$DATASET"
for s in "$TUMOUR_ID" "$NORMAL_ID"; do
    slice sarek        "$s" "$(find_aln "$SP" "$s.recal.cram" "$s.md.cram" "$s.recal.bam" "$s.md.bam")"
    slice oncoanalyser "$s" "$(find_aln "$OP" "$s.redux.bam" "$s.bam")"
done

# ---- IGV batch (relative paths) ---------------------------------------------------
cp "$WGS_SCRIPTS/igv_SPRTN_exons.bed" "$OUT/SPRTN_exons.bed" 2>/dev/null || true
{
    echo "new"
    echo "genome hg38"
    for s in "$TUMOUR_ID" "$NORMAL_ID"; do for pl in sarek oncoanalyser; do
        echo "load ${s}_${pl}.slice.bam"
    done; done
    [[ -f "$OUT/SPRTN_exons.bed" ]] && echo "load SPRTN_exons.bed"
    echo "snapshotDirectory igv_snapshots"
    echo "maxPanelHeight 400"
    echo "preference SAM.SHOW_SOFT_CLIPPED false"
    tail -n +2 "$OUT/sites.tsv" | while IFS=$'\t' read -r n locus chrom start end; do
        if [[ "$start" == "$end" ]]; then echo "goto $chrom:$(( start - 60 ))-$(( end + 60 ))"; else echo "goto $chrom:$start-$end"; fi
        echo "sort base $chrom:$start"
        echo "collapse"
        printf 'snapshot %02d_%s_%s.png\n' "$n" "$NAME" "${locus//[:-]/_}"
    done
    echo "exit"
} > "$OUT/igv_batch.txt"

ls -lh "$OUT"
log "copy to your laptop (a local folder, NOT OneDrive/iCloud - these are patient reads):"
echo "  rsync -av pstancl@ssi-access.chem.pmf.hr:$OUT/ ~/igv_RJALS/$NAME/"
echo "  then IGV (genome hg38) > Tools > Run Batch Script > ~/igv_RJALS/$NAME/igv_batch.txt"
echo "  or File > Load from File > the four *.slice.bam, and go to a locus from sites.tsv"
