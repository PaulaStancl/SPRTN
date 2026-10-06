#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 10 - igv-reports: interactive HTML pages (igv.js) of the reads at SPRTN, tumour and
#      normal from both pipelines. Runs on the server; no IGV, no display needed.
#
#   micromamba create -y -p /common/WORK/pstancl/envs/igvreports -c conda-forge -c bioconda igv-reports   # once
#   ./09_igv_slices.sh                                            # slices for Y117C
#   NAME=SPRTN_gene PAD=200 ./09_igv_slices.sh <gene + exon loci>   # slices for the gene (see 09)
#   ./10_igv_reports.sh
#
# Output: $RESULTS_BASE/igv_reports/
#   SPRTN_Y117C.html   the Y117C site (chr1:231347825 A>G), as called in the normal by
#                      Strelka2 germline - its genotype and read counts are in the table
#   SPRTN_gene.html    the whole gene and each exon (regions from 09's sites.tsv)
# Tracks: 09's slices (RJALS_Tm / RJALS_N x sarek / oncoanalyser) + the SPRTN exon track.
# The HTML files embed the patient's reads: keep them local, outside synced folders.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ./00_config.sh
activate_tools                                       # bcftools
umask 077

IGVR="${IGVR:-$ENV_ROOT/igvreports/bin/create_report}"
[[ -x "$IGVR" ]] || die "no create_report at $IGVR - micromamba create -y -p $ENV_ROOT/igvreports -c conda-forge -c bioconda igv-reports"
FASTA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
SLICES="$RESULTS_BASE/igv_slices"
OUT="$RESULTS_BASE/igv_reports"; mkdir -p "$OUT"
GERM="$RESULTS_BASE/sarek/$DATASET/variant_calling/strelka/$NORMAL_ID/$NORMAL_ID.strelka.variants.vcf.gz"
EXONS="$WGS_SCRIPTS/igv_SPRTN_exons.bed"

tracks() {     # <slice dir> -> the four slice BAMs (+ exon track), tumour first
    local d="$1" s pl
    for s in "$TUMOUR_ID" "$NORMAL_ID"; do for pl in sarek oncoanalyser; do
        [[ -f "$d/${s}_${pl}.slice.bam" ]] || die "missing $d/${s}_${pl}.slice.bam - run ./09_igv_slices.sh first"
        echo "$d/${s}_${pl}.slice.bam"
    done; done
    [[ -f "$EXONS" ]] && grep -v '^track' "$EXONS" > "$OUT/SPRTN_exons.bed" && echo "$OUT/SPRTN_exons.bed"
    return 0
}

# ---- 1. Y117C: the germline call itself as the site (genotype, AD, DP in the table) -------
d="$SLICES/SPRTN_Y117C"
if [[ -d "$d" ]]; then
    [[ -f "$GERM" ]] || die "missing $GERM"
    bcftools view -r chr1:231347825 "$GERM" -Oz -o "$OUT/SPRTN_Y117C.site.vcf.gz"
    bcftools index -t -f "$OUT/SPRTN_Y117C.site.vcf.gz"
    (( $(bcftools view -H "$OUT/SPRTN_Y117C.site.vcf.gz" | wc -l) )) || die "Y117C not in $GERM"
    mapfile -t tr < <(tracks "$d")
    "$IGVR" "$OUT/SPRTN_Y117C.site.vcf.gz" --fasta "$FASTA" --tracks "${tr[@]}" \
        --flanking 100 --sample-columns GT AD DP --title "SPRTN p.Tyr117Cys (c.350A>G, rs527236213) - RJALS" \
        --output "$OUT/SPRTN_Y117C.html"
    ok "$OUT/SPRTN_Y117C.html"
else warn "no slices in $d - ./09_igv_slices.sh first; Y117C report skipped"; fi

# ---- 2. the whole gene + each exon, from 09's sites.tsv ---------------------------------------
d="$SLICES/SPRTN_gene"
if [[ -d "$d" ]]; then
    awk -F'\t' 'NR > 1 { printf "%s\t%d\t%d\tsite%s %s\n", $3, $4 - 1, $5, $1, $2 }' "$d/sites.tsv" > "$OUT/SPRTN_gene.sites.bed"
    mapfile -t tr < <(tracks "$d")
    "$IGVR" "$OUT/SPRTN_gene.sites.bed" --fasta "$FASTA" --tracks "${tr[@]}" \
        --flanking 50 --title "SPRTN gene and exons - RJALS" --output "$OUT/SPRTN_gene.html"
    ok "$OUT/SPRTN_gene.html"
else warn "no slices in $d - NAME=SPRTN_gene ./09_igv_slices.sh ... first; gene report skipped"; fi

ls -lh "$OUT"/*.html 2>/dev/null || true
log "copy to your laptop (a local folder, NOT OneDrive/iCloud - the pages embed patient reads):"
echo "  rsync -av pstancl@ssi-access.chem.pmf.hr:$OUT/*.html ~/igv_RJALS/"
