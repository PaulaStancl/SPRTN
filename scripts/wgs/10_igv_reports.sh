#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 10 - igv-reports: interactive HTML pages (igv.js) of the reads in SPRTN, tumour and
#      normal from both pipelines. Runs on the server; no IGV, no display needed.
#
#   micromamba create -y -p /common/WORK/pstancl/envs/igvreports -c conda-forge -c bioconda igv-reports   # once
#   NAME=SPRTN_gene PAD=200 ./09_igv_slices.sh chr1:231337104-231375416 \
#       chr1:231338243-231338654 chr1:231339719-231339918 chr1:231347747-231347975 \
#       chr1:231351254-231351621 chr1:231352560-231355073          # slices of the whole gene
#   ./10_igv_reports.sh
#
# Output: $RESULTS_BASE/igv_reports/
#   SPRTN_patient_variants.html   every PASS variant the patient carries in SPRTN: germline
#                                 (Strelka2, normal; GT / AD / DP in the table) and somatic
#                                 (Mutect2, Strelka2, SAGE/PURPLE), each a row
#   SPRTN_clinvar_pathogenic.html ClinVar pathogenic / likely pathogenic SPRTN variants
#                                 (sprtn_clinvar_variants.tsv): Y117C, c.723del (p.Lys241fs,
#                                 the C-terminal truncation), c.718_718+3del, c.1246_1247del -
#                                 the reads at each position, carried or not
#   SPRTN_clinvar_vus.html        the same for ClinVar's variants of uncertain significance
#   SPRTN_compound_het.html       the patient's two variants side by side: Y117C + c.718_718+3del
#   SPRTN_gene.html               the whole gene and each exon
# Tracks: 09's SPRTN_gene slices (RJALS_Tm / RJALS_N x sarek / oncoanalyser) + exon track.
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
GENE="chr1:231337104-231375416"                      # SPRTN, Ensembl ENSG00000010072
SL="$RESULTS_BASE/igv_slices/SPRTN_gene"
OUT="$RESULTS_BASE/igv_reports"; mkdir -p "$OUT"
TMP="$OUT/tmp"; mkdir -p "$TMP"
CLINVAR="$WGS_SCRIPTS/sprtn_clinvar_variants.tsv"
EXONS="$WGS_SCRIPTS/igv_SPRTN_exons.bed"
VC="$RESULTS_BASE/sarek/$DATASET/variant_calling"
PAIR="${TUMOUR_ID}_vs_${NORMAL_ID}"

# ---- tracks: the four whole-gene slices (tumour first) + exons ------------------------------
TRACKS=()
for s in "$TUMOUR_ID" "$NORMAL_ID"; do for pl in sarek oncoanalyser; do
    f="$SL/${s}_${pl}.slice.bam"
    [[ -f "$f" ]] || die "missing $f - run NAME=SPRTN_gene PAD=200 ./09_igv_slices.sh <gene + exon loci> first (see header)"
    TRACKS+=("$f")
done; done
[[ -f "$EXONS" ]] && { grep -v '^track' "$EXONS" > "$TMP/SPRTN_exons.bed"; TRACKS+=("$TMP/SPRTN_exons.bed"); }
report() {     # <sites> <output name> <title> [extra create_report options...]
    local sites="$1" name="$2" title="$3"; shift 3
    "$IGVR" "$sites" --fasta "$FASTA" --tracks "${TRACKS[@]}" --flanking 100 --title "$title" "$@" \
        --output "$OUT/$name.html"
    ok "$OUT/$name.html"
}

# ---- 1. every PASS variant the patient carries in SPRTN ----------------------------------------
# germline: Strelka2 on the normal; somatic: Mutect2, Strelka2, SAGE/PURPLE (none expected in
# SPRTN so far, but shown if there are any). One VCF per origin, as igv-reports takes one file.
GERM="$VC/strelka/$NORMAL_ID/$NORMAL_ID.strelka.variants.vcf.gz"
[[ -f "$GERM" ]] || die "missing $GERM"
bcftools view -f PASS -r "$GENE" "$GERM" -Oz -o "$TMP/germline.vcf.gz"; bcftools index -t -f "$TMP/germline.vcf.gz"
n=$(bcftools view -H "$TMP/germline.vcf.gz" | wc -l)
log "germline PASS variants in SPRTN (Strelka2, normal): $n"
(( n )) && report "$TMP/germline.vcf.gz" SPRTN_patient_variants \
    "SPRTN - germline PASS variants of RJALS (Strelka2, normal $NORMAL_ID)" --sample-columns GT AD DP

som=()
for f in "$VC/mutect2/$PAIR/$PAIR.mutect2.filtered.vcf.gz" \
         "$VC/strelka/$PAIR/$PAIR.strelka.somatic_snvs.vcf.gz" "$VC/strelka/$PAIR/$PAIR.strelka.somatic_indels.vcf.gz" \
         "$RESULTS_BASE/oncoanalyser/$DATASET/$DATASET/purple/$TUMOUR_ID.purple.somatic.vcf.gz"; do
    [[ -f "$f" ]] || continue
    k=$(bcftools view -H -f PASS -r "$GENE" "$f" | wc -l)
    log "somatic PASS in SPRTN, $(basename "$f"): $k"
    (( k )) && som+=("$f")
done
if (( ${#som[@]} )); then
    : > "$TMP/somatic.bed"
    for f in "${som[@]}"; do
        bcftools query -i 'FILTER="PASS"' -r "$GENE" -f '%CHROM\t%POS0\t%END\tsomatic %REF>%ALT ('"$(basename "$f" .vcf.gz)"')\n' "$f" >> "$TMP/somatic.bed"
    done
    report "$TMP/somatic.bed" SPRTN_somatic_variants "SPRTN - somatic PASS variants (tumour $TUMOUR_ID)"
fi

# ---- 2. ClinVar positions: pathogenic / likely pathogenic, and uncertain -----------------------
[[ -f "$CLINVAR" ]] || die "missing $CLINVAR"
awk -F'\t' '!/^#/ && $1 != "chrom" && $4 ~ /athogenic/ {printf "%s\t%d\t%d\t%s: %s %s\n", $1, $2 - 1, $3, $4, $5, $6}' "$CLINVAR" > "$TMP/clinvar_plp.bed"
awk -F'\t' '!/^#/ && $1 != "chrom" && $4 !~ /athogenic/ {printf "%s\t%d\t%d\t%s: %s %s\n", $1, $2 - 1, $3, $4, $5, $6}' "$CLINVAR" > "$TMP/clinvar_vus.bed"
log "ClinVar sites: $(wc -l < "$TMP/clinvar_plp.bed") pathogenic/likely pathogenic, $(wc -l < "$TMP/clinvar_vus.bed") uncertain"
report "$TMP/clinvar_plp.bed" SPRTN_clinvar_pathogenic \
    "SPRTN - ClinVar pathogenic / likely pathogenic positions (Y117C, p.Lys241fs = C-terminal truncation, ...) - RJALS reads"
[[ -s "$TMP/clinvar_vus.bed" ]] && report "$TMP/clinvar_vus.bed" SPRTN_clinvar_vus \
    "SPRTN - ClinVar variants of uncertain significance - RJALS reads"

# ---- 2b. the patient's two SPRTN variants side by side (compound heterozygous) ----------------
# Y117C (missense, exon 3) and c.718_718+3del (4-bp deletion AGGT at the exon 4 splice donor;
# pileup places it after 231351569, ClinVar/HGVS at 231351571-574 - the same deletion).
printf 'chr1\t231347824\t231347825\tc.350A>G p.Tyr117Cys (missense, exon 3)\n' >  "$TMP/compound_het.bed"
printf 'chr1\t231351568\t231351574\tc.718_718+3del (4-bp deletion, exon 4 splice donor)\n' >> "$TMP/compound_het.bed"
report "$TMP/compound_het.bed" SPRTN_compound_het \
    "SPRTN - the two germline variants of RJALS: c.350A>G (p.Tyr117Cys) and c.718_718+3del"

# ---- 3. the whole gene and each exon ------------------------------------------------------------
awk -F'\t' 'NR > 1 { printf "%s\t%d\t%d\tsite%s %s\n", $3, $4 - 1, $5, $1, $2 }' "$SL/sites.tsv" > "$TMP/gene.bed"
report "$TMP/gene.bed" SPRTN_gene "SPRTN gene and exons - RJALS"

ls -lh "$OUT"/*.html
log "copy to your laptop (a local folder, NOT OneDrive/iCloud - the pages embed patient reads):"
echo "  rsync -av \"pstancl@ssi-access.chem.pmf.hr:$OUT/*.html\" ~/igv_RJALS/"
