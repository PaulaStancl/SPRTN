#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 11 - SPRTN second-allele checks, all in one report: is Y117C (heterozygous in the normal)
#      accompanied by a truncating SPRTN variant on the other allele (compound heterozygous,
#      as in Lessel et al. 2014), and are the samples / lanes consistent?
#
#   ./11_sprtn_checks.sh
#
# Needs 08_region_variants.sh (same folder) and 09's whole-gene slices:
#   NAME=SPRTN_gene PAD=200 ./09_igv_slices.sh <gene + exon loci>   (see 10's header)
# Output: $RESULTS_BASE/region_variants/SPRTN_checks.txt  (+ one <name>.txt per 08 search)
#   1. every VCF record at the ClinVar truncation positions and in the whole gene (08)
#   2. reads with a deletion starting near each truncation, normal and tumour (pileup)
#   3. Y117C reads per lane (read group), normal and tumour - a mixed-up lane would differ
#   4. PURPLE copy number of SPRTN (LOH -> VAF shifts in the tumour can phase two variants)
# Truncations (ClinVar, GRCh38): c.723del p.Lys241fs 231352612 (the C-terminal truncation);
#   c.718_718+3del 231351571-4; c.1246_1247del p.Val416fs 231353136-7. Y117C: 231347825 A>G.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ./00_config.sh
activate_tools
umask 077

OUT_DIR="$RESULTS_BASE/region_variants"; mkdir -p "$OUT_DIR"
REPORT="$OUT_DIR/SPRTN_checks.txt"
S="$RESULTS_BASE/igv_slices/SPRTN_gene"
FA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
CNV_GENE="$RESULTS_BASE/oncoanalyser/$DATASET/$DATASET/purple/$TUMOUR_ID.purple.cnv.gene.tsv"
for s in "$NORMAL_ID" "$TUMOUR_ID"; do
    [[ -f "$S/${s}_sarek.slice.bam" ]] || die "missing $S/${s}_sarek.slice.bam - run 09 for SPRTN_gene first (see 10's header)"
done

{
echo "# SPRTN second-allele checks - $DATASET - $(date '+%F %T')"
echo

echo "## 1. VCF records (all VCFs; 08_region_variants.sh) ------------------------------"
for x in "chr1:231352600-231352620 SPRTN_K241fs" "chr1:231351560-231351580 SPRTN_c718del" \
         "chr1:231353125-231353145 SPRTN_V416fs" "chr1:231337104-231375416 SPRTN_gene"; do
    set -- $x
    echo; echo "### $2  ($1)"
    ./08_region_variants.sh "$1" "$2" | grep -vE '^\[|^# ' | grep -E 'PASS: [1-9]|all: [1-9]|^  ' || echo "  (no records in any VCF)"
done

echo; echo "## 2. reads with a deletion starting near each truncation (sarek alignments) -------"
echo "#    ~half of the normal's reads = heterozygous carrier; 1-2 reads = noise"
for x in "chr1:231352600-231352620 p.Lys241fs" "chr1:231351560-231351580 c.718_718+3del" "chr1:231353125-231353145 p.Val416fs"; do
    set -- $x
    for s in "$NORMAL_ID" "$TUMOUR_ID"; do
        samtools mpileup -Q 20 -q 20 -f "$FA" -r "$1" "$S/${s}_sarek.slice.bam" 2>/dev/null \
          | awk -v v="$2" -v s="$s" '{
                b = $5; n = 0; delete c
                while (match(b, /-[0-9]+/)) {                     # each -<len><bases> in the pileup
                    len = substr(b, RSTART + 1, RLENGTH - 1) + 0
                    seq = toupper(substr(b, RSTART + RLENGTH, len)); c[len "bp " seq]++; n++
                    b = substr(b, RSTART + RLENGTH + len)
                }
                if (n > 0) { out = ""; for (k in c) out = out " " k "=" c[k]
                    printf "  %-16s %-9s %s:%s  depth %d  deletion reads %d  (deleted:%s)\n", v, s, $1, $2, $4, n, out; hit = 1 } }
              END { if (!hit) printf "  %-16s %-9s no deletion reads\n", v, s }'
    done
done

echo; echo "## 3. Y117C (chr1:231347825 A>G) per lane / read group ---------------------------"
echo "#    every lane of a sample should give a similar VAF; one deviating lane = possible mix-up"
for s in "$NORMAL_ID" "$TUMOUR_ID"; do
    f="$S/${s}_sarek.slice.bam"
    for rg in $(samtools view -H "$f" | awk '$1 == "@RG" { for (i = 2; i <= NF; i++) if ($i ~ /^ID:/) print substr($i, 4) }'); do
        # region query on the indexed slice first; piped input has no index, so no -r in mpileup
        samtools view -b -r "$rg" "$f" chr1:231347825-231347825 | samtools mpileup -Q 20 -q 20 -f "$FA" - 2>/dev/null \
          | awk -v s="$s" -v rg="$rg" '$2 == 231347825 {b = toupper($5); r = gsub(/[.,]/, "", b); g = gsub(/G/, "", b)
                                      printf "  %-9s %-34s A=%-4d G=%-4d VAF=%.2f\n", s, rg, r, g, (r + g ? g / (r + g) : 0) }'
    done
done

echo; echo "## 3b. both variants, all lanes: normal vs tumour (VAF shifts + LOH -> phase) --------"
for s in "$NORMAL_ID" "$TUMOUR_ID"; do
    samtools mpileup -Q 20 -q 20 -f "$FA" -r chr1:231347825-231347825 "$S/${s}_sarek.slice.bam" 2>/dev/null \
      | awk -v s="$s" '{b = toupper($5); r = gsub(/[.,]/, "", b); g = gsub(/G/, "", b)
                        printf "  Y117C          %-9s A=%-4d G=%-4d VAF=%.2f\n", s, r, g, (r + g ? g / (r + g) : 0) }'
    samtools mpileup -Q 20 -q 20 -f "$FA" -r chr1:231351569-231351569 "$S/${s}_sarek.slice.bam" 2>/dev/null \
      | awk -v s="$s" '{n = gsub(/-[0-9]+[ACGTNacgtn]+/, "", $5); printf "  c.718_718+3del %-9s deletion %d of %d reads, VAF=%.2f\n", s, n, $4, ($4 ? n / $4 : 0) }'
done

echo; echo "## 4. PURPLE copy number of SPRTN (tumour) -----------------------------------------"
echo "#    minor vs major allele copies (= copy number - minor): LOH (minor ~0) or any imbalance"
echo "#    (e.g. 4.7 vs 1.9 copies) phases two heterozygous variants: one on the more-amplified copy"
echo "#    rises in tumour VAF, one on the other copy stays / falls - opposite moves = in trans"
if [[ -f "$CNV_GENE" ]]; then
    paste <(head -1 "$CNV_GENE" | tr '\t' '\n') <(awk -F'\t' '$0 ~ /\tSPRTN\t/' "$CNV_GENE" | head -1 | tr '\t' '\n') \
      | grep -iE 'gene|copyNumber|minorAllele|somaticRegions|germline' | sed 's/^/  /'
else echo "  (no $CNV_GENE)"; fi
} 2>&1 | tee "$REPORT"
log "saved: $REPORT"
