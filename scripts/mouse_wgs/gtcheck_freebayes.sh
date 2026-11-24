#!/usr/bin/env bash
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24
#
# gtcheck_freebayes.sh - tumour / KI normal sample identity with bcftools gtcheck, from the FreeBayes VCFs
# of setting A (no BAMs needed).
#
# Steps
#   1. per FreeBayes file: biallelic SNVs, autosomes 1-19, QUAL >= MIN_QUAL
#   2. low-quality genotypes -> missing: FORMAT/DP < MIN_DP or FORMAT/GQ < MIN_GQ
#        bcftools >= 1.10: +setGT masks only that sample's genotype
#        older bcftools:   the site is dropped unless both samples pass (MIN(FMT/DP), MIN(FMT/GQ))
#   3. merge the three filtered files (-m none: no joining of different alleles)
#   4. gtcheck twice:
#        all_sites    - every site with a genotype in both samples of a pair (labelled pairs, which share a
#                       file, are compared at more sites than pairs from different files)
#        common_sites - only sites genotyped in ALL six samples (-g ^miss): every pair on the same sites.
#                       Use this one for the comparison.
#      Missing genotypes stay missing (never --missing-to-ref: the files were called per pair, so a site
#      absent from a file is unknown, not hom-ref).
#   5. table: discordance per 1000 compared sites, short names, labelled pairs marked
#
# gtcheck options by version: >= 1.19 "-u GT,GT -E 0"; 1.11-1.18 "-u GT,GT -e 0"; older "-G 1".
# Output lines: DCv2 / DC (query, genotyped, discordance, -logP(HWE), n_sites ...) or CN (old:
# discordance, n_sites, avg min depth, sample i, sample j).
# Note: this uses the same FreeBayes calls as the IBS0 analysis, so it is a second method, not
# independent evidence (that would need BAM-level pileups or the sample records).
#
# Run on ssi (r-mutation env; bcftools >= 1.10 recommended):
#   bash scripts/mouse_wgs/gtcheck_freebayes.sh
# Output: results/mouse_wgs/gtcheck/
set -euo pipefail

# ---- EDIT THIS BLOCK ---------------------------------------------------------
REPO=/common/WORK/pstancl/projects/MariaBoskovic/SPRTN
DIR_A=$REPO/data/mouse/wgs/tumor_vs_normal/snv_indel
OUT=$REPO/results/mouse_wgs/gtcheck
MIN_QUAL=30
MIN_DP=30
MIN_GQ=30
# labelled pairs (tumour normal), used only to mark rows in the table
LABELLED="HCC_1:KI1 HCC_2:KI2 HCC_3:KI3"

# ---- setup -------------------------------------------------------------------
command -v bcftools >/dev/null || PATH=/common/WORK/pstancl/envs/r-mutation/bin:$PATH
mkdir -p "$OUT/filtered"
cd "$OUT"
ver=$(bcftools --version | head -1 | awk '{print $2}')
major=${ver%%.*}; rest=${ver#*.}; minor=${rest%%[!0-9]*}
echo "bcftools $ver"
newer() { [ "$major" -gt "$1" ] || { [ "$major" -eq "$1" ] && [ "$minor" -ge "$2" ]; }; }
if newer 1 10; then MASK=setgt; else MASK=site; echo "NOTE: bcftools < 1.10 - low-quality genotypes drop the whole site (upgrade for per-sample masking)"; fi
if newer 1 19; then GT_OPT="-u GT,GT -E 0"; elif newer 1 11; then GT_OPT="-u GT,GT -e 0"; else GT_OPT="-G 1"; fi
echo "gtcheck options: $GT_OPT; genotype masking: $MASK"

# ---- 1-2. filter each FreeBayes file -------------------------------------------
inputs=()
for f in "$DIR_A"/*/*freebayes*.bcf "$DIR_A"/*freebayes*.bcf; do
  [ -e "$f" ] || continue
  case "$f" in *union*|*merged*) continue ;; esac
  b=$(basename "$f" .bcf)
  chr=""; if bcftools view -h "$f" | grep '^##contig=<ID=chr1,' >/dev/null; then chr=chr; fi
  autos=$(seq -s, -f "${chr}%g" 1 19)
  out="filtered/$b.filtered.bcf"
  if [ "$MASK" = setgt ]; then
    bcftools view -v snps -m2 -M2 -t "$autos" -i "QUAL>=$MIN_QUAL" -Ou "$f" |
      bcftools +setGT -Ob -o "$out" -- -t q -n . -i "FMT/DP<$MIN_DP | FMT/GQ<$MIN_GQ"
  else
    bcftools view -v snps -m2 -M2 -t "$autos" -i "QUAL>=$MIN_QUAL && MIN(FMT/DP)>=$MIN_DP && MIN(FMT/GQ)>=$MIN_GQ" -Ob -o "$out" "$f"
  fi
  bcftools index -f "$out"
  echo "$b: $(bcftools view -H "$out" | wc -l) sites"
  inputs+=("$out")
done
[ "${#inputs[@]}" -eq 3 ] || { echo "expected 3 FreeBayes files, found ${#inputs[@]}"; exit 1; }

# ---- 3. merge ----------------------------------------------------------------
bcftools merge -m none -Ob -o freebayes_merged.bcf "${inputs[@]}"
bcftools index -f freebayes_merged.bcf
bcftools view -g ^miss -Ob -o freebayes_merged_common.bcf freebayes_merged.bcf
bcftools index -f freebayes_merged_common.bcf
echo "merged sites: $(bcftools view -H freebayes_merged.bcf | wc -l); genotyped in all 6 samples: $(bcftools view -H freebayes_merged_common.bcf | wc -l)"

# ---- 4. gtcheck ----------------------------------------------------------------
bcftools gtcheck $GT_OPT freebayes_merged.bcf        > gtcheck_all_sites.txt
bcftools gtcheck $GT_OPT freebayes_merged_common.bcf > gtcheck_common_sites.txt

# ---- 5. table ----------------------------------------------------------------
# columns: site_set s1 s2 discordance n_sites per_1000 labelled_pair
{ printf 'site_set\ts1\ts2\tdiscordance\tn_sites\tper_1000\tlabelled_pair\n'
for set in all_sites common_sites; do
  awk -F'\t' -v set="$set" -v lab="$LABELLED" '
    function short(s) { if (s ~ /KI_HCC_/) { sub(/.*KI_HCC_/, "HCC_", s) } else { sub(/.*_/, "", s) } return s }
    BEGIN { n = split(lab, L, " "); for (i = 1; i <= n; i++) { split(L[i], p, ":"); pair[p[1] ":" p[2]] = 1; pair[p[2] ":" p[1]] = 1 } }
    $1 == "DCv2" || $1 == "DC" { a = short($2); b = short($3); d = $4; n = $6 }
    $1 == "CN"                 { a = short($5); b = short($6); d = $2; n = $3 }
    ($1 == "DCv2" || $1 == "DC" || $1 == "CN") && n > 0 {
      printf "%s\t%s\t%s\t%s\t%s\t%.1f\t%s\n", set, a, b, d, n, 1000 * d / n, ((a ":" b) in pair) ? "yes" : "no" }
  ' "gtcheck_$set.txt"
done | sort -t$'\t' -k1,1 -k6,6g; } > gtcheck_table.tsv

echo; echo "== common sites (same sites for every pair; lower per_1000 = more alike) =="
grep -E '^(site_set|common_sites)' gtcheck_table.tsv | column -t
echo; echo "Tables: $OUT/gtcheck_table.tsv, raw: gtcheck_all_sites.txt, gtcheck_common_sites.txt"
rm -rf filtered
