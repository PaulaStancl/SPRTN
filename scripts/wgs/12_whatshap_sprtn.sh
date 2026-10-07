#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 12 - read-based phasing (WhatsHap) of the two germline SPRTN variants of the normal:
#      are c.350A>G (p.Tyr117Cys, chr1:231347825) and c.718_718+3del (4-bp deletion AGGT
#      after chr1:231351569) on the same or on opposite chromosome copies?
#
#   micromamba create -y -p /common/WORK/pstancl/envs/whatshap -c conda-forge -c bioconda whatshap   # once
#   ./12_whatshap_sprtn.sh
#
# The variants are ~3.7 kb apart and a read pair spans ~300 bp, so no read covers both.
# WhatsHap can still link them through the patient's heterozygous SNPs that lie BETWEEN
# them (each read pair links neighbouring SNPs; the chain bridges the gap). So the input is
# the two variants + the heterozygous PASS germline SNPs/indels between them (Strelka2,
# normal) - nothing else - and the report shows only the two variants.
# Reads: normal alone, and normal + tumour (germline variants are in both; more linking reads). Docs: https://whatshap.readthedocs.io/en/latest/guide.html
# Output: $RESULTS_BASE/phasing/SPRTN/ input.vcf.gz, {normal_only,normal_tumour}/ phased.vcf.gz, whatshap.log, result.txt
#   same PS (phase set) and opposite haplotypes (0|1 vs 1|0)  -> in trans
#   same PS and same haplotype (both 0|1 or both 1|0)          -> in cis
#   different PS, or a variant left unphased                   -> the reads cannot tell
# plus haplotype_imbalance.txt: tumour read fractions of the SNPs phased with the deletion (its
#   parental copy) vs Y117C - with SPRTN gained unequally in the tumour, the two copies differ
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ./00_config.sh
activate_tools                                       # bcftools, samtools
umask 077

WH="${WHATSHAP:-$ENV_ROOT/whatshap/bin/whatshap}"
[[ -x "$WH" ]] || die "no whatshap at $WH - micromamba create -y -p $ENV_ROOT/whatshap -c conda-forge -c bioconda whatshap"
FA="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
GERM="$RESULTS_BASE/sarek/$DATASET/variant_calling/strelka/$NORMAL_ID/$NORMAL_ID.strelka.variants.vcf.gz"
# reads: 09's whole-gene slices of sarek's CRAMs - every read in SPRTN, and fast (the full CRAMs
# would make WhatsHap read through much of chr1)
SL="$RESULTS_BASE/igv_slices/SPRTN_gene"
CRAM="$SL/${NORMAL_ID}_sarek.slice.bam"
TCRAM="$SL/${TUMOUR_ID}_sarek.slice.bam"
[[ -f "$GERM" ]] || die "missing $GERM"
for f in "$CRAM" "$TCRAM"; do
    [[ -f "$f" ]] || die "missing $f - run NAME=SPRTN_gene PAD=200 ./09_igv_slices.sh <gene + exon loci> first (see 10's header)"
done
OUT="$RESULTS_BASE/phasing/SPRTN"; mkdir -p "$OUT"

Y117C=231347825                       # A>G
DEL_ANCHOR=231351569                  # deletion of the 4 bases after it (AGGT), VCF-style
REGION="chr1:${Y117C}-$(( DEL_ANCHOR + 5 ))"
SAMPLE=$(bcftools query -l "$GERM" | head -1)

# ---- input VCF: heterozygous PASS germline variants between (and including) the two -----------
bcftools view -f PASS -g het -r "$REGION" "$GERM" -Oz -o "$OUT/between.vcf.gz"
bcftools index -t -f "$OUT/between.vcf.gz"
# the deletion, added from the reference if Strelka2 did not call it
ndel=$(bcftools view -H -r "chr1:$(( DEL_ANCHOR - 4 ))-$(( DEL_ANCHOR + 4 ))" "$OUT/between.vcf.gz" | awk '$5 != "." && length($4) > length($5) { c++ } END { print c + 0 }')
if (( ndel )); then
    ok "c.718_718+3del: called by Strelka2"
    cp "$OUT/between.vcf.gz" "$OUT/input.vcf.gz"
else
    ref=$(samtools faidx "$FA" "chr1:$DEL_ANCHOR-$(( DEL_ANCHOR + 4 ))" | grep -v '^>' | tr -d '\n' | tr a-z A-Z)
    [[ "${ref:1}" == "AGGT" ]] || warn "reference after $DEL_ANCHOR is ${ref:1}, not AGGT - check the position"
    warn "c.718_718+3del not in Strelka2's germline VCF - added as chr1:$DEL_ANCHOR $ref>${ref:0:1} (0/1) from the read counts"
    { bcftools view -h "$OUT/between.vcf.gz" | sed '$d'
      echo '##INFO=<ID=ADDED,Number=0,Type=Flag,Description="added by 12_whatshap_sprtn.sh: c.718_718+3del seen in reads, not called by Strelka2">'
      bcftools view -h "$OUT/between.vcf.gz" | tail -1
      printf 'chr1\t%s\t.\t%s\t%s\t.\tPASS\tADDED\tGT\t0/1\n' "$DEL_ANCHOR" "$ref" "${ref:0:1}"
      bcftools view -H "$OUT/between.vcf.gz" | awk -F'\t' -v OFS='\t' '{ split($9, f, ":"); split($10, v, ":"); for (i in f) if (f[i] == "GT") g = v[i]; $8 = "."; $9 = "GT"; $10 = g; print }'
    } | bcftools sort -Oz -o "$OUT/input.vcf.gz"
fi
bcftools index -t -f "$OUT/input.vcf.gz"
n=$(bcftools view -H "$OUT/input.vcf.gz" | wc -l)
log "heterozygous sites used (the two variants + $(( n - 2 )) between them): $n"
bcftools query -f '  %POS\t%REF>%ALT\n' "$OUT/input.vcf.gz"
(( n > 2 )) || warn "no heterozygous SNPs between the two variants - nothing for the reads to chain through"

# ---- phase ---------------------------------------------------------------------------------------
help=$("$WH" phase --help 2>&1 || true); IND=(); [[ "$help" == *"--indels"* ]] && IND=(--indels)

# phase_run <name> <label> <alignment files...>: phase, then report the two variants
phase_run() {
    local name="$1" label="$2" od="$OUT/$1"; shift 2; mkdir -p "$od"
    "$WH" phase --reference "$FA" ${IND[@]+"${IND[@]}"} --ignore-read-groups -o "$od/phased.vcf" \
        "$OUT/input.vcf.gz" "$@" > "$od/whatshap.log" 2>&1 || { tail -20 "$od/whatshap.log"; die "whatshap phase failed"; }
    bcftools view "$od/phased.vcf" -Oz -o "$od/phased.vcf.gz"; bcftools index -t -f "$od/phased.vcf.gz"

    # ---- result: the two variants only -------------------------------------------------------------
    {
    echo "# SPRTN phasing (WhatsHap, reads: $label, sarek alignments, SPRTN slices) - $(date '+%F %T')"
    echo "# heterozygous sites in the chain: $n (the two variants + $(( n - 2 )) between them)"
    printf '%-16s %-12s %-14s %-6s %s\n' variant position allele GT PS
    bcftools query -f '%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' \
        -r "chr1:$Y117C-$Y117C,chr1:$(( DEL_ANCHOR - 4 ))-$(( DEL_ANCHOR + 4 ))" "$od/phased.vcf.gz" \
      | awk -F'\t' -v y="$Y117C" '$1 == y || length($2) > length($3)' > "$od/two.tsv"
    awk -F'\t' -v y="$Y117C" '{ printf "%-16s %-12s %-14s %-6s %s\n", ($1 == y ? "c.350A>G" : "c.718_718+3del"), $1, $2 ">" $3, $4, $5 }' "$od/two.tsv"
    awk -F'\t' '
        { gt[NR] = $4; ps[NR] = $5 }
        END {
            if (NR < 2)                                         v = "one variant missing from the output - check whatshap.log"
            else if (gt[1] !~ /\|/ || gt[2] !~ /\|/)           v = "NOT PHASED - the reads do not link the two variants"
            else if (ps[1] == "." || ps[1] != ps[2])            v = "DIFFERENT PHASE BLOCKS - the reads do not link the two variants"
            else if (gt[1] == gt[2])                            v = "SAME HAPLOTYPE (in cis)"
            else                                                v = "OPPOSITE HAPLOTYPES (in trans)"
            print "# result: " v }' "$od/two.tsv"
    } | tee "$od/result.txt"
    log "saved: $od/result.txt"
}

# germline variants are in both samples, on the same chromosome copies - the tumour reads add
# linking evidence (about twice the reads); normal-only is shown as a check
phase_run normal_only "normal $NORMAL_ID" "$CRAM"
phase_run normal_tumour "normal $NORMAL_ID + tumour $TUMOUR_ID" "$CRAM" "$TCRAM"
log "results: $OUT/normal_only/result.txt and $OUT/normal_tumour/result.txt"

# ---- haplotype-level allelic imbalance: the deletion's copy in the tumour ------------------------
# WhatsHap put the deletion and some heterozygous SNPs into one phase block = one parental copy.
# For each SNP in that block: the fraction of reads carrying the allele that sits on the DELETION's
# copy (ALT if its GT is oriented like the deletion's, REF otherwise). Normal: ~0.5. Tumour: the
# deletion's copy is the amplified one (~0.68 expected from PURPLE 4.7 vs 1.9 copies, purity 0.6)
# or the other one (~0.32). Many SNPs give a firmer answer than the deletion's own VAF (which is
# under-counted). Y117C's tumour VAF on top: ~0.65 = on the amplified copy.
# If the deletion's copy is the LESS-amplified one while Y117C is on the amplified one -> in trans.
P="$OUT/normal_tumour/phased.vcf.gz"
H="$OUT/haplotype_imbalance.txt"
dps=$(bcftools query -f '%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' -r "chr1:$(( DEL_ANCHOR - 4 ))-$(( DEL_ANCHOR + 4 ))" "$P" \
      | awk -F'\t' 'length($2) > length($3) && !d { print $4 "\t" $5; d = 1 }')
del_gt="${dps%%$'\t'*}"; del_ps="${dps##*$'\t'}"
{
echo "# SPRTN haplotype-level allelic imbalance (phase block of c.718_718+3del, WhatsHap normal + tumour) - $(date '+%F %T')"
if [[ -z "$dps" || "$del_ps" == "." || "$del_gt" != *"|"* ]]; then
    echo "# the deletion is not phased - nothing to do"
else
    echo "# deletion: GT $del_gt in phase block $del_ps; SNPs of that block:"
    # SNPs in the block: position, REF, ALT, GT
    bcftools query -i 'TYPE="snp"' -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%PS]\n' "$P" \
      | awk -F'\t' -v ps="$del_ps" '$6 == ps && $5 ~ /\|/' > "$OUT/block_snps.tsv"
    cut -f1,2 "$OUT/block_snps.tsv" > "$OUT/block_snps.targets"
    if [[ ! -s "$OUT/block_snps.tsv" ]]; then
        echo "# no SNPs in the deletion's phase block"
    else
        # REF / ALT read counts at the block's SNPs, one sample at a time (MAPQ, BQ >= 20)
        ad_counts() {   # <alignment> <out>: POS  ref_reads  alt_reads (for the SNP's own ALT)
            bcftools mpileup -f "$FA" -T "$OUT/block_snps.targets" -a AD -q 20 -Q 20 -I --max-depth 100000 "$1" 2>/dev/null \
              | bcftools query -f '%POS\t%ALT\t[%AD]\n' \
              | awk -F'\t' 'FNR == NR { alt[$2] = $4; next }
                            ($1 in alt) { n = split($2, A, ","); k = 0; for (i = 1; i <= n; i++) if (A[i] == alt[$1]) k = i
                                          split($3, D, ","); print $1 "\t" D[1] "\t" (k ? D[k + 1] : 0) }' "$OUT/block_snps.tsv" - > "$2"
        }
        ad_counts "$CRAM" "$OUT/block_snps.normal"
        ad_counts "$TCRAM" "$OUT/block_snps.tumour"
        y117c_t=$(samtools mpileup -Q 20 -q 20 -f "$FA" -r chr1:231347825-231347825 "$TCRAM" 2>/dev/null \
                  | awk '{ b = toupper($5); r = gsub(/[.,]/, "", b); g = gsub(/G/, "", b); printf "%.2f", (r + g ? g / (r + g) : 0) }')
        awk -F'\t' -v dgt="$del_gt" -v y117c_t="$y117c_t" '
            function med(a, n,   i, j, t) { for (i = 2; i <= n; i++) { t = a[i]; for (j = i - 1; j >= 1 && a[j] > t; j--) a[j + 1] = a[j]; a[j + 1] = t }
                                          return n % 2 ? a[(n + 1) / 2] : (a[n / 2] + a[n / 2 + 1]) / 2 }
            FILENAME ~ /block_snps.tsv$/    { ref[$2] = $3; alt[$2] = $4; gt[$2] = $5; order[++np] = $2; next }
            FILENAME ~ /block_snps.normal$/ { nr[$1] = $2; na[$1] = $3; next }
            FILENAME ~ /block_snps.tumour$/ { tr[$1] = $2; ta[$1] = $3; next }
            END {
                for (q = 1; q <= np; q++) {
                    p = order[q]
                    if (nr[p] + na[p] < 10 || tr[p] + ta[p] < 10) continue
                    same = (substr(gt[p], 1, 1) == substr(dgt, 1, 1))           # ALT on the deletion copy?
                    cn = same ? na[p] : nr[p]; ct = same ? ta[p] : tr[p]
                    hn = cn / (nr[p] + na[p]); ht = ct / (tr[p] + ta[p])
                    printf "  %s  %s>%s  GT %s  deletion-copy reads: normal %d/%d (%.2f)  tumour %d/%d (%.2f)\n",
                           p, ref[p], alt[p], gt[p], cn, nr[p] + na[p], hn, ct, tr[p] + ta[p], ht
                    m++; HN[m] = hn; HT[m] = ht
                }
                if (!m) { print "# no SNP with enough reads (>= 10 in each sample)"; exit }
                mn = med(HN, m); mt = med(HT, m)
                printf "# %d SNPs on the deletion copy: median fraction normal %.2f, tumour %.2f\n", m, mn, mt
                printf "# Y117C tumour VAF: %s   (expected ~0.68 on the amplified copy, ~0.32 on the other)\n", y117c_t
                if (mt < 0.45 && y117c_t + 0 > 0.55)      v = "deletion copy LESS amplified, Y117C on the MORE amplified copy -> supports IN TRANS"
                else if (mt > 0.55 && y117c_t + 0 > 0.55) v = "deletion copy and Y117C both on the MORE amplified copy -> suggests IN CIS"
                else                                      v = "no clear imbalance - inconclusive"
                print "# result: " v
            }' "$OUT/block_snps.tsv" "$OUT/block_snps.normal" "$OUT/block_snps.tumour"
    fi
fi
} | tee "$H"
log "saved: $H"
