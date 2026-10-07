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
