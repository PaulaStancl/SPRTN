#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 04 - Stage reference data ONCE into the shared references folder.
#
#   ./04_download_references.sh sarek          # iGenomes GATK.GRCh38, ~47 GB
#   ./04_download_references.sh oncoanalyser   # GRCh38_hmf + bwa-mem2 + WiGiTS, ~31 GB
#   ./04_download_references.sh vep            # VEP 115 GRCh38 cache (tumourevo), ~28 GB
#   ./04_download_references.sh all
#
# Safe to re-run: aws s3 sync skips finished files, nextflow uses -resume.
# Run inside tmux/screen.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
make_dirs

WHAT="${1:-}"
[[ "$WHAT" =~ ^(sarek|oncoanalyser|vep|all)$ ]] || die "usage: $0 sarek|oncoanalyser|vep|all"

# ---- sarek: iGenomes GATK.GRCh38 -------------------------------------------
# Full prefix is 63 GB. Only GermlineResource is skipped - it is the one large
# directory sarek's igenomes.config never references.
#
# BWAIndex (5.6 GB) and dragmap (6.8 GB) are downloaded even though we align with
# bwa-mem2: sarek's schema marks every genome path `exists: true` and validates
# ALL of them before the run, so missing indexes abort the pipeline with
# "--bwa ... does not exist" even when that aligner is never used.
download_sarek() {
    activate_tools
    local dest="$IGENOMES_BASE/Homo_sapiens/GATK/GRCh38"
    log "iGenomes GATK.GRCh38 -> $dest"
    mkdir -p "$dest"
    aws s3 sync --no-sign-request --region eu-west-1 \
        s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/ "$dest/" \
        --exclude 'Annotation/GermlineResource/*'
    for f in Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
             Sequence/BWAmem2Index \
             Annotation/ASCAT/G1000_loci_hg38.zip \
             Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz \
             Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed; do
        [[ -e "$dest/$f" ]] && ok "$f" || die "missing after sync: $dest/$f"
    done
    ok "sarek references: $(du -sh "$dest" | cut -f1)   use --igenomes_base $IGENOMES_BASE"
}

# ---- tumourevo: VEP cache ---------------------------------------------------
# sarek 3.10 ships VEP 116 but tumourevo ships VEP 115.2 - the cache must match
# tumourevo. Layout after sync: $VEP_CACHE/homo_sapiens/115_GRCh38/
download_vep() {
    activate_tools
    log "VEP cache ${VEP_CACHE_VERSION}_GRCh38 -> $VEP_CACHE"
    mkdir -p "$VEP_CACHE"
    aws s3 sync --no-sign-request --region eu-west-1 \
        "s3://annotation-cache/vep_cache/${VEP_CACHE_VERSION}_GRCh38/" "$VEP_CACHE/"
    [[ -d "$VEP_CACHE/homo_sapiens/${VEP_CACHE_VERSION}_GRCh38" ]] \
        && ok "VEP cache: $(du -sh "$VEP_CACHE/homo_sapiens/${VEP_CACHE_VERSION}_GRCh38" | cut -f1)" \
        || die "unexpected layout, expected $VEP_CACHE/homo_sapiens/${VEP_CACHE_VERSION}_GRCh38"
}

# ---- oncoanalyser: pipeline's own staging mode, then write a refdata config -
# shortest matching path (the top-level copy, not one inside an index folder)
find_one() { find -L "$1" "${@:2}" 2>/dev/null | awk '{print length, $0}' | sort -n | head -1 | cut -d' ' -f2-; }

download_oncoanalyser() {
    activate_env
    log "oncoanalyser --mode prepare_reference -> $HMF_REF_DIR"
    nf_run oncoanalyser_prepare_reference "$NXF_PROFILE" nf-core/oncoanalyser \
        -r "$ONCOANALYSER_REV" \
        --mode prepare_reference \
        --ref_data_types wgs,dna_alignment \
        --genome GRCh38_hmf \
        --outdir "$HMF_REF_DIR"

    local stage="$HMF_REF_DIR/reference_data/$ONCOANALYSER_REV"
    [[ -d "$stage" ]] || stage="$HMF_REF_DIR"

    # Extract any tarball that was published unextracted
    while IFS= read -r tgz; do
        local out="${tgz%.tar.gz}"
        [[ -d "$out" ]] && continue
        log "extracting $(basename "$tgz")"
        mkdir -p "$out" && tar -xzf "$tgz" -C "$out"
    done < <(find -L "$stage" -name '*.tar.gz' 2>/dev/null)

    local base=GRCh38_masked_exclusions_alts_hlas.fasta
    local fasta fai dict img bwamem2 gridss hmf
    fasta=$(find_one "$stage" -name "$base")
    fai=$(find_one "$stage" -name "$base.fai")
    dict=$(find_one "$stage" -name "$base.dict")
    img=$(find_one "$stage" -name "$base.img")
    bwamem2=$(dirname "$(find_one "$stage" -name '*.bwt.2bit.64')")
    gridss=$(dirname "$(find_one "$stage" -name '*.gridsscache')")
    hmf=$(find_one "$stage" -type d -path '*/dna/copy_number'); hmf="${hmf%/dna/copy_number}"

    for v in fasta fai dict img bwamem2 gridss hmf; do
        [[ -n "${!v}" && "${!v}" != "." && -e "${!v}" ]] && ok "$v = ${!v}" \
            || die "could not locate '$v' under $stage - check: find $stage -maxdepth 3"
    done

    cat > "$ONCO_REFDATA_CONFIG" <<EOF
// Written by 04_download_references.sh on $(date +%F) - locally staged
// oncoanalyser $ONCOANALYSER_REV reference data. Re-run 04 to regenerate.
params {
    genomes {
        GRCh38_hmf {
            fasta         = "$fasta"
            fai           = "$fai"
            dict          = "$dict"
            img           = "$img"
            bwamem2_index = "$bwamem2/"
            gridss_index  = "$gridss/"
        }
    }
    ref_data_hmf_data_path = "$hmf/"
}
EOF
    ok "wrote $ONCO_REFDATA_CONFIG"
    log "the prepare work dir can go:  rm -rf $NXF_WORK_BASE/oncoanalyser_prepare_reference"
}

case "$WHAT" in
    sarek)        download_sarek ;;
    oncoanalyser) download_oncoanalyser ;;
    vep)          download_vep ;;
    all)          download_sarek; download_vep; download_oncoanalyser ;;
esac

log "Next:  ./05_download_test_data.sh"
