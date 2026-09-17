#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 08 - Installation check with each pipeline's own small test profile.
#      Minutes to ~1 h each; also fills the container cache.
#
#   ./08_test_pipelines.sh sarek|oncoanalyser|tumourevo|all
#
#   sarek        -profile test          tiny FASTQ, strelka
#   oncoanalyser -profile test          simulated tumour/normal DNA FASTQ on the
#                                       staged GRCh38_hmf refs (needs 04 oncoanalyser)
#   tumourevo    -profile test          simulated chr17 VCF + ASCAT
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./00_config.sh
activate_env
make_dirs

WHAT="${1:-}"
[[ "$WHAT" =~ ^(sarek|oncoanalyser|tumourevo|all)$ ]] || die "usage: $0 sarek|oncoanalyser|tumourevo|all"
[[ "$CONTAINER_ENGINE" != none ]] || die "no container engine"

OUT="$RESULTS_BASE/_install_tests"

test_sarek() {
    nf_run test_sarek "test,$NXF_PROFILE" nf-core/sarek -r "$SAREK_REV" \
        --outdir "$OUT/sarek"
}

test_oncoanalyser() {
    [[ -f "$ONCO_REFDATA_CONFIG" ]] || die "missing $ONCO_REFDATA_CONFIG - run ./04_download_references.sh oncoanalyser"
    # DNA-only minimal sheet: the default test sheet includes RNA and would
    # pull the 27 GB STAR index.
    nf_run test_oncoanalyser "test,$NXF_PROFILE" nf-core/oncoanalyser -r "$ONCOANALYSER_REV" \
        -c "$ONCO_REFDATA_CONFIG" \
        --input https://raw.githubusercontent.com/nf-core/test-datasets/oncoanalyser/samplesheet/fastq_eval.subject_a.wgts.tndna.minimal.csv \
        --outdir "$OUT/oncoanalyser"
}

test_tumourevo() {
    # tumourevo dev still uses pre-strict Nextflow syntax; 26.04 is strict by default.
    #
    # --tools: the test profile lists sparsesignatures, but the pipeline's own
    # nf-test (tests/default.nf.test) overrides tools to "tinc,mobster,pyclone-vi",
    # i.e. upstream does not test signature extraction on this data either. The
    # chr17 test mutations are too few for SparseSignatures cross-validation: every
    # grid MSE comes back NA, min K is NA, and the module dies in `if (K < 2)`.
    NXF_SYNTAX_PARSER="${NXF_SYNTAX_PARSER:-v1}" \
    nf_run test_tumourevo "test,$NXF_PROFILE" nf-core/tumourevo -r "$TUMOUREVO_REV" \
        --tools tinc,mobster,pyclone-vi \
        --outdir "$OUT/tumourevo"
}

case "$WHAT" in
    sarek)        test_sarek ;;
    oncoanalyser) test_oncoanalyser ;;
    tumourevo)    test_tumourevo ;;
    all)          test_sarek; test_tumourevo; test_oncoanalyser ;;
esac

log "Test(s) passed. Results: $OUT"
log "Next:  ./10_run_sarek.sh  and  ./11_run_oncoanalyser.sh  (can run one after the other)"
