#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24 | Rewritten: 2026-11-24 | Revised: 2026-11-25
#
# mouse_wgs_somalier.R - tumour / KI-normal sample identity from the setting A FreeBayes VCFs, with
# several methods on ONE masked genotype file: discordance at normal-homozygous sites (primary), IBS0,
# somalier (Pedersen et al., Genome Medicine 2020) and bcftools gtcheck. No BAMs.
#
# Principle: samples of one animal share their inherited genotypes. [CHANGE 22] The primary metric is the
# discordance at sites where the NORMAL is homozygous (Conpair-style): tumour genotype != normal genotype.
# LOH only turns het into hom, so a tumour het against a hom normal (or the opposite homozygote) is
# mismatch evidence, and it stays observable when every animal is forced non-ref by the per-pair
# ascertainment (het vs 1/1). IBS0 (0/0 vs 1/1) is kept as a secondary method: it can detect a swap but,
# with this ascertainment, cannot confirm correct labels (IBS0 is ~0 for every pair when all labels are
# right). IBS2 and homozygous concordance are secondary, somalier relatedness is exploratory. [CHANGE 21] (replaces the old claim that relatedness "assumes population allele frequencies":
# somalier v0.3.5 relatedness = 2*(shared_hets - 2*ibs0) / het_ab, het_ab = hets of either sample at sites
# callable in both - KING-robust style, no allele frequencies; it is unstable with the few hets of
# inbred mice.)
#
# Inputs: the 3 FreeBayes BCFs of setting A (nf-core/sarek), one per labelled tumour-normal pair, 2
# samples each. The input folder is read-only; everything is written under OUT.
#
# Site selection (block 2), per FreeBayes file then combined:
#   biallelic SNV, autosomes 1-19, QUAL >= MIN_QUAL; AD must exist (and = RO,AO where RO/AO exist)
#   REF+ALT >= MIN_DP in ALL SIX samples, i.e. the site is in all 3 files (ascertainment fix, see CHANGE 7)
#   genotype bins from VAF = ALT/(REF+ALT): <= 0.05 hom-ref, 0.35-0.65 het, >= 0.95 hom-alt, else ambiguous
#   every KI normal confidently genotyped (none ambiguous), >= 1 KI normal non-ref
#   no sample deeper than MAX_DP_FACTOR x its median depth (collapsed duplications); optional MASK_BED;
#   optionally not a FreeBayes PASS record (provider's somatic prefilter; EXCLUDE_FB_PASS, default FALSE,
#   TRUE run as a sensitivity analysis)
#   tumours: same bins, ambiguous -> MASKED (./.), never hom-ref; tumours are never used to select sites
# Panels on the final site set: A all; B KI1-3 not all identical; C opposite homs (0/0 vs 1/1) in >= 2 KI
# normals. Thinning: at most one site per window (greedy minimum distance; sites closest to their samples'
# median depth first, seeded random tie-break), primary 20 kb, sensitivity 100 kb, 500 kb, 1 Mb. Every method
# runs on every panel x thinning. Linkage-aware checks on the primary combination: paired sign tests,
# leave-one-chromosome-out jackknife, 10 Mb block bootstrap, discordant-segment counts.
#
# Blocks:
#   0  EDIT THIS BLOCK
#   1  setup, tool checks, versions, help texts, logging helpers
#   2  germline site selection, masked 6-sample genotype file, panels, thinning
#   3  somalier extract (per panel x thinning, on the masked file)
#   4  somalier relate; bcftools gtcheck; custom metrics (normal-hom discordance, IBS0); paired tests
#   5  matching (normal-hom discordance primary), assignments, bootstrap, jackknife, segments, sensitivity
#      analyses (ascertainment files 1+2, MIN_DP), plots
#   6  outputs, RUN_LOG.md (changelog, assumptions, limitations, deviations), TEST_MODE checks
#
# Install somalier once (static binary, no env needed):
#   mkdir -p /common/WORK/pstancl/tools/somalier && curl -sSL -o /common/WORK/pstancl/tools/somalier/somalier https://github.com/brentp/somalier/releases/download/v0.3.5/somalier && chmod +x /common/WORK/pstancl/tools/somalier/somalier
#
# Run on ssi (r-mutation env: data.table, ggplot2, bcftools):
#   Rscript scripts/mouse_wgs/mouse_wgs_somalier.R           # real data
#   Rscript scripts/mouse_wgs/mouse_wgs_somalier.R --test    # synthetic data: scenario "swap" and "null"
# Output: results/mouse_wgs/mouse_wgs_somalier/ (TEST_MODE: results/mouse_wgs/mouse_wgs_somalier_TEST/<scenario>/)
#
# CHANGELOG (2026-11-24 rewrite; each item is tagged [CHANGE n] in the code)
#   1  TEST_MODE / --test: synthetic FreeBayes-like VCFs (known HCC_1/HCC_2 swap, ambiguous VAFs, sites
#      missing from one file) and a synthetic FASTA; checks at the end. Without bcftools/somalier (e.g. on
#      a laptop) the VCF text is parsed in R and the two tools are skipped with a DEVIATION line.
#   2  New parameters: genotype bins, THIN_PRIMARY / THIN_SENS (replace SITE_SPACING), MIN_PANEL_SITES,
#      MIN_JOINT_SITES, KNOCKIN_LOCUS / KNOCKIN_FLANK, HYPOTHESIS. MIN_SAMPLES_DP removed (CHANGE 7).
#   3  Logging helpers: DEVIATION lines, a filter-count table, run-generated folders cleared each run.
#   4  Versions: somalier prints "somalier version: x" to stderr on every call ("--version" is not an
#      option and returned "unknown program"); taken from there. bcftools and htslib versions both logged.
#      `somalier extract/relate --help` and `bcftools gtcheck` usage saved to tool_help/.
#   5  FreeBayes read: GT, AD, RO, AO (was AD only); AD required; AD == RO,AO checked and logged.
#   6  Assertions: one tumour + one KI normal per file, the six expected short names, VCF contigs in the
#      FASTA .fai (chr / no-chr).
#   7  Ascertainment fix: REF+ALT >= MIN_DP in ALL SIX samples (was >= 4 of 6). Each file only holds sites
#      variant in its own pair, so with 4 of 6 a labelled pair (same file) was compared at more sites than
#      a cross pair (two files). Number of sites removed by this is logged.
#   8  Genotype bins: KI normals with VAF 0.65-0.95 are now ambiguous too (the old check only caught
#      0.05-0.35); tumour genotypes binned the same way, ambiguous/low-depth ones masked, never hom-ref.
#   9  ONE masked 6-sample file (GT:AD, masked = "./.:.", INFO/AF from the KI normals, no GL/PL) used by
#      all three methods. Replaces the per-file cut copies in OUT/input/ (commit 4bf031a).
#   10 Panels A / B / C.
#   11 Thinning: greedy minimum distance keeping the best minimum depth (replaces fixed bins POS %/% SPACING,
#      which allowed two sites 1 bp apart across a bin boundary); primary + sensitivity distances.
#   12 somalier extract once per panel x thinning on the masked file (6 .somalier files each, checked);
#      replaces one extract per FreeBayes file.
#   13 somalier relate: -s sites (INFO/AF), -d MIN_DP, --min-ab HET_MIN, SOMALIER_AB_HOM_CUTOFF set so its
#      hom calls match our bins; -u is never used.
#   14 Crash fix: somalier v0.3.5 pairs.tsv has no hom_concordance column (it has `concordance`, a different
#      quantity). Columns are checked with stopifnot; ibs0/n and ibs2/n use somalier's n.
#   15 bcftools gtcheck on the same masked file, options checked against the installed version's usage,
#      output parsed from that version's own column description (1.4.1: ERR lines).
#   16 Custom IBS0 in R: joint sites, missing rate, IBS0, IBS2, homozygous concordance, 95% exact CIs.
#   17 Matching: IBS0 rate first, IBS2 rate second; ambiguity flags; all 3! one-to-one assignments scored;
#      labelled pairing and HYPOTHESIS tested; labelled-pair labels from KI_NORMALS (replaces the hard-coded
#      list); normal-normal and tumour-tumour pairs kept as reference.
#   18 Plots: pairwise heatmaps (6 metrics), IBS0 bars per tumour (replace the relatedness bars), IBS0 vs
#      IBS2, robustness tiles, chromosomal distribution of panel B/C sites and IBS0 sites (+ knock-in locus).
#   19 Tables: filter_counts, panel_sizes, pairwise_metrics, assignment_scores, verdict, robustness, ...
#   20 RUN_LOG.md: this changelog, assumptions, limitations, deviations, denominators, tool help files.
#   21 Header corrected (relatedness formula, see above).
# 2026-11-25 revision (review): the per-pair FreeBayes ascertainment forces every correctly paired animal
# non-ref at every retained site, so IBS0 cannot confirm correct labels; the old --test only simulated a swap.
#   22 Primary metric: discordance at normal-homozygous sites (tumour-normal), symmetric either-hom version
#      for normal-normal / tumour-tumour pairs; exact CIs. IBS0 kept as a secondary method everywhere.
#   23 Paired exact sign tests (binom.test on discordant-only counts) for each tumour and pair of candidate
#      normals (paired_tests.tsv). Primary ambiguity = paired test vs the second best p >= PAIRED_ALPHA, or
#      < MIN_JOINT_SITES informative sites, or underpowered combination; the old CI-overlap flag is kept.
#   24 10 Mb block bootstrap (BOOT_BLOCK, BOOT_N, BOOT_SEED) of D = [s(T1,N1)+s(T2,N2)] - [s(T1,N2)+s(T2,N1)]
#      and of the support for the winning 3! assignment; replaces the z~2 rule as "clearly_separated"
#      (z~2 kept as separated_z2).
#   25 Leave-one-chromosome-out jackknife of the verdict and the labelled-vs-hypothesis gap (jackknife.tsv).
#   26 Discordant segments per tumour-normal pair (segments.tsv; "X segments on Y chromosomes" in the
#      verdict); informative / discordant sites within +-KI_FLANK of KNOCKIN_LOCUS (ki_locus.tsv).
#   27 Artefact filters: MAX_DP_FACTOR x per-sample median depth (all records read), optional MASK_BED,
#      count of panel-A sites with all three KI normals het (artefact indicator) - in filter_counts.tsv.
#   28 Thinning order: smallest max |log2(depth / sample median)| first, seeded random tie-break (THIN_SEED),
#      then chromosome and position (replaces "best minimum depth first", which prefers duplications).
#   29 FreeBayes PASS records (provider's somatic prefilter) removed from the germline panel (EXCLUDE_FB_PASS).
#   30 Ascertainment sensitivity: T1/T2 swap test on sites required only in files 1 and 2 (extra verdict rows).
#   31 MIN_DP sensitivity (MIN_DP_SENS = 15, 20) of the primary verdict (mindp_sensitivity.tsv).
#   32 gtcheck: bcftools >= 1.11 run with -u GT,GT (+ -E 0 if listed); output parsed from the column names
#      in the version's own header; 1.4.1 -G 1 / ERR kept; DEVIATION if neither works.
#   33 somalier: relatedness removed from the heatmaps (column relatedness_exploratory in
#      pairwise_metrics.tsv); somalier n must equal the custom n_joint (masked AD "." read as unknown):
#      stop in --test, DEVIATION otherwise.
#   34 --test: FreeBayes-like emission (record only if the pair is non-ref, rare stray hom-ref records),
#      inbred-like animals (segments of donor-strain haplotype), paralog artefacts with inflated depth,
#      somatic PASS calls; two scenarios: "swap" (must be detected) and "null" (labels must be confirmed;
#      IBS0 ~0 for every pair documented); new assertions (artefact filter, jackknife/bootstrap/paired outputs).
#   35 RUN_LOG.md: new tables and definitions; limitations on ascertainment asymmetry, linkage, and that
#      genetics cannot tell whether the tumour or the normal label is wrong.
# 2026-11-25 second review:
#   36 EXCLUDE_FB_PASS default FALSE (supersedes the default of CHANGE 29): the PASS prefilter (normal 0/0,
#      tumour 0/1) selects germline differences between mice in a swapped pair, so removing PASS sites
#      deletes swap evidence and makes the site set depend on the true labelling (--test: 2,699 removed in
#      the swap vs 462 in the null scenario); tumour-private somatic sites are already excluded because a KI
#      normal must be non-ref. The opposite setting is run as a sensitivity analysis (pass_sensitivity.tsv).
#   37 Jackknife: column ambiguity_flipped and a DEVIATION when dropping a chromosome flips a tumour's
#      ambiguity flag, not only when its best normal changes.
#   38 paired_tests.tsv: n_hom_j / n_hom_k = common sites where that normal is homozygous (its chance to
#      show discordance).
#   39 RUN_LOG notes: panels B/C expected to be underpowered when most animals are forced non-ref; the
#      max-depth filter does not catch paralog artefacts at normal depth (MASK_BED does).
# END CHANGELOG

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
DIR_A    <- file.path(REPO, "data/mouse/wgs/tumor_vs_normal/snv_indel")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_somalier")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_somalier.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
SOMALIER <- Sys.which("somalier"); if (!nzchar(SOMALIER)) SOMALIER <- "/common/WORK/pstancl/tools/somalier/somalier"
# GRCm38 FASTA (+ .fai) with the same contig names as the VCFs ("1", "2", ...)
FASTA    <- "/common/WORK/pstancl/references/cross_species/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa"

TUMOURS    <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1", HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
                HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing
# [CHANGE 2] hypothesis tested explicitly (from the IBS0 analysis), never assumed
HYPOTHESIS <- c(HCC_1 = "KI2", HCC_2 = "KI1", HCC_3 = "KI3")

MIN_QUAL <- 30
MIN_DP   <- 30            # REF+ALT reads, every one of the six samples
# [CHANGE 2] genotype bins on VAF = ALT / (REF + ALT); anything between the bins is ambiguous
HOM_REF_MAX <- 0.05
HET_MIN     <- 0.35
HET_MAX     <- 0.65
HOM_ALT_MIN <- 0.95
# [CHANGE 2] thinning distances (bp): primary + sensitivity (replaces SITE_SPACING = 20000)
THIN_PRIMARY <- 20000
THIN_SENS    <- c(1e5, 5e5, 1e6)
MIN_PANEL_SITES <- 50     # fewer sites in a panel x thinning -> marked underpowered (not stopped)
MIN_JOINT_SITES <- 50     # fewer jointly callable sites for a tumour-normal pair -> ambiguous
# [CHANGE 2] optional knock-in locus (chrom, position in bp); NA = skipped. Flank used for "near".
KNOCKIN_LOCUS <- list(chrom = NA_character_, pos = NA_real_)
KNOCKIN_FLANK <- 10e6
# [CHANGE 26] flank for the knock-in locus table (informative / discordant sites near the knock-in)
KI_FLANK <- 20e6
# [CHANGE 27] artefact filters
MAX_DP_FACTOR <- 2.5      # drop a site if any sample's REF+ALT > 2.5 x that sample's median REF+ALT
MASK_BED      <- NA_character_   # optional BED of low-mappability / segmental-duplication regions; NA = none
# [CHANGE 28] thinning tie-break seed
THIN_SEED <- 1
# [CHANGE 29] [CHANGE 36] remove FreeBayes PASS records (provider's somatic prefilter) from the germline panel?
# Default FALSE: in a swapped pair the PASS records are largely germline differences between the mice (swap
# evidence). The opposite setting is always run as a sensitivity analysis (pass_sensitivity.tsv).
EXCLUDE_FB_PASS <- FALSE
# [CHANGE 23] paired sign test level for "ambiguous"
PAIRED_ALPHA <- 0.05
# [CHANGE 24] block bootstrap: block length (bp), replicates, seed
BOOT_BLOCK <- 10e6
BOOT_N     <- 1000
BOOT_SEED  <- 1
# [CHANGE 31] depth sensitivity of the primary verdict
MIN_DP_SENS <- c(15, 20)

# [CHANGE 1] synthetic test: TEST_MODE <- TRUE here, or `Rscript ... --test` (runs scenarios "swap" and "null")
TEST_MODE <- FALSE

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# Console as wide as the terminal (falls back to 80 when there is no terminal, e.g. batch Rscript)
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
started <- Sys.time()

# [CHANGE 1] [CHANGE 34] test mode: own output folder; synthetic inputs are generated there in block 2.
# `--test` alone runs this script twice (scenario "swap", then "null") and reports both.
args <- commandArgs(trailingOnly = TRUE)
if ("--test" %in% args) TEST_MODE <- TRUE
SCENARIO <- sub("^--scenario=", "", grep("^--scenario=", args, value = TRUE))
script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
SCRIPT_RUN <- if (length(script_arg) && file.exists(script_arg[1])) normalizePath(script_arg[1]) else SCRIPT
if (TEST_MODE) {
  base <- if (dir.exists(REPO)) REPO else getwd()
  TEST_BASE <- Sys.getenv("SOMALIER_TEST_OUT", file.path(base, "results/mouse_wgs/mouse_wgs_somalier_TEST"))
  if (!length(SCENARIO)) {
    dir.create(TEST_BASE, recursive = TRUE, showWarnings = FALSE)
    status <- c()
    for (sc in c("swap", "null")) {
      message("\n========== TEST scenario: ", sc, " ==========")
      status[sc] <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(SCRIPT_RUN), "--test", paste0("--scenario=", sc)))
    }
    print(status)
    if (any(status != 0)) stop("TEST_MODE: scenario(s) failed: ", paste(names(status)[status != 0], collapse = ", "))
    message("TEST_MODE: both scenarios passed (outputs in ", TEST_BASE, "/swap and /null)")
    quit(save = "no", status = 0)
  }
  stopifnot(SCENARIO %in% c("swap", "null"))
  OUT   <- file.path(TEST_BASE, SCENARIO)
  DIR_A <- file.path(OUT, "test_input")
  FASTA <- file.path(TEST_BASE, "test_genome.fa")          # shared by both scenarios (written once)
  KNOCKIN_LOCUS <- list(chrom = "8", pos = 10e6)
  KNOCKIN_FLANK <- 2e6
  KI_FLANK      <- 5e6
  message("TEST_MODE: synthetic data, scenario ", SCENARIO, ", output ", OUT)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# [CHANGE 3] run-generated folders/files are cleared each run (RUN_LOG.md is kept); also the files the
# previous version of this script wrote
unlink(file.path(OUT, c("extracted", "input", "masked", "sites", "somalier", "gtcheck", "tool_help")), recursive = TRUE)
unlink(file.path(OUT, c("sites.vcf.gz", "sites.vcf.gz.tbi", "somalier.pairs.tsv", "somalier.samples.tsv",
                        "somalier.groups.tsv", "somalier.html", "tumour_vs_normal.tsv", "sites_used.tsv")))
for (d in c("masked", "sites", "somalier", "gtcheck", "tool_help", "extracted")) dir.create(file.path(OUT, d), showWarnings = FALSE)

# [CHANGE 3] logging helpers
DEVIATIONS <- character()
deviation <- function(...) { msg <- paste0(...); DEVIATIONS <<- c(DEVIATIONS, msg); message("DEVIATION: ", msg) }
FILTER_COUNTS <- data.table()
count_step <- function(step, n, scope = "all files", note = "") {
  FILTER_COUNTS <<- rbind(FILTER_COUNTS, data.table(order = nrow(FILTER_COUNTS) + 1L, scope = scope, step = step, n_sites = n, note = note))
  message(sprintf("[filter] %-14s %-70s %s", scope, step, format(n, big.mark = ",")))
}
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150, limitsize = FALSE)
  invisible(p)
}
short_name <- function(s) fifelse(grepl("KI_HCC_[0-9]", s), sub(".*KI_(HCC_[0-9]).*", "\\1", s), sub(".*_(KI[0-9])$", "\\1", s))
run <- function(cmd) { message("$ ", cmd); st <- system(cmd); if (st != 0) stop("failed (exit ", st, "): ", cmd) }
# stdout + stderr of a command as text, whatever its exit status
capture <- function(bin, args = character()) suppressWarnings(system2(bin, args, stdout = TRUE, stderr = TRUE))

TUMOUR_IDS <- names(TUMOURS)
NORMAL_IDS <- unname(KI_NORMALS)
SAMPLES    <- c(TUMOUR_IDS, NORMAL_IDS)
stopifnot(!anyDuplicated(SAMPLES), setequal(HYPOTHESIS, NORMAL_IDS), setequal(names(HYPOTHESIS), TUMOUR_IDS))

# [CHANGE 1] tools: required for real data; in TEST_MODE a missing tool is a DEVIATION and is skipped
HAVE_BCFTOOLS <- file.exists(BCFTOOLS)
HAVE_SOMALIER <- file.exists(SOMALIER)
if (!TEST_MODE) {
  if (!HAVE_BCFTOOLS) stop("bcftools not found at ", BCFTOOLS)
  if (!HAVE_SOMALIER) stop("somalier not found at ", SOMALIER, " - see the install line in the header")
  if (!file.exists(paste0(FASTA, ".fai"))) stop("FASTA index missing: ", FASTA, ".fai")
} else {
  if (!HAVE_BCFTOOLS) deviation("TEST_MODE: bcftools not found - VCF text parsed in R, gtcheck skipped")
  if (!HAVE_SOMALIER) deviation("TEST_MODE: somalier not found - somalier extract/relate skipped")
}
B <- shQuote(BCFTOOLS)

# [CHANGE 4] versions and help texts (somalier writes "somalier version: x" to stderr on every call)
somalier_version <- if (HAVE_SOMALIER) {
  v <- grep("somalier version", capture(SOMALIER), value = TRUE)[1]
  if (is.na(v)) { deviation("somalier version line not found in its output"); "somalier version: unknown" } else v
} else "somalier not available"
bcftools_version <- if (HAVE_BCFTOOLS) capture(BCFTOOLS, "--version")[1:2] else c("bcftools not available", "")
message(somalier_version, " | ", paste(bcftools_version, collapse = " | "))
help_files <- character()
if (HAVE_SOMALIER) for (sub_cmd in c("extract", "relate")) {
  f <- file.path(OUT, "tool_help", paste0("somalier_", sub_cmd, "_help.txt"))
  writeLines(capture(SOMALIER, c(sub_cmd, "--help")), f); help_files <- c(help_files, f)
}
gtcheck_usage <- character()
if (HAVE_BCFTOOLS) {
  f <- file.path(OUT, "tool_help", "bcftools_version.txt"); writeLines(capture(BCFTOOLS, "--version"), f); help_files <- c(help_files, f)
  # [CHANGE 15] "-h", not no arguments: with no file and stdin not a terminal (Rscript, pipes) 1.4.1 reads
  # stdin instead of printing its usage (vcfgtcheck.c: optind==argc && !isatty -> fname = "-")
  gtcheck_usage <- capture(BCFTOOLS, c("gtcheck", "-h"))  # prints its usage, exit status 1
  f <- file.path(OUT, "tool_help", "bcftools_gtcheck_usage.txt"); writeLines(gtcheck_usage, f); help_files <- c(help_files, f)
  # 1.4.1 prints "Version: 1.4.1 (using htslib 1.4.1)": bcftools and htslib versions are separate numbers
  print(grep("Version|htslib", c(bcftools_version, gtcheck_usage), value = TRUE))
}
if (HAVE_SOMALIER) {
  rel_help <- readLines(file.path(OUT, "tool_help", "somalier_relate_help.txt"))
  for (opt in c("--sites", "--min-depth", "--min-ab", "--output-prefix"))
    if (!any(grepl(opt, rel_help, fixed = TRUE))) deviation("somalier relate has no ", opt, " option in this version")
}

# ---- 2. Germline site selection --------------------------------------------------
# [CHANGE 1] [CHANGE 34] synthetic test data. Three animals A, B, C from an inbred-like colony: mostly the
# reference strain, with donor-strain haplotype segments (a shared pool; each animal carries each segment
# at dosage 0, 1 or 2), a donor segment around the knock-in carried by all animals, and sites fixed
# non-reference in all animals. Plus paralog artefacts (all samples ~3.5x depth, random VAF), tumour LOH,
# ambiguous VAFs, low depth, somatic calls (FILTER PASS like the provider's prefilter), indels /
# multiallelics / low QUAL the filters must drop. FreeBayes-like emission: a record goes in a pair's file
# only if one of the pair's two samples is non-ref (plus rare stray hom-ref records).
#   scenario "swap": HCC_1 is animal B filed with KI1 (A), HCC_2 is A filed with KI2 (B), HCC_3 + KI3 = C
#   scenario "null": every tumour is filed with its own animal's normal
# Chromosomes are 19 x ~20 Mb (segments 4-12 Mb), scaled down from the real genome to keep the test small.
make_test_data <- function(dir, scenario) {
  set.seed(if (scenario == "swap") 20261124 else 20261125)
  dir.create(file.path(dir, "calls"), recursive = TRUE, showWarnings = FALSE)
  chroms <- as.character(1:19); L <- 333333L * 60L        # integer length (a double is written as 2e+07)
  line60 <- "ACGTTGCAACGTGCATTACGGATCCAGTTCAGCATGCATCGATCGTAGCTAGCTTACGAC"
  if (!file.exists(paste0(FASTA, ".fai"))) {              # FASTA: one repeated 60-base line, written once
    fa <- file(FASTA, "w"); fai <- data.table(); off <- 0
    for (ch in chroms) {
      hdr <- paste0(">", ch); writeLines(hdr, fa); off <- off + nchar(hdr) + 1
      writeLines(rep(line60, L / 60L), fa)
      fai <- rbind(fai, data.table(ch, L, as.integer(off), 60L, 61L)); off <- off + (L / 60L) * 61
    }
    close(fa); fwrite(fai, paste0(FASTA, ".fai"), sep = "\t", col.names = FALSE)
  }
  base_at <- function(pos) substr(rep(line60, length(pos)), (pos - 1L) %% 60L + 1L, (pos - 1L) %% 60L + 1L)
  st <- rbindlist(lapply(chroms, function(ch) data.table(CHROM = ch, POS = sort(sample(2:(L - 1), 1500)))))
  st[, REF := base_at(POS)]
  st[, ALT := vapply(REF, function(r) sample(setdiff(c("A", "C", "G", "T"), r), 1), "")]
  st[, key := paste(CHROM, POS, REF, ALT, sep = ":")]
  # donor segments: a pool shared by the colony; dosage per animal
  pool <- data.table(CHROM = sample(chroms, 26, replace = TRUE))[, start := runif(.N, 0, L - 4e6)][, end := pmin(start + runif(.N, 4e6, 12e6), L)]
  dose <- sapply(c("A", "B", "C"), function(x) sample(0:2, nrow(pool), replace = TRUE, prob = c(0.35, 0.4, 0.25)))
  ki <- data.table(CHROM = "8", start = 6e6, end = 14e6)                 # donor flank of the knock-in
  geno_of <- function(animal) {
    g <- integer(nrow(st))
    for (i in seq_len(nrow(pool))) { k <- st$CHROM == pool$CHROM[i] & st$POS > pool$start[i] & st$POS <= pool$end[i]; g[k & g == 0L] <- dose[i, animal] }
    k <- st$CHROM == ki$CHROM & st$POS > ki$start & st$POS <= ki$end; g[k] <- pmax(g[k], c(A = 1L, B = 2L, C = 1L)[[animal]])
    g
  }
  truth <- lapply(c(A = "A", B = "B", C = "C"), geno_of)
  fixed <- runif(nrow(st)) < 0.04                                        # non-ref in every animal
  for (x in names(truth)) truth[[x]][fixed] <- 2L
  artefact <- runif(nrow(st)) < 0.01 & !fixed                            # paralog: inflated depth, random VAF
  who <- if (scenario == "swap") c(HCC_1 = "B", KI1 = "A", HCC_2 = "A", KI2 = "B", HCC_3 = "C", KI3 = "C") else
                                 c(HCC_1 = "A", KI1 = "A", HCC_2 = "B", KI2 = "B", HCC_3 = "C", KI3 = "C")
  long <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1", HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
            HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3", KI1 = "sfb_dikic_2025_05_05_KI1",
            KI2 = "sfb_dikic_2025_05_06_KI2", KI3 = "sfb_dikic_2025_05_07_KI3")
  sim_sample <- function(s) {
    g <- truth[[who[[s]]]]
    som <- rep(FALSE, length(g))
    if (grepl("^HCC", s)) {
      k <- g == 1L & runif(length(g)) < 0.05; g[k] <- sample(c(0L, 2L), sum(k), replace = TRUE)   # LOH
      som <- g == 0L & runif(length(g)) < 0.01                                                     # somatic
    }
    dp <- rpois(length(g), 45); low <- runif(length(g)) < 0.05; dp[low] <- rpois(sum(low), 15)
    dp[artefact] <- rpois(sum(artefact), 160)
    p <- c(0.002, 0.5, 0.998)[g + 1L]
    p[som] <- 0.3
    amb <- runif(length(g)) < 0.03; p[amb] <- sample(c(0.2, 0.8), sum(amb), replace = TRUE)
    p[artefact] <- runif(sum(artefact), 0.3, 0.7)
    ao <- rbinom(length(g), dp, p)
    data.table(g = g, som = som, DP = dp, RO = dp - ao, AO = ao)
  }
  for (k in 1:3) {
    pr <- c(paste0("HCC_", k), paste0("KI", k))
    a <- sim_sample(pr[1]); b <- sim_sample(pr[2])
    stray <- runif(nrow(st)) < 0.005
    keep <- ((a$g > 0L | b$g > 0L | a$som | artefact) | stray) & runif(nrow(st)) > 0.02
    gt_txt <- function(x) fifelse(x$AO / pmax(x$DP, 1) < 0.2, "0/0", fifelse(x$AO / pmax(x$DP, 1) > 0.8, "1/1", "0/1"))
    fmt <- function(x) paste(gt_txt(x), x$DP, paste(x$RO, x$AO, sep = ","), x$RO, x$RO * 30, x$AO, x$AO * 30, "-10,-1,0", sep = ":")
    # provider-like prefilter: normal 0/0, tumour 0/1, depth >= 30 in both -> PASS, else "."
    filt <- fifelse(gt_txt(b) == "0/0" & gt_txt(a) == "0/1" & a$DP >= 30 & b$DP >= 30, "PASS", ".")
    rec <- data.table(st[, .(CHROM, POS, REF, ALT)], ID = ".", QUAL = fifelse(runif(nrow(st)) < 0.03, 15, 100), FILTER = filt,
                      INFO = "TYPE=snp", FORMAT = "GT:DP:AD:RO:QR:AO:QA:GL", S1 = fmt(a), S2 = fmt(b))[keep]
    junk <- rec[sample(.N, 30)][, `:=`(ALT = fifelse(seq_len(.N) %% 2 == 0, paste0(ALT, "T"), paste0(ALT, ",G")), INFO = "TYPE=other")]
    rec <- rbind(rec, junk)[order(as.integer(CHROM), POS)]
    setcolorder(rec, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2"))
    f <- file.path(dir, "calls", paste0(long[[pr[1]]], "_vs_", long[[pr[2]]], ".freebayes.ann.norm.vcf"))
    writeLines(c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s,length=%d>", chroms, L),
                 '##FILTER=<ID=PASS,Description="All filters passed">',
                 '##INFO=<ID=TYPE,Number=A,Type=String,Description="type">',
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                 '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
                 '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
                 '##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Ref count">',
                 '##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Ref quality sum">',
                 '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alt count">',
                 '##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Alt quality sum">',
                 '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods">',
                 paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", long[[pr[1]]], long[[pr[2]]]), collapse = "\t")), f)
    fwrite(rec, f, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE)
  }
  tt <- data.table(artefact = artefact, fixed = fixed)[, key := st$key]   # (key = ... inside data.table() would set a key)
  fwrite(tt, file.path(dir, "truth_sites.tsv"), sep = "\t")
  invisible(TRUE)
}
if (TEST_MODE) { unlink(DIR_A, recursive = TRUE); dir.create(DIR_A, recursive = TRUE); make_test_data(DIR_A, SCENARIO) }

fb_files <- list.files(DIR_A, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf|vcf\\.gz)$")
fb_files <- fb_files[grepl("freebayes", fb_files, ignore.case = TRUE) & !grepl("union|merged", fb_files, ignore.case = TRUE)]
basename(fb_files)
stopifnot(length(fb_files) == length(TUMOURS))
fai <- fread(paste0(FASTA, ".fai"), header = FALSE, select = 1:2, col.names = c("ID", "length"), colClasses = list(character = 1))

# Header of a VCF/BCF: bcftools if available, else plain text (TEST_MODE only)
read_header <- function(path) {
  if (HAVE_BCFTOOLS) return(capture(BCFTOOLS, c("view", "-h", shQuote(path))))
  con <- file(path, "r"); on.exit(close(con)); h <- character()
  repeat { l <- readLines(con, n = 1000); h <- c(h, l[startsWith(l, "#")]); if (length(l) < 1000 || !all(startsWith(l, "#"))) break }
  h
}

# [CHANGE 5] [CHANGE 6] one FreeBayes file -> long table (site x sample) with GT, AD, RO, AO
read_fb <- function(path) {
  hdr <- read_header(path)
  fmt_ids <- sub("^##FORMAT=<ID=([^,]+),.*", "\\1", grep("^##FORMAT=<ID=", hdr, value = TRUE))
  if (!"AD" %in% fmt_ids) stop("FORMAT/AD missing in ", basename(path), " - somalier and the VAF bins need it")
  has_roao <- all(c("RO", "AO") %in% fmt_ids)
  if (!has_roao) deviation(basename(path), ": no FORMAT RO/AO - AD == RO,AO check skipped")
  smp_long <- strsplit(tail(grep("^#CHROM", hdr, value = TRUE), 1), "\t", fixed = TRUE)[[1]][-(1:9)]
  smp <- short_name(smp_long)
  # [CHANGE 6] exactly one tumour and one KI normal per file
  if (length(smp) != 2 || sum(smp %chin% TUMOUR_IDS) != 1 || sum(smp %chin% NORMAL_IDS) != 1)
    stop(basename(path), ": expected one tumour and one KI normal, found ", paste(smp_long, collapse = ", "))
  # [CHANGE 6] contig names: VCF autosomes must be in the FASTA index (chr / no-chr)
  vcf_contigs <- sub("^##contig=<ID=([^,>]+).*", "\\1", grep("^##contig=<ID=", hdr, value = TRUE))
  autos <- intersect(c(as.character(1:19), paste0("chr", 1:19)), vcf_contigs)
  if (length(autos) != 19) stop(basename(path), ": expected 19 autosome contigs in the header, found ", length(autos))
  if (!all(autos %chin% fai$ID)) stop(basename(path), ": contig names (", autos[1], "...) not in the FASTA .fai (", fai$ID[1], "...)")
  fields <- c("GT", "AD", if (has_roao) c("RO", "AO"))
  # [CHANGE 29] FILTER read as well, to find the provider's PASS (somatic prefilter) records
  cn <- c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", as.vector(t(outer(smp, fields, paste, sep = "."))))
  if (HAVE_BCFTOOLS) {
    q <- paste0("%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER[", paste0("\\t%", fields, collapse = ""), "]\\n")
    cmd <- paste(B, "view -v snps -m2 -M2 -i", shQuote(paste0("QUAL>=", MIN_QUAL)), "-t", paste(autos, collapse = ","),
                 shQuote(path), "|", B, "query -f", shQuote(q))
    x <- fread(cmd = cmd, sep = "\t", header = FALSE, col.names = cn, colClasses = "character")
  } else {
    # [CHANGE 1] TEST_MODE without bcftools: same filters in R on the VCF text
    v <- fread(path, skip = "#CHROM", sep = "\t", header = TRUE, colClasses = "character")
    setnames(v, 1, "CHROM")
    v <- v[CHROM %chin% autos & nchar(REF) == 1 & grepl("^[ACGT]$", ALT) & as.numeric(QUAL) >= MIN_QUAL]
    fk <- strsplit(v$FORMAT[1], ":", fixed = TRUE)[[1]]
    x <- v[, .(CHROM, POS, REF, ALT, QUAL, FILTER)]
    for (j in seq_along(smp)) {
      parts <- tstrsplit(v[[9 + j]], ":", fixed = TRUE)
      for (fld in fields) set(x, j = paste(smp[j], fld, sep = "."), value = parts[[match(fld, fk)]])
    }
  }
  count_step(paste0("biallelic SNV, autosome, QUAL>=", MIN_QUAL, " (records read)"), nrow(x), scope = paste(smp, collapse = "+"))
  rbindlist(lapply(smp, function(s) {
    ad <- tstrsplit(x[[paste0(s, ".AD")]], ",", fixed = TRUE)
    d <- data.table(CHROM = x$CHROM, POS = as.integer(x$POS), REF = x$REF, ALT = x$ALT, QUAL = as.numeric(x$QUAL), FILTER = x$FILTER,
                    sample = s, file = basename(path), GT_fb = x[[paste0(s, ".GT")]],
                    REF_N = suppressWarnings(as.numeric(ad[[1]])),
                    ALT_N = suppressWarnings(as.numeric(if (length(ad) > 1) ad[[2]] else NA)),
                    AD_n_values = lengths(strsplit(x[[paste0(s, ".AD")]], ",", fixed = TRUE)))
    if (has_roao) d[, `:=`(RO = suppressWarnings(as.numeric(x[[paste0(s, ".RO")]])), AO = suppressWarnings(as.numeric(x[[paste0(s, ".AO")]])))]
    else d[, `:=`(RO = NA_real_, AO = NA_real_)]
    d
  }))
}
ad <- rbindlist(lapply(fb_files, read_fb))
# [CHANGE 6] the six expected short names, each in exactly one file
smp_files <- unique(ad[, .(sample, file)])
if (!setequal(smp_files$sample, SAMPLES) || anyDuplicated(smp_files$sample))
  stop("sample names after short_name(): ", paste(smp_files$sample, collapse = ", "), "; expected ", paste(SAMPLES, collapse = ", "))
print(smp_files)
# [CHANGE 5] AD consistency: biallelic records must have 2 AD values equal to RO,AO
ad_check <- ad[, .(records = .N, AD_not_2_values = sum(AD_n_values != 2, na.rm = TRUE),
                   AD_missing = sum(is.na(REF_N) | is.na(ALT_N)),
                   AD_ne_RO_AO = sum(!is.na(RO) & (REF_N != RO | ALT_N != AO), na.rm = TRUE)), by = .(file, sample)]
print(ad_check)
if (any(ad_check$AD_ne_RO_AO > 0)) deviation("AD differs from RO,AO in ", sum(ad_check$AD_ne_RO_AO), " sample records (AD used; see ad_check.tsv)")
ad[, AD_n_values := NULL]

# A duplicated site within one file (rare): keep the record with the higher depth
ad[, DPA := REF_N + ALT_N]
n_before <- nrow(ad)
ad <- unique(ad[order(-DPA)], by = c("CHROM", "POS", "REF", "ALT", "sample"))
if (nrow(ad) < n_before) message(n_before - nrow(ad), " duplicated sample records dropped (higher depth kept)")

# [CHANGE 27] per-sample median REF+ALT over all records read (for the max-depth filter and thinning order)
MED <- ad[, .(med = median(DPA, na.rm = TRUE)), by = sample][, setNames(med, sample)]
print(MED)

# [CHANGE 8] genotype bins: 0 hom-ref, 1 het, 2 hom-alt, NA = low depth / ambiguous / absent.
# vbin is the VAF bin alone; status adds the MIN_DP depth rule (low_depth)
ad[, VAF := fifelse(DPA > 0, ALT_N / DPA, NA_real_)]
ad[, vbin := fcase(is.na(VAF), "ambiguous", VAF <= HOM_REF_MAX, "hom_ref", VAF >= HET_MIN & VAF <= HET_MAX, "het",
                   VAF >= HOM_ALT_MIN, "hom_alt", default = "ambiguous")]
ad[, status := fifelse(is.na(DPA) | DPA < MIN_DP, "low_depth", vbin)]
ad[, g := fcase(status == "hom_ref", 0L, status == "het", 1L, status == "hom_alt", 2L, default = NA_integer_)]

# Site x sample, all six samples (absent from a file = "absent")
ad[, key := paste(CHROM, POS, REF, ALT, sep = ":")]
# [CHANGE 29] FreeBayes PASS = the provider's somatic prefilter (normal 0/0, tumour 0/1, ...); record-level
pass_keys <- unique(ad[FILTER == "PASS", key])
if (!length(pass_keys)) deviation("no FILTER=PASS records in the FreeBayes files: the PASS exclusion (sensitivity or primary) has no effect")
sites_all <- unique(ad[, .(key, CHROM, POS, REF, ALT)])
full <- CJ(key = sites_all$key, sample = SAMPLES)[ad[, .(key, sample, REF_N, ALT_N, DPA, VAF, vbin, status, g)], on = .(key, sample),
                                                  `:=`(REF_N = i.REF_N, ALT_N = i.ALT_N, DPA = i.DPA, VAF = i.VAF, vbin = i.vbin, status = i.status, g = i.g)]
full[is.na(status), `:=`(status = "absent", vbin = "absent")]
count_step("union of sites over the 3 files", nrow(sites_all))

# [CHANGE 27] optional mask of low-mappability / segmental-duplication regions (BED, 0-based starts)
mask_bed <- NULL
if (!is.na(MASK_BED)) {
  if (!file.exists(MASK_BED)) stop("MASK_BED not found: ", MASK_BED)
  mask_bed <- fread(MASK_BED, header = FALSE, select = 1:3, col.names = c("CHROM", "start", "end"), colClasses = list(character = 1))
  mask_bed[, CHROM := sub("^chr", "", CHROM)]
}
in_mask <- function(chrom, pos) {
  if (is.null(mask_bed)) return(rep(FALSE, length(pos)))
  x <- data.table(CHROM = sub("^chr", "", chrom), POS = pos, i = seq_along(pos))
  hit <- x[mask_bed, on = .(CHROM, POS > start, POS <= end), nomatch = NULL, unique(x.i)]
  seq_along(pos) %in% hit
}

# [CHANGE 7] [CHANGE 27] [CHANGE 29] germline site selection as a function (reused by the sensitivity analyses).
# tum = tumours whose pair files are required (their labelled KI normals must be confident); log = write
# each step to filter_counts.tsv under `scope`.
select_sites <- function(min_dp, tum = TUMOUR_IDS, log = FALSE, scope = "all files", exclude_pass = EXCLUDE_FB_PASS) {
  nrm <- unname(KI_NORMALS[tum]); smp <- c(tum, nrm); k <- length(smp)
  f <- full[sample %chin% smp]
  f[, `:=`(low = is.na(DPA) | DPA < min_dp, is_nrm = sample %chin% nrm, med = MED[sample])]
  s <- f[, .(n_present = sum(status != "absent"), n_dp = sum(!low),
             too_deep = any(!is.na(DPA) & DPA > MAX_DP_FACTOR * med),
             dev = suppressWarnings(max(abs(log2(pmax(DPA, 1) / med)), na.rm = TRUE)),
             min_dp = suppressWarnings(min(DPA, na.rm = TRUE)),
             normal_bad = any(is_nrm & (low | !vbin %chin% c("hom_ref", "het", "hom_alt"))),
             normal_nonref = any(is_nrm & !low & vbin %chin% c("het", "hom_alt")),
             normals_all_het = all(vbin[is_nrm] == "het")), by = key]
  s <- sites_all[s, on = "key"]
  s1 <- s[n_present == k]
  s2 <- s1[n_dp == k]
  s3 <- s2[too_deep == FALSE]
  s4 <- s3[!in_mask(CHROM, POS)]
  s5 <- if (exclude_pass) s4[!key %chin% pass_keys] else s4
  s6 <- s5[normal_bad == FALSE]
  s7 <- s6[normal_nonref == TRUE]
  if (log) {
    count_step(paste0("present in the files of ", paste(tum, collapse = ", "), " (all ", k, " samples have a record)"), nrow(s1), scope)
    n_old <- if (k == 6) s[n_dp >= 4, .N] else NA
    count_step(paste0("REF+ALT >= ", min_dp, " in all ", k, " samples"), nrow(s2), scope,
               note = if (k == 6) sprintf("old rule (>= 4 of 6 samples) would keep %s; the new rule removes %s of those", n_old, n_old - nrow(s2)) else "")
    count_step(paste0("+ no sample deeper than ", MAX_DP_FACTOR, " x its median REF+ALT (collapsed duplications)"), nrow(s3), scope,
               note = paste("medians:", paste(names(MED), round(MED), sep = "=", collapse = " ")))
    count_step("+ outside MASK_BED regions", nrow(s4), scope, note = if (is.null(mask_bed)) "no MASK_BED given: step not applied" else MASK_BED)
    count_step("+ not a FreeBayes PASS record (provider's somatic prefilter)", nrow(s5), scope,
               note = if (!exclude_pass) sprintf("EXCLUDE_FB_PASS = FALSE: not applied (%d PASS keys in the files; TRUE run as sensitivity)", length(pass_keys)) else sprintf("%d PASS keys in the files", length(pass_keys)))
    count_step("+ no KI normal ambiguous (VAF in 0.05-0.35 or 0.65-0.95) or low depth", nrow(s6), scope)
    count_step("+ >= 1 KI normal confidently non-ref (= panel A before thinning)", nrow(s7), scope)
    count_step("  (info) of these, all KI normals het - artefact indicator in an inbred colony", s7[normals_all_het == TRUE, .N], scope,
               note = "not a filter")
  }
  s7
}

# Masked genotype table and matrix for selected sites: tumours ambiguous / low depth -> NA, never hom-ref
build_G <- function(sel, tum = TUMOUR_IDS, min_dp = MIN_DP) {
  smp <- c(tum, unname(KI_NORMALS[tum]))
  geno <- full[key %chin% sel$key & sample %chin% smp]
  geno[, g := fifelse(is.na(DPA) | DPA < min_dp, NA_integer_,
                      fcase(vbin == "hom_ref", 0L, vbin == "het", 1L, vbin == "hom_alt", 2L, default = NA_integer_))]
  geno[, masked := is.na(g)]
  G <- dcast(geno, key ~ sample, value.var = "g")[sel[, .(key, CHROM, POS, REF, ALT, min_dp, dev)], on = "key"]
  setcolorder(G, c("key", "CHROM", "POS", "REF", "ALT", "min_dp", "dev", smp))
  G[, chr_n := as.integer(sub("^chr", "", CHROM))]
  setorder(G, chr_n, POS)
  list(G = G, geno = geno)
}

final <- select_sites(MIN_DP, log = TRUE)
if (!nrow(final)) stop("no site passes the germline selection - see filter_counts in the log")
bg <- build_G(final)
geno <- bg$geno; G <- bg$G
# [CHANGE 8] masked genotypes on the final sites: tumours ambiguous -> NA (masked), never hom-ref
masked_summary <- geno[, .(sites = .N, masked = sum(masked), masked_ambiguous = sum(status == "ambiguous"),
                           hom_ref = sum(g == 0, na.rm = TRUE), het = sum(g == 1, na.rm = TRUE), hom_alt = sum(g == 2, na.rm = TRUE)), by = sample]
print(masked_summary)
stopifnot(geno[sample %chin% NORMAL_IDS, !any(masked)])          # normals are confident by construction
G[, AF := round(rowSums(.SD) / (2 * length(NORMAL_IDS)), 4), .SDcols = NORMAL_IDS]   # in-sample, KI normals only

# [CHANGE 10] panels
G[, `:=`(n_min = do.call(pmin, .SD), n_max = do.call(pmax, .SD)), .SDcols = NORMAL_IDS]
G[, `:=`(panel_B = n_min != n_max, panel_C = n_min == 0L & n_max == 2L)]      # B: not all identical; C: 0/0 and 1/1 present
G[, c("n_min", "n_max") := NULL]
panel_sites <- list(A = G$key, B = G[panel_B == TRUE, key], C = G[panel_C == TRUE, key])

# [CHANGE 11] [CHANGE 28] greedy minimum-distance thinning: a site is kept only if no kept site on the same
# chromosome is closer than `dist`. Order: smallest max |log2(depth / sample median)| first (sites with
# typical depth; "best minimum depth first" preferred collapsed duplications), seeded random tie-break,
# then chromosome and position.
thin_greedy <- function(d, dist) {
  set.seed(THIN_SEED)
  d <- copy(d)[, rnd := runif(.N)][order(dev, rnd, chr_n, POS)]
  keep <- logical(nrow(d)); acc <- new.env()
  for (i in seq_len(nrow(d))) {
    ch <- as.character(d$chr_n[i]); p <- d$POS[i]
    v <- acc[[ch]]
    if (is.null(v) || min(abs(v - p)) >= dist) { keep[i] <- TRUE; acc[[ch]] <- c(v, p) }
  }
  d$key[keep]
}
THINS <- c(THIN_PRIMARY, THIN_SENS)
thin_label <- function(x) fifelse(x >= 1e6, paste0(x / 1e6, "Mb"), paste0(x / 1e3, "kb"))
combos <- CJ(panel = c("A", "B", "C"), thin = THINS, sorted = FALSE)
combos[, combo := paste0(panel, "_", thin_label(thin))]
combo_sites <- setNames(lapply(seq_len(nrow(combos)), function(i)
  thin_greedy(G[key %chin% panel_sites[[combos$panel[i]]]], combos$thin[i])), combos$combo)
combos[, n_sites := lengths(combo_sites[combo])]
combos[, panel_unthinned := lengths(panel_sites[panel])]
combos[, underpowered := n_sites < MIN_PANEL_SITES]
combos[, primary := panel == "A" & thin == THIN_PRIMARY]
print(combos)
for (i in seq_len(nrow(combos))) count_step(paste0("panel ", combos$panel[i], ", thinned ", thin_label(combos$thin[i])),
                                            combos$n_sites[i], note = if (combos$underpowered[i]) "UNDERPOWERED" else "")
for (i in which(combos$underpowered)) deviation("panel x thinning ", combos$combo[i], " has ", combos$n_sites[i],
                                                " sites (< MIN_PANEL_SITES = ", MIN_PANEL_SITES, "): results marked underpowered")

# [CHANGE 9] ONE masked 6-sample file: GT:AD, masked = "./.:.", INFO/AF, no GL/PL
gt_str <- function(g) fcase(g == 0L, "0/0", g == 1L, "0/1", g == 2L, "1/1", default = "./.")
geno[, ad_str := fifelse(masked, ".", paste(REF_N, ALT_N, sep = ","))]
AD_W <- dcast(geno, key ~ sample, value.var = "ad_str")
masked_vcf <- file.path(OUT, "masked", "masked_6samples.vcf")
# [CHANGE 9] lengths written as plain integers (a length like "1e+06" is read by bcftools as 1 and breaks the index)
contig_lines <- fai[, sprintf("##contig=<ID=%s,length=%s>", ID, format(as.numeric(length), scientific = FALSE, trim = TRUE))]
writeLines(c("##fileformat=VCFv4.2", contig_lines,
             '##INFO=<ID=AF,Number=A,Type=Float,Description="ALT allele frequency among the KI normal genotypes (in-sample)">',
             '##FILTER=<ID=PASS,Description="All filters passed">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype from VAF bins; ./. = masked (ambiguous or low depth)">',
             '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (REF,ALT) from FreeBayes; . when masked">',
             paste0("##mouse_wgs_somalier=bins hom_ref<=", HOM_REF_MAX, " het ", HET_MIN, "-", HET_MAX, " hom_alt>=", HOM_ALT_MIN, " min_dp ", MIN_DP),
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", SAMPLES), collapse = "\t")), masked_vcf)
body <- G[, .(CHROM, POS, ID = ".", REF, ALT, QUAL = ".", FILTER = "PASS", INFO = paste0("AF=", AF), FORMAT = "GT:AD")]
gA <- AD_W[G[, .(key)], on = "key"]
for (s in SAMPLES) set(body, j = s, value = paste(gt_str(G[[s]]), gA[[s]], sep = ":"))
fwrite(body, masked_vcf, sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE)
masked_file <- masked_vcf
if (HAVE_BCFTOOLS) {
  masked_file <- file.path(OUT, "masked", "masked_6samples.bcf")
  run(paste(B, "view -Ob -o", shQuote(masked_file), shQuote(masked_vcf)))
  run(paste(B, "index -f", shQuote(masked_file)))
  unlink(masked_vcf)
}
# [CHANGE 8] check: a masked genotype is never written as hom-ref
stopifnot(all(body[, unlist(lapply(.SD, function(x) !startsWith(x, "0/0:."))), .SDcols = SAMPLES]))

# Per combo: somalier sites VCF (with INFO/AF) and a gtcheck targets file, both subsets of the masked file
write_combo_sites <- function(cb) {
  d <- G[key %chin% combo_sites[[cb]]]
  vcf <- file.path(OUT, "sites", paste0("sites_", cb, ".vcf"))
  writeLines(c("##fileformat=VCFv4.2", contig_lines,
               '##INFO=<ID=AF,Number=A,Type=Float,Description="ALT allele frequency among the KI normal genotypes (in-sample)">',
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), vcf)
  fwrite(d[, .(CHROM, POS, ID = ".", REF, ALT, QUAL = ".", FILTER = "PASS", INFO = paste0("AF=", AF))], vcf,
         sep = "\t", col.names = FALSE, append = TRUE, quote = FALSE)
  fwrite(d[, .(CHROM, POS)], file.path(OUT, "sites", paste0("targets_", cb, ".tsv")), sep = "\t", col.names = FALSE)
  if (HAVE_BCFTOOLS) {
    run(paste(B, "view -Oz -o", shQuote(paste0(vcf, ".gz")), shQuote(vcf)))
    run(paste(B, "index -t -f", shQuote(paste0(vcf, ".gz"))))
    unlink(vcf); return(paste0(vcf, ".gz"))
  }
  vcf
}
combos[, sites_file := vapply(combo, function(cb) if (length(combo_sites[[cb]])) write_combo_sites(cb) else NA_character_, "")]

# ---- 3. somalier extract ---------------------------------------------------------
# [CHANGE 12] per combo, on the masked file; exactly the 6 samples must come out
somalier_extract <- function(cb) {
  d <- file.path(OUT, "extracted", cb); unlink(d, recursive = TRUE); dir.create(d, recursive = TRUE)
  run(paste(shQuote(SOMALIER), "extract -d", shQuote(d), "-s", shQuote(combos[combo == cb, sites_file]),
            "-f", shQuote(FASTA), shQuote(masked_file)))
  got <- list.files(d, pattern = "\\.somalier$")
  if (!setequal(got, paste0(SAMPLES, ".somalier")) || length(got) != 6)
    stop("somalier extract ", cb, ": expected 6 files ", paste0(SAMPLES, ".somalier", collapse = ", "), "; got ", paste(got, collapse = ", "))
  file.path(d, got)
}
som_files <- list()
if (HAVE_SOMALIER) for (cb in combos[n_sites > 0, combo]) som_files[[cb]] <- somalier_extract(cb)
for (cb in combos[n_sites == 0, combo]) deviation("panel x thinning ", cb, " is empty: no method run on it")

# ---- 4. somalier relate, bcftools gtcheck, custom metrics ---------------------------
pair_dt <- as.data.table(t(combn(SAMPLES, 2)))
setnames(pair_dt, c("s1", "s2"))
pair_type <- function(a, b) {
  ta <- a %chin% TUMOUR_IDS; tb <- b %chin% TUMOUR_IDS
  fcase(ta & tb, "tumour-tumour", !ta & !tb, "normal-normal", default = "tumour-normal")
}
# [CHANGE 17] labelled / hypothesis flags from KI_NORMALS and HYPOTHESIS (no hard-coded pair list)
is_pair <- function(a, b, map) (a %chin% names(map) & map[a] == b) | (b %chin% names(map) & map[b] == a)
ci_lo <- function(x, n) mapply(function(x, n) if (!is.na(n) && n > 0) binom.test(x, n)$conf.int[1] else NA_real_, x, n)
ci_hi <- function(x, n) mapply(function(x, n) if (!is.na(n) && n > 0) binom.test(x, n)$conf.int[2] else NA_real_, x, n)

# [CHANGE 13] [CHANGE 14] somalier relate; SOMALIER_AB_HOM_CUTOFF makes its hom calls match our bins
# (v0.3.5 source: hom-ref if AB < cutoff, hom-alt if AB > 1 - cutoff, het if min-ab <= AB <= 1 - min-ab,
# else unknown; depth < -d unknown; default cutoff 0.01). The tiny offset keeps AB = 0.05 / 0.95 as hom.
SOM_PAIRS_COLS <- c("sample_a", "sample_b", "relatedness", "ibs0", "ibs2", "concordance", "hets_a", "hets_b", "hets_ab",
                    "shared_hets", "hom_alts_a", "hom_alts_b", "shared_hom_alts", "n", "x_ibs0", "x_ibs2", "expected_relatedness")
AB_HOM_CUTOFF <- HOM_REF_MAX + 1e-6
stopifnot(abs((1 - HOM_ALT_MIN) - HOM_REF_MAX) < 1e-9, abs((1 - HET_MAX) - HET_MIN) < 1e-9)   # symmetric bins needed
somalier_relate <- function(cb) {
  d <- file.path(OUT, "somalier", cb); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(d, "somalier")
  run(paste0("SOMALIER_AB_HOM_CUTOFF=", format(AB_HOM_CUTOFF, scientific = FALSE), " ", shQuote(SOMALIER), " relate",
             " -s ", shQuote(combos[combo == cb, sites_file]), " -d ", MIN_DP, " --min-ab ", HET_MIN,
             " -o ", shQuote(prefix), " ", paste(shQuote(som_files[[cb]]), collapse = " ")))
  p <- fread(paste0(prefix, ".pairs.tsv"))
  setnames(p, 1, sub("^#", "", names(p)[1]))
  message("somalier pairs.tsv columns (", cb, "): ", paste(names(p), collapse = " "))
  if (!all(SOM_PAIRS_COLS %in% names(p))) stop("somalier pairs.tsv lacks: ", paste(setdiff(SOM_PAIRS_COLS, names(p)), collapse = ", "))
  # [CHANGE 33] relatedness kept only as an exploratory column
  p[, .(combo = cb, s1 = sample_a, s2 = sample_b, method = "somalier", n_joint = n, ibs0, ibs2,
        relatedness_exploratory = relatedness, somalier_concordance = concordance, hets_a, hets_b, shared_hets)]
}
som <- if (HAVE_SOMALIER) rbindlist(lapply(names(som_files), somalier_relate)) else data.table()

# [CHANGE 15] [CHANGE 32] gtcheck: only options listed by the installed version's usage.
# bcftools >= 1.11: -u GT,GT (genotypes only) + -E 0 if listed (discordance = number of differing genotypes);
# 1.4.1: -G 1 (ERR lines). Output parsed from the column names in the version's own header line.
has_gt_opt <- function(o) any(grepl(o, gtcheck_usage, fixed = TRUE))
GT_ARGS <- if (!HAVE_BCFTOOLS || !has_gt_opt("-T, --targets-file")) NULL else
  if (has_gt_opt("-u, --use")) c("-u GT,GT", if (has_gt_opt("-E, --error-probability")) "-E 0", if (has_gt_opt("--no-HWE-prob")) "--no-HWE-prob") else
  if (has_gt_opt("-G, --GTs-only")) "-G 1" else NULL
if (HAVE_BCFTOOLS && is.null(GT_ARGS)) deviation("bcftools gtcheck usage lists neither -u/--use nor -G/--GTs-only (or no -T): gtcheck skipped (usage in tool_help/)")
gtcheck_run <- function(cb) {
  raw <- file.path(OUT, "gtcheck", paste0("gtcheck_", cb, ".txt"))
  run(paste(B, "gtcheck", paste(GT_ARGS, collapse = " "), "-T", shQuote(file.path(OUT, "sites", paste0("targets_", cb, ".tsv"))),
            shQuote(masked_file), ">", shQuote(raw), "2>", shQuote(paste0(raw, ".stderr"))))
  l <- readLines(raw)
  hdr <- grep("^#\\s*(DCv2|DC)\\b.*\\[2\\]", l, value = TRUE, perl = TRUE)[1]
  if (!is.na(hdr)) {
    tag <- sub("^#\\s*(DCv2|DC)\\b.*", "\\1", hdr, perl = TRUE)
    m <- regmatches(hdr, gregexpr("\\[([0-9]+)\\]([^\\t\\[]+)", hdr, perl = TRUE))[[1]]
    nm <- tolower(trimws(sub("[,;:]\\s*$", "", trimws(sub("^\\[[0-9]+\\]", "", m)))))
    idx <- setNames(as.integer(sub("^\\[([0-9]+)\\].*", "\\1", m)), nm)
    find <- function(p) { i <- grep(p, names(idx)); if (length(i)) idx[[i[1]]] else NA_integer_ }
    iq <- find("query"); ig <- find("genotyped"); id <- find("^discordance"); isites <- find("sites compared|number of sites")
    if (anyNA(c(iq, ig, id, isites))) { deviation("gtcheck ", cb, ": header columns not recognised (", hdr, ") - raw kept in ", raw); return(NULL) }
    rows <- grep(paste0("^", tag, "\t"), l, value = TRUE)
    if (!length(rows)) { deviation("gtcheck ", cb, ": no ", tag, " lines - raw kept in ", raw); return(NULL) }
    x <- as.data.table(do.call(rbind, strsplit(rows, "\t", fixed = TRUE)))
    x <- data.table(s1 = x[[iq]], s2 = x[[ig]], disc = as.numeric(x[[id]]), n = as.numeric(x[[isites]]))[s1 != s2]
    if (!has_gt_opt("-E, --error-probability")) deviation("gtcheck ", cb, ": no -E option, discordance is the default error-model score")
    x <- unique(x[, c("s1", "s2") := .(pmin(s1, s2), pmax(s1, s2))], by = c("s1", "s2"))
    return(x[, .(combo = cb, s1, s2, method = "gtcheck", n_joint = n, gtcheck_discordance = disc / n, gtcheck_n_discordant = round(disc))])
  }
  # 1.4.1 cross-check output: "# ERR, error rate\t[2]Pairwise error rate\t[3]Number of sites compared\t[4]Sample i\t[5]Sample j"
  if (any(grepl("^# ERR.*Pairwise error rate.*Number of sites compared.*Sample i.*Sample j", l))) {
    e <- grep("^ERR\t", l, value = TRUE)
    if (!length(e)) { deviation("gtcheck ", cb, ": no ERR lines"); return(NULL) }
    x <- as.data.table(tstrsplit(e, "\t", fixed = TRUE))[, .(err = as.numeric(V2), n = as.integer(V3), s1 = V4, s2 = V5)]
    return(x[, .(combo = cb, s1, s2, method = "gtcheck", n_joint = n, gtcheck_discordance = err, gtcheck_n_discordant = round(err * n))])
  }
  deviation("gtcheck ", cb, ": output format not recognised - raw output kept in ", raw, ", not parsed"); NULL
}
gtc <- if (!is.null(GT_ARGS)) rbindlist(lapply(combos[n_sites > 0, combo], gtcheck_run)) else data.table()

# [CHANGE 16] [CHANGE 22] custom metrics on the masked genotypes, per sample pair.
# Tumour-normal: informative = normal homozygous (0/0 or 1/1) and tumour unmasked; discordant = tumour
# genotype != normal genotype. Normal-normal / tumour-tumour: symmetric version, informative = both
# unmasked and at least one homozygous; discordant = genotypes differ.
pair_counts <- function(d, s1, s2) {
  a <- d[[s1]]; b <- d[[s2]]; j <- !is.na(a) & !is.na(b)
  if (xor(s1 %chin% TUMOUR_IDS, s2 %chin% TUMOUR_IDS)) {
    nrm <- if (s1 %chin% NORMAL_IDS) a else b; tum <- if (s1 %chin% NORMAL_IDS) b else a
    inf <- j & nrm != 1L; dsc <- inf & tum != nrm; def <- "normal-hom (tumour-normal)"
  } else { inf <- j & (a != 1L | b != 1L); dsc <- inf & a != b; def <- "either-hom (symmetric)" }
  hb <- j & a != 1L & b != 1L
  list(n_joint = sum(j), ibs0 = sum(j & abs(a - b) == 2L), ibs2 = sum(j & a == b), n_hom_both = sum(hb),
       hom_identical = sum(hb & a == b), n_inf = sum(inf), disc = sum(dsc), disc_def = def,
       missing_s1 = mean(is.na(a)), missing_s2 = mean(is.na(b)))
}
custom_metrics <- function(cb, Gx = G, keys = combo_sites[[cb]]) {
  d <- Gx[key %chin% keys]
  rbindlist(lapply(seq_len(nrow(pair_dt)), function(i)
    as.data.table(c(list(combo = cb, s1 = pair_dt$s1[i], s2 = pair_dt$s2[i], method = "custom", n_sites = nrow(d)),
                    pair_counts(d, pair_dt$s1[i], pair_dt$s2[i])))))
}
cust <- rbindlist(lapply(combos[n_sites > 0, combo], custom_metrics))

# One long table for every method; rates and exact 95% CIs on counts
canon <- function(x) { if (!nrow(x)) return(x); sw <- match(x$s1, SAMPLES) > match(x$s2, SAMPLES)
  x[sw, c("s1", "s2") := .(s2, s1)]; x }
cust <- canon(cust); som <- canon(som); gtc <- canon(gtc)
add_rates <- function(x) { if (!nrow(x)) return(x)
  x[, `:=`(ibs0_rate = ibs0 / n_joint, ibs2_rate = ibs2 / n_joint, ibs0_lo = ci_lo(ibs0, n_joint), ibs0_hi = ci_hi(ibs0, n_joint))] }
cust <- add_rates(cust); som <- add_rates(som)
cust[, `:=`(hom_concordance = fifelse(n_hom_both > 0, hom_identical / n_hom_both, NA_real_),
            disc_rate = fifelse(n_inf > 0, disc / n_inf, NA_real_), disc_lo = ci_lo(disc, n_inf), disc_hi = ci_hi(disc, n_inf))]
if (nrow(gtc)) gtc[, `:=`(gt_lo = ci_lo(gtcheck_n_discordant, n_joint), gt_hi = ci_hi(gtcheck_n_discordant, n_joint))]
pairwise <- rbindlist(list(cust, som, gtc), fill = TRUE)
for (cl in c("ibs0_rate", "ibs2_rate", "ibs0_lo", "ibs0_hi", "hom_concordance", "relatedness_exploratory", "somalier_concordance",
             "gtcheck_discordance", "gt_lo", "gt_hi", "n_inf", "disc", "disc_rate", "disc_lo", "disc_hi")) if (!cl %in% names(pairwise)) pairwise[, (cl) := NA_real_]
pairwise[, `:=`(type = pair_type(s1, s2), labelled_pair = is_pair(s1, s2, KI_NORMALS), hypothesis_pair = is_pair(s1, s2, HYPOTHESIS))]
pairwise <- combos[, .(combo, panel, thin, underpowered, primary)][pairwise, on = "combo"]
pairwise[primary == TRUE & method == "custom", .(s1, s2, type, n_inf, disc_rate = round(disc_rate, 4), n_joint, ibs0,
                                                 ibs0_rate = round(ibs0_rate, 4), ibs2_rate = round(ibs2_rate, 3), labelled_pair)][order(type, disc_rate)]

# [CHANGE 33] masked AD "." must be read by somalier as unknown: its n must equal the custom n_joint
# (same thresholds; a masked genotype read as hom-ref would raise somalier's n)
som_n_check <- if (nrow(som)) merge(cust[, .(combo, s1, s2, n_custom = n_joint)], som[, .(combo, s1, s2, n_somalier = n_joint)],
                                    by = c("combo", "s1", "s2")) else data.table()
if (nrow(som_n_check)) {
  bad <- som_n_check[n_custom != n_somalier]
  message("somalier n vs custom n_joint: ", nrow(som_n_check) - nrow(bad), " / ", nrow(som_n_check), " pairs equal")
  if (nrow(bad)) {
    msg <- paste0("somalier n differs from the custom n_joint for ", nrow(bad), " pair x combo rows (max |diff| ",
                  max(abs(bad$n_custom - bad$n_somalier)), "): masked genotypes may not be read as unknown")
    if (TEST_MODE) stop(msg) else deviation(msg)
  }
}

# [CHANGE 23] paired exact sign tests, tumour T vs candidate normals Nj, Nk, on sites where T, Nj and Nk are
# all unmasked and at least one of the two normals is homozygous. T counts as discordant with a normal
# only where THAT normal is homozygous (the primary metric, judged per normal on common sites):
# b = discordant with Nj only, c = with Nk only; exact two-sided binom.test(b, b + c).
# The review spec asked for sites where BOTH normals are homozygous; that version is kept as columns
# (*_bothhom) for comparison. It drops the het-vs-hom sites that carry most of the difference between
# animals of an inbred-like colony (dosage 1 vs 2 in a donor segment), so it is underpowered: in --test it
# gave p = 0.25 for 0/1127 vs 80/938 discordant sites.
paired_counts <- function(d, t, nj, nk) {
  x <- d[[t]]; a <- d[[nj]]; b <- d[[nk]]
  ok <- !is.na(x) & !is.na(a) & !is.na(b) & (a != 1L | b != 1L)
  dj <- ok & a != 1L & x != a; dk <- ok & b != 1L & x != b
  n_b <- sum(dj & !dk); n_c <- sum(!dj & dk)
  bh <- ok & a != 1L & b != 1L; bj <- bh & x != a; bk <- bh & x != b
  b2 <- sum(bj & !bk); c2 <- sum(!bj & bk)
  # [CHANGE 38] n_hom_j / n_hom_k: common sites where that normal is homozygous (its chance to show discordance)
  list(n_common = sum(ok), n_hom_j = sum(ok & a != 1L), n_hom_k = sum(ok & b != 1L), disc_j_only = n_b, disc_k_only = n_c,
       p_value = if (n_b + n_c > 0) binom.test(n_b, n_b + n_c)$p.value else 1,
       n_common_bothhom = sum(bh), disc_j_only_bothhom = b2, disc_k_only_bothhom = c2,
       p_value_bothhom = if (b2 + c2 > 0) binom.test(b2, b2 + c2)$p.value else 1)
}
paired_table <- function(d, cb, tum = TUMOUR_IDS, nrms = NORMAL_IDS) {
  np <- t(combn(nrms, 2))
  rbindlist(lapply(tum, function(t) rbindlist(lapply(seq_len(nrow(np)), function(i)
    as.data.table(c(list(combo = cb, tumour = t, normal_j = np[i, 1], normal_k = np[i, 2]), paired_counts(d, t, np[i, 1], np[i, 2])))))))
}
paired <- rbindlist(lapply(combos[n_sites > 0, combo], function(cb) paired_table(G[key %chin% combo_sites[[cb]]], cb)))
paired_p <- function(tb, cb, t, x, y) {
  r <- tb[combo == cb & tumour == t & ((normal_j == x & normal_k == y) | (normal_j == y & normal_k == x))]
  if (nrow(r)) r$p_value[1] else NA_real_
}

# ---- 5. Matching, assignments, plots ----------------------------------------------
# [CHANGE 17] [CHANGE 22] [CHANGE 23] per tumour and method, normals ranked by score:
#   disc_nhom (primary): normal-hom discordance rate; tie-break IBS0 rate, then higher IBS2 rate
#   custom (IBS0, secondary) and somalier: IBS0 rate, then higher IBS2 rate; gtcheck: discordance
# ambiguous (disc_nhom): paired sign test best vs second p >= PAIRED_ALPHA, or < MIN_JOINT_SITES informative
# sites, or underpowered; other methods: CI overlap, or < MIN_JOINT_SITES, or underpowered
tn <- pairwise[type == "tumour-normal"]
tn[, `:=`(tumour = fifelse(s1 %chin% TUMOUR_IDS, s1, s2), normal = fifelse(s1 %chin% TUMOUR_IDS, s2, s1))]
tn <- rbind(tn, copy(tn[method == "custom"])[, method := "disc_nhom"])
tn[, score := fcase(method == "disc_nhom", disc_rate, method == "gtcheck", gtcheck_discordance, default = ibs0_rate)]
tn[, score_lo := fcase(method == "disc_nhom", disc_lo, method == "gtcheck", gt_lo, default = ibs0_lo)]
tn[, score_hi := fcase(method == "disc_nhom", disc_hi, method == "gtcheck", gt_hi, default = ibs0_hi)]
tn[, n_eff := fifelse(method == "disc_nhom", n_inf, n_joint)]
tn[, tie1 := fcase(method == "disc_nhom", ibs0_rate, method == "gtcheck", 0, default = -ibs2_rate)]
tn[, tie2 := fifelse(method == "disc_nhom", -ibs2_rate, 0)]
setorder(tn, combo, method, tumour, score, tie1, tie2)
robust <- tn[, .(best_normal = normal[1], second_normal = normal[2], best_score = score[1], second_score = score[2],
                 best_hi = score_hi[1], second_lo = score_lo[2], n_eff_best = n_eff[1], underpowered = underpowered[1]),
             by = .(combo, panel, thin, method, tumour)]
robust[, ci_overlap := best_hi >= second_lo]
robust[, paired_p_best_vs_second := fifelse(method == "disc_nhom", mapply(function(cb, t, x, y) paired_p(paired, cb, t, x, y),
                                                                           combo, tumour, best_normal, second_normal), NA_real_)]
robust[, ambiguous := fifelse(method == "disc_nhom", is.na(paired_p_best_vs_second) | paired_p_best_vs_second >= PAIRED_ALPHA,
                              ci_overlap %in% TRUE) | n_eff_best < MIN_JOINT_SITES | underpowered]
robust[, `:=`(labelled_normal = KI_NORMALS[tumour], hypothesis_normal = HYPOTHESIS[tumour])]
robust[, `:=`(matches_labelled = best_normal == labelled_normal, matches_hypothesis = best_normal == hypothesis_normal)]
robust[, call := fifelse(ambiguous, paste0(best_normal, "?"), best_normal)]
METHOD_ORDER <- intersect(c("disc_nhom", "custom", "somalier", "gtcheck"), unique(robust$method))
robust[, method := factor(method, METHOD_ORDER)]
dcast(robust, method + tumour ~ combo, value.var = "call")

# [CHANGE 17] all 3! one-to-one tumour -> normal assignments, scored by the summed score of their 3 pairs
perm3 <- function(v) { r <- list(); for (i in v) for (j in setdiff(v, i)) r[[length(r) + 1]] <- c(i, j, setdiff(v, c(i, j))); r }
PERMS <- perm3(NORMAL_IDS)
perm_label <- function(pm) paste(paste0(TUMOUR_IDS, "-", pm), collapse = " ")
assign_scores <- rbindlist(lapply(split(tn, by = c("combo", "method")), function(d) {
  rbindlist(lapply(PERMS, function(pm) {
    m <- setNames(pm, TUMOUR_IDS)
    sel <- d[normal == m[tumour]]
    data.table(combo = d$combo[1], method = d$method[1], assignment = perm_label(pm),
               total_score = sum(sel$score), var_sum = sum(sel[, score * (1 - score) / pmax(n_eff, 1)]),
               is_labelled = all(m == KI_NORMALS[TUMOUR_IDS]), is_hypothesis = all(m == HYPOTHESIS[TUMOUR_IDS]))
  }))
}))
setorder(assign_scores, combo, method, total_score)
assign_scores[, rank := seq_len(.N), by = .(combo, method)]
assign_scores[, gap_to_next := shift(total_score, -1) - total_score, by = .(combo, method)]
# [CHANGE 24] z~2 rule kept as a column; "clearly_separated" now comes from the block bootstrap (primary only)
assign_scores[, separated_z2 := gap_to_next > 2 * sqrt(var_sum + shift(var_sum, -1)), by = .(combo, method)]

PRIMARY <- combos[primary == TRUE, combo]
prim_d <- G[key %chin% combo_sites[[PRIMARY]]]
TN_PAIRS <- CJ(tumour = TUMOUR_IDS, normal = NORMAL_IDS, sorted = FALSE)
TN_PAIRS[, pid := paste(tumour, normal, sep = "|")]
# site-level informative / discordant vectors for one tumour-normal pair (primary metric)
tn_site <- function(d, t, n) { x <- d[[t]]; y <- d[[n]]; inf <- !is.na(x) & !is.na(y) & y != 1L; list(inf = inf, disc = inf & x != y) }

# [CHANGE 24] block bootstrap on the primary combination (BOOT_BLOCK bp blocks, BOOT_N replicates)
blk <- factor(paste(prim_d$chr_n, floor(prim_d$POS / BOOT_BLOCK)))
INF <- DISC <- matrix(0, nlevels(blk), nrow(TN_PAIRS), dimnames = list(levels(blk), TN_PAIRS$pid))
for (i in seq_len(nrow(TN_PAIRS))) {
  v <- tn_site(prim_d, TN_PAIRS$tumour[i], TN_PAIRS$normal[i])
  INF[, i] <- rowsum(as.numeric(v$inf), blk)[, 1]; DISC[, i] <- rowsum(as.numeric(v$disc), blk)[, 1]
}
set.seed(BOOT_SEED)
W <- t(vapply(seq_len(BOOT_N), function(r) tabulate(sample.int(nlevels(blk), nlevels(blk), replace = TRUE), nlevels(blk)), numeric(nlevels(blk))))
RATE_B <- (W %*% DISC) / (W %*% INF)
rate_full <- colSums(DISC) / colSums(INF)
T1 <- TUMOUR_IDS[1]; T2 <- TUMOUR_IDS[2]; N1 <- KI_NORMALS[[T1]]; N2 <- KI_NORMALS[[T2]]
d_of <- function(r) r[, paste(T1, N1, sep = "|")] + r[, paste(T2, N2, sep = "|")] - r[, paste(T1, N2, sep = "|")] - r[, paste(T2, N1, sep = "|")]
D_boot <- d_of(RATE_B); D_point <- d_of(matrix(rate_full, 1, dimnames = list(NULL, names(rate_full))))
ASSIGN_B <- sapply(PERMS, function(pm) rowSums(RATE_B[, paste(TUMOUR_IDS, pm, sep = "|"), drop = FALSE]))
colnames(ASSIGN_B) <- vapply(PERMS, perm_label, "")
valid_b <- complete.cases(ASSIGN_B)
win_b <- colnames(ASSIGN_B)[apply(ASSIGN_B[valid_b, , drop = FALSE], 1, which.min)]
boot_D <- data.table(T1 = T1, T2 = T2, N1 = N1, N2 = N2, D = D_point, D_lo = quantile(D_boot, 0.025, na.rm = TRUE),
                     D_hi = quantile(D_boot, 0.975, na.rm = TRUE), frac_D_gt_0 = mean(D_boot > 0, na.rm = TRUE),
                     n_blocks = nlevels(blk), block_bp = BOOT_BLOCK, reps = BOOT_N, valid_reps = sum(!is.na(D_boot)),
                     reading = "D > 0: the labelled pairs (T1-N1, T2-N2) are more discordant than the swapped ones")
boot_support <- data.table(assignment = colnames(ASSIGN_B))[, support := vapply(assignment, function(a) mean(win_b == a), 0)]
boot_support <- assign_scores[combo == PRIMARY & method == "disc_nhom", .(assignment, total_score, rank, is_labelled, is_hypothesis)][boot_support, on = "assignment"]
setorder(boot_support, rank)
print(boot_D); print(boot_support)
WINNER <- boot_support[rank == 1, assignment]
assign_scores[, clearly_separated := NA]
assign_scores[combo == PRIMARY & method == "disc_nhom" & rank == 1, clearly_separated := boot_support[rank == 1, support] >= 0.95]
assign_scores[combo == PRIMARY & method == "disc_nhom", bootstrap_support := boot_support[.SD, on = "assignment", x.support]]
assign_scores[rank <= 2 & combo == PRIMARY]

# Verdict on the primary combo (panel A, THIN_PRIMARY), primary metric = normal-hom discordance
best_of <- function(d, tum = TUMOUR_IDS, nrms = NORMAL_IDS) {
  rbindlist(lapply(tum, function(t) {
    r <- rbindlist(lapply(nrms, function(n) { pc <- pair_counts(d, t, n)
      data.table(normal = n, n_inf = pc$n_inf, disc = pc$disc, rate = pc$disc / pc$n_inf, ibs0_rate = pc$ibs0 / pc$n_joint,
                 ibs2_rate = pc$ibs2 / pc$n_joint) }))
    setorder(r, rate, ibs0_rate, -ibs2_rate)
    pp <- paired_counts(d, t, r$normal[1], r$normal[2])
    data.table(tumour = t, best_normal = r$normal[1], second_normal = r$normal[2], best_rate = r$rate[1], second_rate = r$rate[2],
               n_inf_best = r$n_inf[1], paired_p = pp$p_value,
               ambiguous = pp$p_value >= PAIRED_ALPHA | r$n_inf[1] < MIN_JOINT_SITES)
  }))
}
# [CHANGE 26] discordant segments: consecutive discordant informative sites on one chromosome with no
# concordant informative site between them
segments <- rbindlist(lapply(seq_len(nrow(TN_PAIRS)), function(i) {
  t <- TN_PAIRS$tumour[i]; n <- TN_PAIRS$normal[i]
  v <- tn_site(prim_d, t, n)
  x <- data.table(chr_n = prim_d$chr_n, disc = v$disc)[v$inf]
  seg <- x[, .(n_seg = { r <- rle(disc); sum(r$values) }, any_disc = any(disc)), by = chr_n]
  data.table(tumour = t, normal = n, informative = nrow(x), discordant = sum(x$disc), n_segments = sum(seg$n_seg),
             n_chromosomes = sum(seg$any_disc))
}))
seg_txt <- function(t, n) { r <- segments[tumour == t & normal == n]; sprintf("%d segments on %d chromosomes", r$n_segments, r$n_chromosomes) }

vr <- robust[combo == PRIMARY & method == "disc_nhom"]
verdict <- vr[, .(analysis = "primary", tumour, labelled_normal, best_normal, second_normal,
                  best_disc_rate = best_score, best_disc_ci_hi = best_hi, second_disc_rate = second_score, second_disc_ci_lo = second_lo,
                  n_inf_best = n_eff_best, paired_p_best_vs_second, ambiguous, ci_overlap,
                  consistent_with_label = matches_labelled, matches_hypothesis)]
lab_rate <- tn[combo == PRIMARY & method == "disc_nhom" & normal == KI_NORMALS[tumour],
               .(tumour, labelled_disc_rate = disc_rate, labelled_n_inf = n_inf, labelled_ibs0_rate = ibs0_rate)]
verdict <- lab_rate[verdict, on = "tumour"]
verdict[, evidence_vs_second := mapply(seg_txt, tumour, second_normal)]
verdict[, evidence_vs_labelled := fifelse(best_normal == labelled_normal, "", mapply(seg_txt, tumour, labelled_normal))]
for (m in c("custom", "somalier", "gtcheck")) {
  r <- robust[combo == PRIMARY & method == m, .(tumour, b = best_normal, a = ambiguous)]
  lab <- if (m == "custom") "ibs0" else m
  if (nrow(r)) verdict[r, on = "tumour", c(paste0(lab, "_best"), paste0(lab, "_ambiguous")) := .(i.b, i.a)]
}

# [CHANGE 30] ascertainment sensitivity: sites required only in the files of T1 and T2 (file 3 not
# required), their labelled normals confident; T1/T2 against N1/N2 only
asc <- data.table()
if (length(TUMOUR_IDS) >= 2) {
  tum12 <- TUMOUR_IDS[1:2]
  sel12 <- select_sites(MIN_DP, tum12, log = TRUE, scope = "files 1+2 only")
  if (nrow(sel12)) {
    G12 <- build_G(sel12, tum12)$G
    k12 <- thin_greedy(G12, THIN_PRIMARY)
    b12 <- best_of(G12[key %chin% k12], tum12, unname(KI_NORMALS[tum12]))
    asc <- b12[, .(analysis = sprintf("ascertainment: files 1+2 only (%d sites, %s thinning)", length(k12), thin_label(THIN_PRIMARY)),
                   tumour, labelled_normal = KI_NORMALS[tumour], best_normal, second_normal, best_disc_rate = best_rate,
                   second_disc_rate = second_rate, n_inf_best, paired_p_best_vs_second = paired_p, ambiguous,
                   consistent_with_label = best_normal == KI_NORMALS[tumour], matches_hypothesis = best_normal == HYPOTHESIS[tumour])]
  } else deviation("ascertainment sensitivity (files 1+2): no sites")
}
verdict <- rbind(verdict, asc, fill = TRUE)
setorder(verdict, analysis, best_disc_rate)
verdict[, analysis := factor(analysis, unique(c("primary", analysis)))]; setorder(verdict, analysis, best_disc_rate)

# [CHANGE 25] leave-one-chromosome-out jackknife on the primary combination
full_best <- setNames(verdict[analysis == "primary", best_normal], verdict[analysis == "primary", tumour])
full_amb  <- setNames(verdict[analysis == "primary", ambiguous], verdict[analysis == "primary", tumour])
total_of <- function(d, map) sum(vapply(names(map), function(t) { pc <- pair_counts(d, t, map[[t]]); pc$disc / pc$n_inf }, 0))
jackknife <- rbindlist(lapply(sort(unique(prim_d$chr_n)), function(ch) {
  d <- prim_d[chr_n != ch]
  b <- best_of(d)
  b[, `:=`(dropped_chr = ch, n_sites = nrow(d), changed = best_normal != full_best[tumour],
           ambiguity_flipped = ambiguous != full_amb[tumour],             # [CHANGE 37]
           gap_labelled_minus_hypothesis = total_of(d, KI_NORMALS) - total_of(d, HYPOTHESIS))]
}))
setcolorder(jackknife, c("dropped_chr", "n_sites", "tumour"))
jk_flag <- jackknife[changed == TRUE, sort(unique(dropped_chr))]
if (length(jk_flag)) deviation("jackknife: the verdict changes when chromosome(s) ", paste(jk_flag, collapse = ", "), " are dropped (see jackknife.tsv)")
# [CHANGE 37] ambiguity flips
jk_amb <- jackknife[ambiguity_flipped == TRUE, .(chrs = paste(sort(dropped_chr), collapse = ",")), by = tumour]
for (i in seq_len(nrow(jk_amb))) deviation("jackknife: ", jk_amb$tumour[i], " changes its ambiguity flag (",
                                           if (full_amb[[jk_amb$tumour[i]]]) "ambiguous -> clear" else "clear -> ambiguous",
                                           ") when chromosome(s) ", jk_amb$chrs[i], " are dropped (see jackknife.tsv)")
jackknife[, .(changed = sum(changed), ambiguity_flipped = sum(ambiguity_flipped), max_paired_p = max(paired_p),
              min_gap = min(gap_labelled_minus_hypothesis), max_gap = max(gap_labelled_minus_hypothesis)), by = tumour]

# [CHANGE 31] MIN_DP sensitivity of the primary verdict
mindp_sens <- rbindlist(lapply(c(MIN_DP, MIN_DP_SENS), function(md) {
  sel <- select_sites(md)
  if (!nrow(sel)) { deviation("MIN_DP sensitivity: no sites at MIN_DP = ", md); return(NULL) }
  Gm <- build_G(sel, TUMOUR_IDS, md)$G
  km <- thin_greedy(Gm, THIN_PRIMARY)
  best_of(Gm[key %chin% km])[, `:=`(min_dp = md, n_sites = length(km), primary_setting = md == MIN_DP)]
}))
mindp_sens[, same_as_primary := best_normal == full_best[tumour]]
setcolorder(mindp_sens, c("min_dp", "primary_setting", "n_sites", "tumour"))
mindp_sens

# [CHANGE 36] FreeBayes PASS exclusion as a sensitivity analysis (primary panel and thinning)
pass_sens <- rbindlist(lapply(c(EXCLUDE_FB_PASS, !EXCLUDE_FB_PASS), function(ex) {
  sel <- select_sites(MIN_DP, exclude_pass = ex)
  if (!nrow(sel)) { deviation("PASS sensitivity: no sites with EXCLUDE_FB_PASS = ", ex); return(NULL) }
  Gp <- build_G(sel)$G
  kp <- thin_greedy(Gp, THIN_PRIMARY)
  best_of(Gp[key %chin% kp])[, `:=`(exclude_fb_pass = ex, primary_setting = ex == EXCLUDE_FB_PASS,
                                    sites_unthinned = nrow(sel), n_sites = length(kp))]
}))
pass_sens[, same_as_primary := best_normal == full_best[tumour] & ambiguous == full_amb[tumour]]
setcolorder(pass_sens, c("exclude_fb_pass", "primary_setting", "sites_unthinned", "n_sites", "tumour"))
pass_sens

# [CHANGE 26] informative and discordant sites near the knock-in (+- KI_FLANK)
ki_locus <- if (!is.na(KNOCKIN_LOCUS$chrom)) rbindlist(lapply(seq_len(nrow(TN_PAIRS)), function(i) {
  v <- tn_site(prim_d, TN_PAIRS$tumour[i], TN_PAIRS$normal[i])
  near <- prim_d$chr_n == as.integer(sub("^chr", "", KNOCKIN_LOCUS$chrom)) & abs(prim_d$POS - KNOCKIN_LOCUS$pos) <= KI_FLANK
  data.table(tumour = TN_PAIRS$tumour[i], normal = TN_PAIRS$normal[i], informative = sum(v$inf), informative_near_ki = sum(v$inf & near),
             discordant = sum(v$disc), discordant_near_ki = sum(v$disc & near))[
    , `:=`(frac_informative_near_ki = round(informative_near_ki / informative, 3),
           frac_discordant_near_ki = round(fifelse(discordant > 0, discordant_near_ki / discordant, NA_real_), 3))]
})) else data.table()

print(verdict)

# Plots ------------------------------------------------------------------------------
supported <- verdict[analysis == "primary" & ambiguous == FALSE, .(tumour, normal = best_normal)]
heat_metric <- function(meth, col, label, digits = 3, lower_better = TRUE, ncol = "n_joint", subtitle = NULL) {
  d <- pairwise[combo == PRIMARY & method == meth & !is.na(get(col)), .(s1, s2, v = get(col), n = get(ncol))]
  if (!nrow(d)) return(invisible(NULL))
  d <- rbind(d, d[, .(s1 = s2, s2 = s1, v, n)])
  d[, `:=`(s1 = factor(s1, SAMPLES), s2 = factor(s2, rev(SAMPLES)))]
  ol <- rbind(data.table(tumour = names(KI_NORMALS), normal = unname(KI_NORMALS), outline = "labelled pair"),
              supported[, .(tumour, normal, outline = "genetically supported")])
  ol <- rbind(ol[, .(s1 = tumour, s2 = normal, outline)], ol[, .(s1 = normal, s2 = tumour, outline)])
  ol[, `:=`(s1 = factor(s1, SAMPLES), s2 = factor(s2, rev(SAMPLES)))]
  p <- ggplot(d, aes(s1, s2)) + geom_tile(aes(fill = v), colour = "white") +
    geom_text(aes(label = sprintf(paste0("%.", digits, "f\nn=%d"), v, as.integer(n))), size = 3) +
    geom_tile(data = ol, aes(colour = outline, linetype = outline), fill = NA, linewidth = 1.1) +
    scale_colour_manual(values = c(`labelled pair` = "black", `genetically supported` = "#C62828"), name = NULL) +
    scale_linetype_manual(values = c(`labelled pair` = "dashed", `genetically supported` = "solid"), name = NULL) +
    scale_fill_gradient(low = if (lower_better) "#08519C" else "white", high = if (lower_better) "white" else "#08519C", name = label) +
    labs(x = NULL, y = NULL, title = paste0(label, " - ", meth, ", ", PRIMARY),
         subtitle = if (is.null(subtitle)) "n = sites compared for that metric. Outlines: labelled pairing (dashed) vs best primary-metric match (red)" else subtitle) +
    coord_equal() + theme_bw()
  save_plot(p, paste0("heatmap_", meth, "_", col), 8, 7)
}
heat_metric("custom", "disc_rate", "discordance at hom sites", 3, ncol = "n_inf",
            subtitle = "Tumour-normal: normal homozygous, tumour != normal (primary). Normal-normal / tumour-tumour: either homozygous, genotypes differ. n = informative sites")
heat_metric("custom", "ibs0_rate", "IBS0 / joint", 4)
heat_metric("custom", "ibs2_rate", "IBS2 / joint", 3, lower_better = FALSE)
heat_metric("custom", "hom_concordance", "homozygous concordance", 3, lower_better = FALSE)
heat_metric("somalier", "somalier_concordance", "somalier concordance", 2, lower_better = FALSE)   # [CHANGE 33] relatedness heatmap removed
heat_metric("gtcheck", "gtcheck_discordance", "gtcheck discordance", 3)

# Primary metric per tumour vs each KI normal with 95% CI [CHANGE 22]
pb <- tn[combo == PRIMARY & method == "disc_nhom"]
pb[, label_type := fcase(normal == KI_NORMALS[tumour], "labelled normal", normal == HYPOTHESIS[tumour], "hypothesis normal", default = "other normal")]
p_disc <- ggplot(pb, aes(normal, disc_rate, fill = label_type)) +
  geom_col(width = 0.7) + geom_errorbar(aes(ymin = disc_lo, ymax = disc_hi), width = 0.2) +
  geom_text(aes(y = disc_hi, label = sprintf("%d / %d", as.integer(disc), as.integer(n_inf))), vjust = -0.4, size = 3.2) +
  facet_wrap(~ tumour, labeller = labeller(tumour = function(x) paste("Tumour", x))) +
  scale_fill_manual(values = c(`labelled normal` = "#C62828", `hypothesis normal` = "#1565C0", `other normal` = "grey60"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "KI normal liver", y = "discordant / normal-homozygous sites (95% CI)",
       title = "Which KI normal matches each tumour? Lowest discordance at normal-homozygous sites = same animal",
       subtitle = sprintf("Panel A, thinned %s: %s sites. Masked tumour genotypes excluded.", thin_label(THIN_PRIMARY), format(combos[primary == TRUE, n_sites], big.mark = ","))) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
save_plot(p_disc, "disc_per_tumour", 10, 5)

# IBS0 rate per tumour vs each KI normal (secondary) [CHANGE 18]
pi0 <- tn[combo == PRIMARY & method == "custom"]
pi0[, label_type := fcase(normal == KI_NORMALS[tumour], "labelled normal", normal == HYPOTHESIS[tumour], "hypothesis normal", default = "other normal")]
p_ibs0 <- ggplot(pi0, aes(normal, ibs0_rate, fill = label_type)) +
  geom_col(width = 0.7) + geom_errorbar(aes(ymin = ibs0_lo, ymax = ibs0_hi), width = 0.2) +
  geom_text(aes(y = ibs0_hi, label = sprintf("%d / %d", as.integer(ibs0), as.integer(n_joint))), vjust = -0.4, size = 3.2) +
  facet_wrap(~ tumour, labeller = labeller(tumour = function(x) paste("Tumour", x))) +
  scale_fill_manual(values = c(`labelled normal` = "#C62828", `hypothesis normal` = "#1565C0", `other normal` = "grey60"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "KI normal liver", y = "IBS0 / jointly callable sites (95% CI)",
       title = "IBS0 per tumour (secondary): detects a swap, cannot confirm correct labels under this ascertainment",
       subtitle = sprintf("Panel A, thinned %s", thin_label(THIN_PRIMARY))) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
save_plot(p_ibs0, "ibs0_per_tumour", 10, 5)

# IBS0 vs IBS2, all pairs [CHANGE 18]
sc <- pairwise[combo == PRIMARY & method == "custom"]
sc[, pair := paste(s1, "-", s2)]
sc[, colour_type := fcase(labelled_pair, "labelled pair", type == "tumour-normal", "other tumour-normal", default = type)]
lab_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) ggrepel::geom_text_repel(aes(label = pair), size = 3) else
  geom_text(aes(label = pair), size = 3, vjust = -0.8)
p_ibs <- ggplot(sc, aes(ibs2_rate, ibs0_rate, colour = colour_type)) + geom_point(size = 3) + lab_layer +
  scale_colour_manual(values = c(`labelled pair` = "#C62828", `other tumour-normal` = "grey30", `normal-normal` = "#1565C0", `tumour-tumour` = "#2E7D32"), name = NULL) +
  labs(x = "IBS2 / joint (identical genotypes)", y = "IBS0 / joint (opposite homozygotes)",
       title = paste("IBS0 vs IBS2, all sample pairs -", PRIMARY),
       subtitle = "Same animal: bottom right. Normal-normal and tumour-tumour = different-animal baseline") + theme_bw()
save_plot(p_ibs, "ibs0_vs_ibs2", 9, 7)

# Robustness: best normal per tumour across panels x thinning x methods [CHANGE 18]
rb <- copy(robust)
rb[, combo := factor(combo, combos$combo)]
p_rob <- ggplot(rb, aes(combo, tumour, fill = best_normal)) + geom_tile(colour = "white") +
  geom_text(aes(label = paste0(best_normal, fifelse(ambiguous, "?", ""))), size = 3) +
  facet_wrap(~ method, ncol = 1) +
  labs(x = "panel_thinning", y = NULL, fill = "best normal",
       title = "Best-matching KI normal per tumour, every panel x thinning x method (disc_nhom = primary)",
       subtitle = "? = ambiguous (disc_nhom: paired sign test vs second best p >= 0.05; others: CI overlap; or too few sites / underpowered). Panels B/C are defined from normal differences.") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_rob, "robustness", 12, 2 + 1.8 * uniqueN(rb$method))

# Bootstrap distribution of D [CHANGE 24]
p_boot <- ggplot(data.table(D = D_boot), aes(D)) + geom_histogram(bins = 50, fill = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed") + geom_vline(xintercept = D_point, colour = "#C62828") +
  labs(x = sprintf("D = [s(%s,%s) + s(%s,%s)] - [s(%s,%s) + s(%s,%s)]", T1, N1, T2, N2, T1, N2, T2, N1), y = "bootstrap replicates",
       title = sprintf("%s-block bootstrap of the labelled-vs-swapped difference (primary metric)", thin_label(BOOT_BLOCK)),
       subtitle = "D > 0: labelled pairs more discordant than swapped pairs; red = point estimate") + theme_bw()
save_plot(p_boot, "bootstrap_D", 8, 5)

# Chromosomal distribution: panel B/C sites, and IBS0 sites per tumour-normal comparison [CHANGE 18]
bc <- rbind(G[panel_B == TRUE, .(chr_n, POS, set = "panel B")], G[panel_C == TRUE, .(chr_n, POS, set = "panel C")])
ibs0_sites <- rbindlist(lapply(TUMOUR_IDS, function(t) rbindlist(lapply(NORMAL_IDS, function(n) {
  prim_d[!is.na(get(t)) & abs(get(t) - get(n)) == 2, .(chr_n, POS, comparison = paste(t, "-", n), labelled = KI_NORMALS[t] == n)]
}))))
ki_layer <- if (!is.na(KNOCKIN_LOCUS$chrom))
  geom_vline(data = data.table(chr_n = as.integer(sub("^chr", "", KNOCKIN_LOCUS$chrom)), x = KNOCKIN_LOCUS$pos / 1e6),
             aes(xintercept = x), colour = "#C62828", linetype = "dashed") else NULL
p_chr <- ggplot() +
  geom_point(data = bc, aes(POS / 1e6, set), shape = "|", size = 3, alpha = 0.6) +
  geom_point(data = ibs0_sites, aes(POS / 1e6, comparison, colour = labelled), shape = "|", size = 3) +
  ki_layer + facet_wrap(~ chr_n, scales = "free_x", ncol = 5) +
  scale_colour_manual(values = c(`TRUE` = "#C62828", `FALSE` = "grey20"), labels = c(`TRUE` = "labelled pair", `FALSE` = "other pair"), name = "IBS0 sites") +
  labs(x = "position (Mb)", y = NULL, title = paste("Where the evidence lies: panel B/C sites and IBS0 sites (", PRIMARY, ")"),
       subtitle = paste0("Few blocks = evidence from few inherited haplotypes; thinning does not make sites independent.",
                         if (!is.na(KNOCKIN_LOCUS$chrom)) " Dashed line: knock-in locus." else "")) +
  theme_bw(base_size = 9)
save_plot(p_chr, "chromosomal_distribution", 16, 12)
knockin_tab <- if (!is.na(KNOCKIN_LOCUS$chrom) && nrow(ibs0_sites)) ibs0_sites[, .(ibs0_sites = .N,
  near_knockin = sum(chr_n == as.integer(sub("^chr", "", KNOCKIN_LOCUS$chrom)) & abs(POS - KNOCKIN_LOCUS$pos) <= KNOCKIN_FLANK)),
  by = .(comparison, labelled)][, frac_near_knockin := round(near_knockin / ibs0_sites, 3)] else data.table()
ibs0_by_chr <- ibs0_sites[, .N, by = .(comparison, chr_n)][order(comparison, chr_n)]

# ---- 6. Outputs, RUN_LOG.md ------------------------------------------------------
fwrite(FILTER_COUNTS, file.path(OUT, "filter_counts.tsv"), sep = "\t")
fwrite(ad_check, file.path(OUT, "ad_check.tsv"), sep = "\t")
fwrite(masked_summary, file.path(OUT, "masked_genotypes_per_sample.tsv"), sep = "\t")
fwrite(combos[, .(combo, panel, thin, panel_unthinned, n_sites, underpowered, primary)], file.path(OUT, "panel_sizes.tsv"), sep = "\t")
fwrite(pairwise, file.path(OUT, "pairwise_metrics.tsv"), sep = "\t")
fwrite(assign_scores, file.path(OUT, "assignment_scores.tsv"), sep = "\t")
fwrite(verdict, file.path(OUT, "verdict.tsv"), sep = "\t")
fwrite(robust, file.path(OUT, "robustness.tsv"), sep = "\t")
fwrite(paired, file.path(OUT, "paired_tests.tsv"), sep = "\t")                       # [CHANGE 23]
fwrite(boot_D, file.path(OUT, "bootstrap_D.tsv"), sep = "\t")                        # [CHANGE 24]
fwrite(boot_support, file.path(OUT, "bootstrap_assignment_support.tsv"), sep = "\t") # [CHANGE 24]
fwrite(jackknife, file.path(OUT, "jackknife.tsv"), sep = "\t")                       # [CHANGE 25]
fwrite(segments, file.path(OUT, "segments.tsv"), sep = "\t")                         # [CHANGE 26]
if (nrow(ki_locus)) fwrite(ki_locus, file.path(OUT, "ki_locus.tsv"), sep = "\t")     # [CHANGE 26]
fwrite(mindp_sens, file.path(OUT, "mindp_sensitivity.tsv"), sep = "\t")              # [CHANGE 31]
fwrite(pass_sens, file.path(OUT, "pass_sensitivity.tsv"), sep = "\t")                # [CHANGE 36]
if (nrow(som_n_check)) fwrite(som_n_check, file.path(OUT, "somalier_n_check.tsv"), sep = "\t")   # [CHANGE 33]
fwrite(ibs0_by_chr, file.path(OUT, "ibs0_sites_per_chromosome.tsv"), sep = "\t")
if (nrow(knockin_tab)) fwrite(knockin_tab, file.path(OUT, "ibs0_near_knockin.tsv"), sep = "\t")
fwrite(G[key %chin% combo_sites[[PRIMARY]], c("CHROM", "POS", "REF", "ALT", "min_dp", "dev", "AF", SAMPLES), with = FALSE],
       file.path(OUT, "sites_used_primary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
git_hash <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(dirname(SCRIPT_RUN)), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)),
                     error = function(e) "unknown")
if (!length(git_hash)) git_hash <- "unknown"
git_dirty <- tryCatch(length(suppressWarnings(system2("git", c("-C", shQuote(dirname(SCRIPT_RUN)), "status", "--porcelain", shQuote(SCRIPT_RUN)),
                                                      stdout = TRUE, stderr = FALSE))) > 0, error = function(e) NA)

# [CHANGE 20] [CHANGE 35] RUN_LOG.md: changelog (copied from the header), assumptions, limitations, deviations
hdr_lines <- readLines(SCRIPT_RUN)
changelog <- hdr_lines[seq(grep("^# CHANGELOG", hdr_lines)[1], grep("^# END CHANGELOG", hdr_lines)[1])]
md_table <- function(x) { x <- as.data.table(x); if (!nrow(x)) return("(none)")
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))) }
fmt_num <- function(d) { d <- copy(d); for (cl in names(d)) if (is.double(d[[cl]])) set(d, j = cl, value = signif(d[[cl]], 4)); d }
log_entry <- c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_somalier", if (TEST_MODE) paste0(" (TEST_MODE, synthetic data, scenario ", SCENARIO, ")")),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Tools:      ", somalier_version, "; ", paste(bcftools_version, collapse = " / ")),
  paste0("- Tool help:  ", if (length(help_files)) paste(help_files, collapse = ", ") else "(not saved: tools missing)"),
  paste0("- Script:     ", SCRIPT_RUN, "  git: ", git_hash, if (isTRUE(git_dirty)) " (+ uncommitted changes)" else ""),
  sprintf("- Parameters: MIN_QUAL=%s MIN_DP=%s (sensitivity %s) bins hom_ref<=%s het %s-%s hom_alt>=%s; MAX_DP_FACTOR=%s MASK_BED=%s EXCLUDE_FB_PASS=%s; thinning %s (primary) + %s, THIN_SEED=%s; MIN_PANEL_SITES=%s MIN_JOINT_SITES=%s PAIRED_ALPHA=%s; bootstrap %s blocks x %s (seed %s); KNOCKIN=%s (KI_FLANK %s)",
          MIN_QUAL, MIN_DP, paste(MIN_DP_SENS, collapse = ","), HOM_REF_MAX, HET_MIN, HET_MAX, HOM_ALT_MIN, MAX_DP_FACTOR, MASK_BED, EXCLUDE_FB_PASS,
          thin_label(THIN_PRIMARY), paste(thin_label(THIN_SENS), collapse = ", "), THIN_SEED, MIN_PANEL_SITES, MIN_JOINT_SITES, PAIRED_ALPHA,
          thin_label(BOOT_BLOCK), BOOT_N, BOOT_SEED,
          if (is.na(KNOCKIN_LOCUS$chrom)) "NA" else paste0(KNOCKIN_LOCUS$chrom, ":", format(KNOCKIN_LOCUS$pos, scientific = FALSE)),
          format(KI_FLANK, scientific = FALSE)),
  paste0("- Labelled:   ", paste(names(KI_NORMALS), KI_NORMALS, sep = "-", collapse = ", "), "; hypothesis tested: ", paste(names(HYPOTHESIS), HYPOTHESIS, sep = "-", collapse = ", ")),
  paste0("- Inputs:     ", paste(basename(fb_files), collapse = ", "), "; FASTA ", FASTA),
  paste0("- Outputs:    ", OUT),
  "", "### Changelog of the script", "", sub("^# ?", "", changelog),
  "", "### Filter counts", "", md_table(FILTER_COUNTS),
  "", "### AD check (AD vs RO,AO)", "", md_table(ad_check),
  "", "### Masked genotypes per sample (final sites)", "", md_table(masked_summary),
  "", "### Panel sizes", "", md_table(combos[, .(combo, panel_unthinned, n_sites, underpowered, primary)]),
  "", paste0("### Verdict (primary: panel A, thinned ", thin_label(THIN_PRIMARY), ", discordance at normal-homozygous sites; plus the files-1+2 ascertainment rows)"), "",
  md_table(fmt_num(verdict)),
  "", "### Paired sign tests, primary combination", "", md_table(fmt_num(paired[combo == PRIMARY])),
  "", "### One-to-one assignments, primary, normal-hom discordance", "",
  md_table(fmt_num(assign_scores[combo == PRIMARY & method == "disc_nhom", .(rank, assignment, total_score, gap_to_next, bootstrap_support, clearly_separated, separated_z2, is_labelled, is_hypothesis)])),
  "", "### One-to-one assignments, primary, IBS0 (secondary)", "",
  md_table(fmt_num(assign_scores[combo == PRIMARY & method == "custom", .(rank, assignment, total_score, gap_to_next, separated_z2, is_labelled, is_hypothesis)])),
  "", "### Block bootstrap", "", md_table(fmt_num(boot_D)), "", md_table(fmt_num(boot_support)),
  "", "### Jackknife (leave one chromosome out)", "",
  if (length(jk_flag)) paste0("Verdict changes when dropping chromosome(s): ", paste(jk_flag, collapse = ", ")) else "The verdict never changes when one chromosome is dropped.",
  if (nrow(jk_amb)) paste0("Ambiguity flag flips: ", paste(jk_amb$tumour, "(chr", jk_amb$chrs, ")", collapse = "; ")) else "No tumour changes its ambiguity flag when one chromosome is dropped.",
  "", md_table(fmt_num(jackknife[, .(changed = sum(changed), ambiguity_flipped = sum(ambiguity_flipped), max_paired_p = max(paired_p),
                                     min_gap_labelled_minus_hypothesis = min(gap_labelled_minus_hypothesis),
                                     max_gap_labelled_minus_hypothesis = max(gap_labelled_minus_hypothesis)), by = tumour])),
  "", "### Discordant segments (primary combination)", "", md_table(segments),
  if (nrow(ki_locus)) c("", paste0("### Sites near the knock-in (+-", format(KI_FLANK, scientific = FALSE), " bp)"), "", md_table(ki_locus)) else character(),
  "", "### MIN_DP sensitivity", "", md_table(fmt_num(mindp_sens)),
  "", "### FreeBayes PASS exclusion sensitivity", "",
  paste0("Primary run: EXCLUDE_FB_PASS = ", EXCLUDE_FB_PASS, "; the other setting is the sensitivity run (", length(pass_keys), " PASS keys in the files)."),
  "", md_table(fmt_num(pass_sens)),
  "", "### Robustness (best normal; ? = ambiguous)", "", md_table(dcast(robust, method + tumour ~ combo, value.var = "call")),
  "", "### Denominators and definitions",
  "- Primary metric (disc_nhom), tumour-normal: informative = sites where the normal is homozygous (0/0 or 1/1) and the tumour is unmasked; discordant = tumour genotype != normal genotype; rate with exact 95% CI. LOH only turns het into hom, so a tumour het (or the opposite homozygote) against a homozygous normal is mismatch evidence; it stays observable when every animal is forced non-ref by the ascertainment (het vs 1/1).",
  "- Normal-normal and tumour-tumour pairs: symmetric version (both unmasked, at least one homozygous, genotypes differ). Reference values only.",
  "- Paired sign test: for tumour T and normals Nj, Nk, sites where T, Nj, Nk are unmasked and at least one normal homozygous; T is discordant with a normal only where that normal is homozygous; b = discordant with Nj only, c = with Nk only; exact two-sided binom.test(b, b + c). Ambiguous (primary) = p >= PAIRED_ALPHA vs the second best, or < MIN_JOINT_SITES informative sites, or underpowered combination. ci_overlap = the old CI-overlap flag.",
  "- DEVIATION FROM THE REVIEW SPEC (paired test): the spec restricted the test to sites where both normals are homozygous; that version is reported as *_bothhom columns in paired_tests.tsv. It drops het-vs-hom sites, which carry most differences between animals of an inbred-like colony, and was underpowered in --test (p = 0.25 for 0/1127 vs 80/938 discordant sites).",
  paste0("- Block bootstrap: ", thin_label(BOOT_BLOCK), " blocks (chromosome x floor(POS / block)) resampled with replacement, ", BOOT_N, " replicates; rates recomputed from the resampled counts. D and the support of each one-to-one assignment (fraction of replicates it has the lowest total). clearly_separated = support of the winning assignment >= 0.95; separated_z2 = the old rule (gap > 2 x sqrt(summed binomial variances))."),
  "- Jackknife: the primary verdict and gap = total(labelled) - total(hypothesis) recomputed with each autosome removed.",
  "- Segments: per tumour-normal pair, informative sites in genome order; a segment = a run of discordant informative sites on one chromosome with no concordant informative site between them.",
  "- custom IBS0 (secondary): n_joint = both unmasked; IBS0 = 0/0 vs 1/1; IBS2 = identical genotype; homozygous concordance = identical among both-homozygous sites. Exact 95% CIs.",
  "- somalier: ibs0/n and ibs2/n with somalier's n = sites known in both (its genotypes, from AD). relatedness_exploratory = 2*(shared_hets - 2*ibs0)/het_ab (v0.3.5 source); `concordance` is somalier's adjusted homozygous concordance. somalier n must equal the custom n_joint (somalier_n_check.tsv).",
  paste0("- somalier genotype calls (v0.3.5 source, relate.nim `alts`): unknown if REF+ALT < -d; hom-ref if ALT = 0 or AB < cutoff; hom-alt if REF = 0 or AB > 1 - cutoff; het if min-ab <= AB <= 1 - min-ab; else unknown. Run with -d ", MIN_DP, ", --min-ab ", HET_MIN, ", SOMALIER_AB_HOM_CUTOFF=", AB_HOM_CUTOFF, "; masked genotypes have AD '.', read as depth 0 -> unknown. -u is never used."),
  paste0("- gtcheck: ", if (is.null(GT_ARGS)) "not run" else paste("bcftools gtcheck", paste(GT_ARGS, collapse = " "), "-T targets masked file (no -g: all sample pairs)"),
         "; discordance per site compared (any genotype difference, so tumour LOH raises it)."),
  "- MIN_DP sensitivity: whole selection rerun at each depth, primary panel and thinning. At depth 20 the probability that a true het lands in a hom bin (<= 1 opposite read) is ~2e-5 (21 / 2^20).",
  "", "### Assumptions",
  "- Each FreeBayes file was called on its own labelled pair; a site missing from a file is unknown, never hom-ref.",
  "- VAF bins give an inferred genotype under a diploid model; FreeBayes' own GT is not used.",
  "- INFO/AF in the sites and masked files is computed in-sample from the 3 KI normals; somalier uses it only for its per-sample hom counts and contamination estimate, not for the pair statistics used here.",
  "- KI normals define the sites; tumours are only required to be callable (depth) and are never used to choose sites.",
  "- FreeBayes PASS records are the provider's somatic prefilter (normal 0/0, tumour 0/1, ...). They are KEPT by default (EXCLUDE_FB_PASS = FALSE): in a swapped pair the prefilter selects germline differences between the two mice, so removing PASS sites deletes swap evidence and makes the site set depend on the true labelling (--test: 2,699 PASS keys in the swap vs 462 in the null scenario). Tumour-private somatic sites are excluded anyway, because a KI normal must be confidently non-ref. The opposite setting is reported in pass_sensitivity.tsv.",
  "", "### Interpretation and limitations",
  "- Ascertainment asymmetry: a site is in a pair's file only if one of the pair's samples is non-ref. Requiring all three files forces every correctly paired animal to be non-ref (het or hom-alt) at every retained site: 0/0 vs 1/1 (IBS0) cannot occur between correctly labelled animals, so IBS0 can detect a swap but cannot confirm correct labels. The primary metric (normal-hom discordance) still sees het vs 1/1 differences. The files-1+2 rows relax the requirement to the two tumours under test.",
  "- Linkage: mouse haplotype blocks span tens of Mb; sites are not independent and site-level CIs are overconfident. Use the paired tests together with the jackknife, the block bootstrap and the segment counts.",
  "- Genetics cannot tell whether the tumour label or the normal label is wrong, or where the swap happened (tissue, tube, library, file naming); only records can.",
  "- There are no human relatedness thresholds for inbred, genetically modified mice; genetic concordance shows which samples share inherited genotypes, it is not definitive proof of individual identity.",
  "- All methods use the same FreeBayes calls and the same masked file: they are complementary, not independent evidence.",
  "- Tumours: purity, LOH and copy number make tumour hets unreliable; ambiguous tumour genotypes are masked. LOH cannot create discordance at normal-homozygous sites.",
  "- Collapsed duplications look like germline hets in all samples; the max-depth filter and the all-normals-het count address this. The max-depth filter does not catch paralog artefacts at normal depth (e.g. reads of a paralog mapped at typical coverage); only a MASK_BED of low-mappability / segmental-duplication regions does. Without a MASK_BED these remain.",
  "- Panels B and C are defined from differences between the KI normals, so normal-normal metrics on them are inflated by construction; only tumour-normal comparisons are interpretable there.",
  "- Panels B and C are expected to be underpowered when most animals are forced non-ref (correct labels): with every animal non-ref at the retained sites, KI normals rarely differ and never as 0/0 vs 1/1 (panel C is near empty in the --test null scenario).",
  "- If results are ambiguous or panels underpowered: genotype the panel A sites from the BAMs in all six samples (e.g. bcftools mpileup -T sites) once BAMs are available.",
  "", "### Deviations", "", if (length(DEVIATIONS)) paste0("- DEVIATION: ", DEVIATIONS) else "- none", "")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(log_entry, old_log), log_file)

# [CHANGE 1] [CHANGE 34] TEST_MODE checks
if (TEST_MODE) {
  vp <- verdict[analysis == "primary"]
  got <- setNames(vp$best_normal, vp$tumour)
  truth_map <- if (SCENARIO == "swap") c(HCC_1 = "KI2", HCC_2 = "KI1", HCC_3 = "KI3") else KI_NORMALS
  tsites <- fread(file.path(DIR_A, "truth_sites.tsv"))
  fc <- FILTER_COUNTS[scope == "all files"]
  # artefact sites reaching the depth step vs removed by the max-depth step
  sel_dbg <- local({ f <- full[sample %chin% SAMPLES]; f[, med := MED[sample]]
    f[, .(n_present = sum(status != "absent"), n_dp = sum(!is.na(DPA) & DPA >= MIN_DP),
          too_deep = any(!is.na(DPA) & DPA > MAX_DP_FACTOR * med)), by = key] })
  art <- sel_dbg[key %chin% tsites[artefact == TRUE, key] & n_present == 6 & n_dp == 6]
  checks <- c(
    primary_matches_truth    = all(got[names(truth_map)] == truth_map),
    primary_not_ambiguous    = !any(vp$ambiguous),
    winning_assignment_truth = assign_scores[combo == PRIMARY & method == "disc_nhom" & rank == 1,
                                             if (SCENARIO == "swap") is_hypothesis else is_labelled],
    ibs0_swap_detected       = if (SCENARIO == "swap") all(verdict[analysis == "primary", setNames(ibs0_best, tumour)][names(truth_map)] == truth_map) else NA,
    ibs0_null_near_zero      = if (SCENARIO == "null") pairwise[primary == TRUE & method == "custom", max(ibs0_rate)] < 0.02 else NA,
    maxdp_removes_artefacts  = nrow(art) == 0 || mean(art$too_deep) >= 0.8,
    masked_never_hom_ref     = all(body[, unlist(lapply(.SD, function(x) !startsWith(x, "0/0:."))), .SDcols = SAMPLES]),
    ambiguous_tumours_masked = geno[sample %chin% TUMOUR_IDS & status == "ambiguous", all(masked)],
    six_somalier_files       = if (HAVE_SOMALIER) all(lengths(som_files) == 6) else NA,
    somalier_masked_unknown  = if (HAVE_SOMALIER) nrow(som_n_check) > 0 && som_n_check[, all(n_custom == n_somalier)] else NA,
    outputs_written          = all(file.exists(file.path(OUT, c("filter_counts.tsv", "panel_sizes.tsv", "pairwise_metrics.tsv",
                                                                "assignment_scores.tsv", "verdict.tsv", "robustness.tsv", "RUN_LOG.md",
                                                                "paired_tests.tsv", "jackknife.tsv", "bootstrap_D.tsv", "bootstrap_assignment_support.tsv",
                                                                "segments.tsv", "ki_locus.tsv", "mindp_sensitivity.tsv", "pass_sensitivity.tsv",
                                                                "heatmap_custom_disc_rate.png", "heatmap_custom_ibs0_rate.png", "disc_per_tumour.png",
                                                                "ibs0_per_tumour.png", "robustness.png", "bootstrap_D.png", "chromosomal_distribution.png")))),
    run_log_sections         = all(vapply(c("### Changelog", "### Filter counts", "### Verdict", "### Paired sign tests", "### Block bootstrap",
                                            "### Jackknife", "### Discordant segments", "### MIN_DP sensitivity", "### FreeBayes PASS exclusion sensitivity", "### Interpretation and limitations",
                                            "### Deviations", "### Denominators"), function(s) any(startsWith(readLines(log_file), s)), TRUE)))
  message("TEST scenario ", SCENARIO, ": IBS0 max rate over all pairs = ", signif(pairwise[primary == TRUE & method == "custom", max(ibs0_rate)], 3),
          "; artefact sites at the depth step: ", nrow(art), ", removed by max-depth: ", sum(art$too_deep))
  print(checks)
  if (!all(checks, na.rm = TRUE)) stop("TEST_MODE checks failed (", SCENARIO, "): ", paste(names(checks)[!is.na(checks) & !checks], collapse = ", "))
  message("TEST_MODE ", SCENARIO, ": all checks passed", if (anyNA(checks)) paste0(" (not run: ", paste(names(checks)[is.na(checks)], collapse = ", "), ")") else "")
}
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Verdict: ", file.path(OUT, "verdict.tsv"),
        "; log: ", log_file)
