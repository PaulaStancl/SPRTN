#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24 | Rewritten: 2026-11-24
#
# mouse_wgs_somalier.R - tumour / KI-normal sample identity from the setting A FreeBayes VCFs, with
# three methods on ONE masked genotype file: somalier (Pedersen et al., Genome Medicine 2020),
# bcftools gtcheck, and a custom IBS0 count in R. No BAMs.
#
# Principle: samples of one animal share their inherited genotypes. A confident hom-ref vs hom-alt
# site (IBS0, "opposite homozygotes") is not expected between them. The primary metric is therefore
# IBS0 / jointly callable sites; IBS2 and homozygous concordance are secondary, somalier relatedness is
# exploratory. [CHANGE 21] (replaces the old claim that relatedness "assumes population allele frequencies":
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
#   tumours: same bins, ambiguous -> MASKED (./.), never hom-ref; tumours are never used to select sites
# Panels on the final site set: A all; B KI1-3 not all identical; C opposite homs (0/0 vs 1/1) in >= 2 KI
# normals. Thinning: at most one site per window (greedy minimum distance, best minimum depth kept),
# primary 20 kb, sensitivity 100 kb, 500 kb, 1 Mb. Every method runs on every panel x thinning.
#
# Blocks:
#   0  EDIT THIS BLOCK
#   1  setup, tool checks, versions, help texts, logging helpers
#   2  germline site selection, masked 6-sample genotype file, panels, thinning
#   3  somalier extract (per panel x thinning, on the masked file)
#   4  somalier relate; bcftools gtcheck; custom IBS0
#   5  matching (IBS0 primary), assignments, hypotheses, plots
#   6  outputs, RUN_LOG.md (changelog, assumptions, limitations, deviations), TEST_MODE checks
#
# Install somalier once (static binary, no env needed):
#   mkdir -p /common/WORK/pstancl/tools/somalier && curl -sSL -o /common/WORK/pstancl/tools/somalier/somalier https://github.com/brentp/somalier/releases/download/v0.3.5/somalier && chmod +x /common/WORK/pstancl/tools/somalier/somalier
#
# Run on ssi (r-mutation env: data.table, ggplot2, bcftools):
#   Rscript scripts/mouse_wgs/mouse_wgs_somalier.R           # real data
#   Rscript scripts/mouse_wgs/mouse_wgs_somalier.R --test    # synthetic data with a known HCC_1/HCC_2 swap
# Output: results/mouse_wgs/mouse_wgs_somalier/ (TEST_MODE: results/mouse_wgs/mouse_wgs_somalier_TEST/)
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

# [CHANGE 1] synthetic test: TEST_MODE <- TRUE here, or `Rscript ... --test`
TEST_MODE <- FALSE

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# Console as wide as the terminal (falls back to 80 when there is no terminal, e.g. batch Rscript)
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
started <- Sys.time()

# [CHANGE 1] test mode: own output folder; synthetic inputs are generated there in block 2
if ("--test" %in% commandArgs(trailingOnly = TRUE)) TEST_MODE <- TRUE
script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
SCRIPT_RUN <- if (length(script_arg) && file.exists(script_arg[1])) normalizePath(script_arg[1]) else SCRIPT
if (TEST_MODE) {
  base <- if (dir.exists(REPO)) REPO else getwd()
  OUT  <- Sys.getenv("SOMALIER_TEST_OUT", file.path(base, "results/mouse_wgs/mouse_wgs_somalier_TEST"))
  DIR_A <- file.path(OUT, "test_input")
  FASTA <- file.path(DIR_A, "test_genome.fa")
  KNOCKIN_LOCUS <- list(chrom = "7", pos = 500000)
  KNOCKIN_FLANK <- 2e5
  message("TEST_MODE: synthetic data, output ", OUT)
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
  gtcheck_usage <- capture(BCFTOOLS, "gtcheck")          # no arguments: prints its usage, exit status 1
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
# [CHANGE 1] synthetic test data: 3 animals (A, B, C); HCC_1 is from B but filed with KI1 (animal A),
# HCC_2 from A filed with KI2 (B), HCC_3 and KI3 both C. Ambiguous VAFs, low depth, tumour LOH, sites
# missing from one file, indels/multiallelics/low QUAL that the filters must drop.
make_test_data <- function(dir) {
  set.seed(20261124)
  dir.create(file.path(dir, "calls"), recursive = TRUE, showWarnings = FALSE)
  chroms <- as.character(1:19); L <- 1000000L       # integer: a double is written as 1e+06 in the .fai
  # FASTA + .fai written by hand (60 bases per line)
  fa <- file(FASTA, "w"); fai <- data.table(); off <- 0
  seqs <- list()
  for (ch in chroms) {
    s <- sample(c("A", "C", "G", "T"), L, replace = TRUE); seqs[[ch]] <- s
    hdr <- paste0(">", ch); writeLines(hdr, fa); off <- off + nchar(hdr) + 1
    lines <- vapply(split(s, ceiling(seq_along(s) / 60)), paste, "", collapse = "")
    writeLines(lines, fa)
    fai <- rbind(fai, data.table(ch, L, as.integer(off), 60L, 61L)); off <- off + sum(nchar(lines) + 1)
  }
  close(fa); fwrite(fai, paste0(FASTA, ".fai"), sep = "\t", col.names = FALSE)
  n_per <- 300
  st <- rbindlist(lapply(chroms, function(ch) data.table(CHROM = ch, POS = sort(sample(2:(L - 1), n_per)))))
  st[, REF := mapply(function(ch, p) seqs[[ch]][p], CHROM, POS)]
  st[, ALT := vapply(REF, function(r) sample(setdiff(c("A", "C", "G", "T"), r), 1), "")]
  g0 <- sample(0:2, nrow(st), replace = TRUE, prob = c(0.5, 0.25, 0.25))
  animal <- function() { g <- g0; k <- runif(length(g)) < 0.3; g[k] <- sample(0:2, sum(k), replace = TRUE); g }
  truth <- list(A = animal(), B = animal(), C = animal())
  who <- c(HCC_1 = "B", KI1 = "A", HCC_2 = "A", KI2 = "B", HCC_3 = "C", KI3 = "C")
  long <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1", HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
            HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3", KI1 = "sfb_dikic_2025_05_05_KI1",
            KI2 = "sfb_dikic_2025_05_06_KI2", KI3 = "sfb_dikic_2025_05_07_KI3")
  sim_sample <- function(s) {
    g <- truth[[who[[s]]]]
    if (grepl("^HCC", s)) { k <- g == 1 & runif(length(g)) < 0.05; g[k] <- sample(c(0L, 2L), sum(k), replace = TRUE) }  # LOH
    dp <- rpois(length(g), 45); low <- runif(length(g)) < 0.05; dp[low] <- rpois(sum(low), 15)
    p <- c(0.002, 0.5, 0.998)[g + 1]
    amb <- runif(length(g)) < 0.03; p[amb] <- sample(c(0.2, 0.8), sum(amb), replace = TRUE)
    ao <- rbinom(length(g), dp, p)
    data.table(g = g, DP = dp, RO = dp - ao, AO = ao, ambiguous_truth = amb)
  }
  pairs <- list(c("HCC_1", "KI1"), c("HCC_2", "KI2"), c("HCC_3", "KI3"))
  for (pr in pairs) {
    a <- sim_sample(pr[1]); b <- sim_sample(pr[2])
    keep <- (a$g > 0 | b$g > 0) & runif(nrow(st)) > 0.10          # variant in the pair, 10% missing from this file
    gt_txt <- function(x) fifelse(x$AO / pmax(x$DP, 1) < 0.2, "0/0", fifelse(x$AO / pmax(x$DP, 1) > 0.8, "1/1", "0/1"))
    fmt <- function(x) paste(gt_txt(x), x$DP, paste(x$RO, x$AO, sep = ","), x$RO, x$RO * 30, x$AO, x$AO * 30, "-10,-1,0", sep = ":")
    rec <- data.table(st, ID = ".", QUAL = fifelse(runif(nrow(st)) < 0.05, 15, 100), FILTER = ".", INFO = "TYPE=snp",
                      FORMAT = "GT:DP:AD:RO:QR:AO:QA:GL", S1 = fmt(a), S2 = fmt(b))[keep]
    junk <- rec[sample(.N, 30)][, `:=`(ALT = fifelse(seq_len(.N) %% 2 == 0, paste0(ALT, "T"), paste0(ALT, ",G")), INFO = "TYPE=other")]
    rec <- rbind(rec, junk)[order(as.integer(CHROM), POS)]
    setcolorder(rec, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2"))
    f <- file.path(dir, "calls", paste0(long[[pr[1]]], "_vs_", long[[pr[2]]], ".freebayes.ann.norm.vcf"))
    writeLines(c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s,length=%d>", chroms, as.integer(L)),
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
  invisible(TRUE)
}
if (TEST_MODE) { unlink(DIR_A, recursive = TRUE); dir.create(DIR_A, recursive = TRUE); make_test_data(DIR_A) }

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
  cn <- c("CHROM", "POS", "REF", "ALT", "QUAL", as.vector(t(outer(smp, fields, paste, sep = "."))))
  if (HAVE_BCFTOOLS) {
    q <- paste0("%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL[", paste0("\\t%", fields, collapse = ""), "]\\n")
    cmd <- paste(B, "view -v snps -m2 -M2 -i", shQuote(paste0("QUAL>=", MIN_QUAL)), "-t", paste(autos, collapse = ","),
                 shQuote(path), "|", B, "query -f", shQuote(q))
    x <- fread(cmd = cmd, sep = "\t", header = FALSE, col.names = cn, colClasses = "character")
  } else {
    # [CHANGE 1] TEST_MODE without bcftools: same filters in R on the VCF text
    v <- fread(path, skip = "#CHROM", sep = "\t", header = TRUE, colClasses = "character")
    setnames(v, 1, "CHROM")
    v <- v[CHROM %chin% autos & nchar(REF) == 1 & grepl("^[ACGT]$", ALT) & as.numeric(QUAL) >= MIN_QUAL]
    fk <- strsplit(v$FORMAT[1], ":", fixed = TRUE)[[1]]
    x <- v[, .(CHROM, POS, REF, ALT, QUAL)]
    for (j in seq_along(smp)) {
      parts <- tstrsplit(v[[9 + j]], ":", fixed = TRUE)
      for (fld in fields) set(x, j = paste(smp[j], fld, sep = "."), value = parts[[match(fld, fk)]])
    }
  }
  count_step(paste0("biallelic SNV, autosome, QUAL>=", MIN_QUAL, " (records read)"), nrow(x), scope = paste(smp, collapse = "+"))
  rbindlist(lapply(smp, function(s) {
    ad <- tstrsplit(x[[paste0(s, ".AD")]], ",", fixed = TRUE)
    d <- data.table(CHROM = x$CHROM, POS = as.integer(x$POS), REF = x$REF, ALT = x$ALT, QUAL = as.numeric(x$QUAL),
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

# [CHANGE 8] genotype bins: 0 hom-ref, 1 het, 2 hom-alt, NA = low depth / ambiguous / absent
ad[, VAF := fifelse(DPA > 0, ALT_N / DPA, NA_real_)]
ad[, status := fcase(is.na(DPA) | DPA < MIN_DP, "low_depth",
                     VAF <= HOM_REF_MAX, "hom_ref", VAF >= HET_MIN & VAF <= HET_MAX, "het", VAF >= HOM_ALT_MIN, "hom_alt",
                     default = "ambiguous")]
ad[, g := fcase(status == "hom_ref", 0L, status == "het", 1L, status == "hom_alt", 2L, default = NA_integer_)]

# Site x sample, all six samples (absent from a file = "absent")
ad[, key := paste(CHROM, POS, REF, ALT, sep = ":")]
sites_all <- unique(ad[, .(key, CHROM, POS, REF, ALT)])
full <- CJ(key = sites_all$key, sample = SAMPLES)[ad[, .(key, sample, REF_N, ALT_N, DPA, VAF, status, g)], on = .(key, sample),
                                                  `:=`(REF_N = i.REF_N, ALT_N = i.ALT_N, DPA = i.DPA, VAF = i.VAF, status = i.status, g = i.g)]
full[is.na(status), status := "absent"]
count_step("union of sites over the 3 files", nrow(sites_all))

st <- full[, .(n_present = sum(status != "absent"), n_dp = sum(!is.na(DPA) & DPA >= MIN_DP),
               min_dp = suppressWarnings(min(DPA, na.rm = TRUE)),
               normal_ambiguous = any(sample %chin% NORMAL_IDS & status %chin% c("ambiguous", "low_depth", "absent")),
               normal_nonref = any(sample %chin% NORMAL_IDS & g %in% 1:2)), by = key]
st <- sites_all[st, on = "key"]
count_step("present in all 3 files (all six samples have a record)", st[n_present == 6, .N])
# [CHANGE 7] ascertainment: depth in ALL six samples (old rule: >= 4 of 6)
n_old_rule <- st[n_dp >= 4, .N]
count_step(paste0("REF+ALT >= ", MIN_DP, " in all six samples"), st[n_dp == 6, .N],
           note = sprintf("old rule (>= 4 of 6 samples) would keep %s; the new rule removes %s of those", n_old_rule, n_old_rule - st[n_dp == 6, .N]))
count_step("+ no KI normal ambiguous (VAF in 0.05-0.35 or 0.65-0.95)", st[n_dp == 6 & !normal_ambiguous, .N])
count_step("+ >= 1 KI normal confidently non-ref (= panel A before thinning)", st[n_dp == 6 & !normal_ambiguous & normal_nonref, .N])
final <- st[n_dp == 6 & !normal_ambiguous & normal_nonref]
if (!nrow(final)) stop("no site passes the germline selection - see filter_counts in the log")

# [CHANGE 8] masked genotypes on the final sites: tumours ambiguous -> NA (masked), never hom-ref
geno <- full[key %chin% final$key]
geno[, masked := is.na(g)]
masked_summary <- geno[, .(sites = .N, masked = sum(masked), masked_ambiguous = sum(status == "ambiguous"),
                           hom_ref = sum(g == 0, na.rm = TRUE), het = sum(g == 1, na.rm = TRUE), hom_alt = sum(g == 2, na.rm = TRUE)), by = sample]
print(masked_summary)
stopifnot(geno[sample %chin% NORMAL_IDS, !any(masked)])          # normals are confident by construction
G <- dcast(geno, key ~ sample, value.var = "g")[final[, .(key, CHROM, POS, REF, ALT, min_dp)], on = "key"]
setcolorder(G, c("key", "CHROM", "POS", "REF", "ALT", "min_dp", SAMPLES))
G[, chr_n := as.integer(sub("^chr", "", CHROM))]
setorder(G, chr_n, POS)
G[, AF := round(rowSums(.SD) / (2 * length(NORMAL_IDS)), 4), .SDcols = NORMAL_IDS]   # in-sample, KI normals only

# [CHANGE 10] panels
G[, `:=`(n_min = do.call(pmin, .SD), n_max = do.call(pmax, .SD)), .SDcols = NORMAL_IDS]
G[, `:=`(panel_B = n_min != n_max, panel_C = n_min == 0L & n_max == 2L)]      # B: not all identical; C: 0/0 and 1/1 present
G[, c("n_min", "n_max") := NULL]
panel_sites <- list(A = G$key, B = G[panel_B == TRUE, key], C = G[panel_C == TRUE, key])

# [CHANGE 11] greedy minimum-distance thinning: best minimum depth first; a site is kept only if no kept
# site on the same chromosome is closer than `dist`
thin_greedy <- function(d, dist) {
  d <- d[order(-min_dp, chr_n, POS)]
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

# ---- 4. somalier relate, bcftools gtcheck, custom IBS0 ------------------------------
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
  p[, .(combo = cb, s1 = sample_a, s2 = sample_b, method = "somalier", n_joint = n, ibs0, ibs2,
        relatedness, somalier_concordance = concordance, hets_a, hets_b, shared_hets)]
}
som <- if (HAVE_SOMALIER) rbindlist(lapply(names(som_files), somalier_relate)) else data.table()

# [CHANGE 15] gtcheck: only options listed by the installed version's usage
gtcheck_ok <- HAVE_BCFTOOLS && any(grepl("-G, --GTs-only", gtcheck_usage, fixed = TRUE)) && any(grepl("-T, --targets-file", gtcheck_usage, fixed = TRUE))
if (HAVE_BCFTOOLS && !gtcheck_ok) deviation("bcftools gtcheck usage lacks -G/--GTs-only or -T/--targets-file: gtcheck skipped (usage in tool_help/)")
gtcheck_run <- function(cb) {
  raw <- file.path(OUT, "gtcheck", paste0("gtcheck_", cb, ".txt"))
  # -G 1: GTs only (the masked file has no PL/GL anyway); masked ./. are skipped by gtcheck's GT code
  run(paste(B, "gtcheck -G 1 -T", shQuote(file.path(OUT, "sites", paste0("targets_", cb, ".tsv"))),
            shQuote(masked_file), ">", shQuote(raw), "2>", shQuote(paste0(raw, ".stderr"))))
  l <- readLines(raw)
  # 1.4.1 cross-check output: "# ERR, error rate\t[2]Pairwise error rate\t[3]Number of sites compared\t[4]Sample i\t[5]Sample j"
  if (!any(grepl("^# ERR.*Pairwise error rate.*Number of sites compared.*Sample i.*Sample j", l))) {
    deviation("gtcheck ", cb, ": output is not the 1.4.1 ERR format - raw output kept in ", raw, ", not parsed"); return(NULL)
  }
  e <- grep("^ERR\t", l, value = TRUE)
  if (!length(e)) { deviation("gtcheck ", cb, ": no ERR lines"); return(NULL) }
  x <- as.data.table(tstrsplit(e, "\t", fixed = TRUE))[, .(err = as.numeric(V2), n = as.integer(V3), s1 = V4, s2 = V5)]
  x[, .(combo = cb, s1, s2, method = "gtcheck", n_joint = n, gtcheck_discordance = err,
        gtcheck_n_discordant = round(err * n))]
}
gtc <- if (gtcheck_ok) rbindlist(lapply(combos[n_sites > 0, combo], gtcheck_run)) else data.table()

# [CHANGE 16] custom IBS0 on the masked genotypes
custom_metrics <- function(cb) {
  d <- G[key %chin% combo_sites[[cb]]]
  rbindlist(lapply(seq_len(nrow(pair_dt)), function(i) {
    a <- d[[pair_dt$s1[i]]]; b <- d[[pair_dt$s2[i]]]
    j <- !is.na(a) & !is.na(b); hb <- j & a != 1 & b != 1
    data.table(combo = cb, s1 = pair_dt$s1[i], s2 = pair_dt$s2[i], method = "custom", n_sites = nrow(d),
               missing_s1 = mean(is.na(a)), missing_s2 = mean(is.na(b)), n_joint = sum(j),
               ibs0 = sum(j & abs(a - b) == 2), ibs2 = sum(j & a == b),
               n_hom_both = sum(hb), hom_identical = sum(hb & a == b))
  }))
}
cust <- rbindlist(lapply(combos[n_sites > 0, combo], custom_metrics))

# One long table for every method; rates and exact 95% CIs on counts
canon <- function(x) { if (!nrow(x)) return(x); sw <- match(x$s1, SAMPLES) > match(x$s2, SAMPLES)
  x[sw, c("s1", "s2") := .(s2, s1)]; x }
cust <- canon(cust); som <- canon(som); gtc <- canon(gtc)
add_rates <- function(x) { if (!nrow(x)) return(x)
  x[, `:=`(ibs0_rate = ibs0 / n_joint, ibs2_rate = ibs2 / n_joint, ibs0_lo = ci_lo(ibs0, n_joint), ibs0_hi = ci_hi(ibs0, n_joint))] }
cust <- add_rates(cust); som <- add_rates(som)
cust[, `:=`(hom_concordance = fifelse(n_hom_both > 0, hom_identical / n_hom_both, NA_real_))]
if (nrow(gtc)) gtc[, `:=`(gt_lo = ci_lo(gtcheck_n_discordant, n_joint), gt_hi = ci_hi(gtcheck_n_discordant, n_joint))]
pairwise <- rbindlist(list(cust, som, gtc), fill = TRUE)
for (cl in c("ibs0_rate", "ibs2_rate", "ibs0_lo", "ibs0_hi", "hom_concordance", "relatedness", "somalier_concordance",
             "gtcheck_discordance", "gt_lo", "gt_hi")) if (!cl %in% names(pairwise)) pairwise[, (cl) := NA_real_]   # methods not run
pairwise[, `:=`(type = pair_type(s1, s2), labelled_pair = is_pair(s1, s2, KI_NORMALS), hypothesis_pair = is_pair(s1, s2, HYPOTHESIS))]
pairwise <- combos[, .(combo, panel, thin, underpowered, primary)][pairwise, on = "combo"]
pairwise[primary == TRUE & method == "custom", .(s1, s2, type, n_joint, ibs0, ibs0_rate = round(ibs0_rate, 4),
                                                 ibs2_rate = round(ibs2_rate, 3), hom_concordance = round(hom_concordance, 3), labelled_pair)][order(type, ibs0_rate)]

# ---- 5. Matching, assignments, plots ----------------------------------------------
# [CHANGE 17] per tumour: lowest IBS0 rate first, highest IBS2 rate second (gtcheck: lowest discordance);
# ambiguous when best and second-best CIs overlap, or too few joint sites, or the panel is underpowered
tn <- pairwise[type == "tumour-normal"]
tn[, `:=`(tumour = fifelse(s1 %chin% TUMOUR_IDS, s1, s2), normal = fifelse(s1 %chin% TUMOUR_IDS, s2, s1))]
tn[, score := fifelse(method == "gtcheck", gtcheck_discordance, ibs0_rate)]
tn[, score_lo := fifelse(method == "gtcheck", gt_lo, ibs0_lo)]
tn[, score_hi := fifelse(method == "gtcheck", gt_hi, ibs0_hi)]
tn[, tiebreak := fifelse(method == "gtcheck", 0, -ibs2_rate)]
setorder(tn, combo, method, tumour, score, tiebreak)
robust <- tn[, .(best_normal = normal[1], second_normal = normal[2], best_score = score[1], second_score = score[2],
                 best_hi = score_hi[1], second_lo = score_lo[2], n_joint_best = n_joint[1], underpowered = underpowered[1]),
             by = .(combo, panel, thin, method, tumour)]
robust[, ci_overlap := best_hi >= second_lo]
robust[, ambiguous := ci_overlap | n_joint_best < MIN_JOINT_SITES | underpowered]
robust[, `:=`(labelled_normal = KI_NORMALS[tumour], hypothesis_normal = HYPOTHESIS[tumour])]
robust[, `:=`(matches_labelled = best_normal == labelled_normal, matches_hypothesis = best_normal == hypothesis_normal)]
robust[, call := fifelse(ambiguous, paste0(best_normal, "?"), best_normal)]
dcast(robust, method + tumour ~ combo, value.var = "call")

# [CHANGE 17] all 3! one-to-one tumour -> normal assignments, scored by the summed score of their 3 pairs
perm3 <- function(v) { r <- list(); for (i in v) for (j in setdiff(v, i)) r[[length(r) + 1]] <- c(i, j, setdiff(v, c(i, j))); r }
assign_scores <- rbindlist(lapply(split(tn, by = c("combo", "method")), function(d) {
  rbindlist(lapply(perm3(NORMAL_IDS), function(pm) {
    m <- setNames(pm, TUMOUR_IDS)
    sel <- d[normal == m[tumour]]
    var_sum <- sum(sel[, score * (1 - score) / pmax(n_joint, 1)])
    data.table(combo = d$combo[1], method = d$method[1], assignment = paste(paste0(TUMOUR_IDS, "-", m), collapse = " "),
               total_score = sum(sel$score), var_sum = var_sum,
               is_labelled = all(m == KI_NORMALS[TUMOUR_IDS]), is_hypothesis = all(m == HYPOTHESIS[TUMOUR_IDS]))
  }))
}))
setorder(assign_scores, combo, method, total_score)
assign_scores[, rank := seq_len(.N), by = .(combo, method)]
# separated: gap to the next assignment > 2 x sqrt(summed binomial variances of both) (approximate z = 2)
assign_scores[, gap_to_next := shift(total_score, -1) - total_score, by = .(combo, method)]
assign_scores[, clearly_separated := gap_to_next > 2 * sqrt(var_sum + shift(var_sum, -1)), by = .(combo, method)]
assign_scores[rank <= 2 & combo == combos[primary == TRUE, combo]]

# Verdict on the primary combo (panel A, THIN_PRIMARY), custom IBS0, ranked by IBS0
PRIMARY <- combos[primary == TRUE, combo]
verdict <- robust[combo == PRIMARY & method == "custom",
                  .(tumour, labelled_normal, best_normal, second_normal, best_ibs0_rate = best_score, best_ibs0_ci_hi = best_hi,
                    second_ibs0_rate = second_score, second_ibs0_ci_lo = second_lo, n_joint_best, ambiguous,
                    consistent_with_label = matches_labelled, matches_hypothesis)]
lab_rate <- tn[combo == PRIMARY & method == "custom" & normal == KI_NORMALS[tumour], .(tumour, labelled_ibs0_rate = ibs0_rate, labelled_n_joint = n_joint)]
verdict <- lab_rate[verdict, on = "tumour"]
for (m in c("somalier", "gtcheck")) {
  r <- robust[combo == PRIMARY & method == m, .(tumour, b = best_normal, a = ambiguous)]
  if (nrow(r)) verdict[r, on = "tumour", c(paste0(m, "_best"), paste0(m, "_ambiguous")) := .(i.b, i.a)]
}
setorder(verdict, best_ibs0_rate)
verdict

# Plots ------------------------------------------------------------------------------
supported <- robust[combo == PRIMARY & method == "custom" & !ambiguous, .(tumour, normal = best_normal)]
heat_metric <- function(meth, col, label, digits = 3, lower_better = TRUE) {
  d <- pairwise[combo == PRIMARY & method == meth & !is.na(get(col)), .(s1, s2, v = get(col), n_joint)]
  if (!nrow(d)) return(invisible(NULL))
  d <- rbind(d, d[, .(s1 = s2, s2 = s1, v, n_joint)])
  d[, `:=`(s1 = factor(s1, SAMPLES), s2 = factor(s2, rev(SAMPLES)))]
  ol <- rbind(data.table(tumour = names(KI_NORMALS), normal = unname(KI_NORMALS), outline = "labelled pair"),
              supported[, .(tumour, normal, outline = "genetically supported")])
  ol <- rbind(ol[, .(s1 = tumour, s2 = normal, outline)], ol[, .(s1 = normal, s2 = tumour, outline)])
  ol[, `:=`(s1 = factor(s1, SAMPLES), s2 = factor(s2, rev(SAMPLES)))]
  p <- ggplot(d, aes(s1, s2)) + geom_tile(aes(fill = v), colour = "white") +
    geom_text(aes(label = sprintf(paste0("%.", digits, "f\nn=%d"), v, n_joint)), size = 3) +
    geom_tile(data = ol, aes(colour = outline, linetype = outline), fill = NA, linewidth = 1.1) +
    scale_colour_manual(values = c(`labelled pair` = "black", `genetically supported` = "#C62828"), name = NULL) +
    scale_linetype_manual(values = c(`labelled pair` = "dashed", `genetically supported` = "solid"), name = NULL) +
    scale_fill_gradient(low = if (lower_better) "#08519C" else "white", high = if (lower_better) "white" else "#08519C", name = label) +
    labs(x = NULL, y = NULL, title = paste0(label, " - ", meth, ", ", PRIMARY),
         subtitle = "n = jointly callable sites for that method. Outlines: labelled pairing (dashed) vs best IBS0 match (red)") +
    coord_equal() + theme_bw()
  save_plot(p, paste0("heatmap_", meth, "_", col), 8, 7)
}
heat_metric("custom", "ibs0_rate", "IBS0 / joint", 4)
heat_metric("custom", "ibs2_rate", "IBS2 / joint", 3, lower_better = FALSE)
heat_metric("custom", "hom_concordance", "homozygous concordance", 3, lower_better = FALSE)
heat_metric("somalier", "relatedness", "somalier relatedness (exploratory)", 2, lower_better = FALSE)
heat_metric("somalier", "somalier_concordance", "somalier concordance", 2, lower_better = FALSE)
heat_metric("gtcheck", "gtcheck_discordance", "gtcheck discordance", 3)

# IBS0 rate per tumour vs each KI normal with 95% CI (replaces the relatedness bars) [CHANGE 18]
pb <- tn[combo == PRIMARY & method == "custom"]
pb[, label_type := fcase(normal == KI_NORMALS[tumour], "labelled normal", normal == HYPOTHESIS[tumour], "hypothesis normal", default = "other normal")]
p_ibs0 <- ggplot(pb, aes(normal, ibs0_rate, fill = label_type)) +
  geom_col(width = 0.7) + geom_errorbar(aes(ymin = ibs0_lo, ymax = ibs0_hi), width = 0.2) +
  geom_text(aes(y = ibs0_hi, label = sprintf("%d / %d", ibs0, n_joint)), vjust = -0.4, size = 3.2) +
  facet_wrap(~ tumour, labeller = labeller(tumour = function(x) paste("Tumour", x))) +
  scale_fill_manual(values = c(`labelled normal` = "#C62828", `hypothesis normal` = "#1565C0", `other normal` = "grey60"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "KI normal liver", y = "IBS0 / jointly callable sites (95% CI)",
       title = "Which KI normal matches each tumour? Lowest IBS0 = same animal",
       subtitle = sprintf("Panel A, thinned %s: %s sites. Masked tumour genotypes excluded.", thin_label(THIN_PRIMARY), format(combos[primary == TRUE, n_sites], big.mark = ","))) +
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
       title = "Best-matching KI normal per tumour, every panel x thinning x method",
       subtitle = "? = ambiguous (CIs of best and second overlap, too few joint sites, or underpowered panel). Panels B/C are defined from normal differences.") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_rob, "robustness", 12, 2 + 1.8 * uniqueN(rb$method))

# Chromosomal distribution: panel B/C sites, and IBS0 sites per tumour-normal comparison [CHANGE 18]
bc <- rbind(G[panel_B == TRUE, .(chr_n, POS, set = "panel B")], G[panel_C == TRUE, .(chr_n, POS, set = "panel C")])
ibs0_sites <- rbindlist(lapply(TUMOUR_IDS, function(t) rbindlist(lapply(NORMAL_IDS, function(n) {
  d <- G[key %chin% combo_sites[[PRIMARY]]]
  d[!is.na(get(t)) & abs(get(t) - get(n)) == 2, .(chr_n, POS, comparison = paste(t, "-", n),
                                                  labelled = KI_NORMALS[t] == n)]
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
knockin_tab <- if (!is.na(KNOCKIN_LOCUS$chrom)) ibs0_sites[, .(ibs0_sites = .N,
  near_knockin = sum(chr_n == as.integer(sub("^chr", "", KNOCKIN_LOCUS$chrom)) & abs(POS - KNOCKIN_LOCUS$pos) <= KNOCKIN_FLANK)),
  by = .(comparison, labelled)][, frac_near_knockin := round(near_knockin / ibs0_sites, 3)] else data.table()
ibs0_by_chr <- ibs0_sites[, .N, by = .(comparison, chr_n)][order(comparison, chr_n)]
knockin_tab

# ---- 6. Outputs, RUN_LOG.md ------------------------------------------------------
fwrite(FILTER_COUNTS, file.path(OUT, "filter_counts.tsv"), sep = "\t")
fwrite(ad_check, file.path(OUT, "ad_check.tsv"), sep = "\t")
fwrite(masked_summary, file.path(OUT, "masked_genotypes_per_sample.tsv"), sep = "\t")
fwrite(combos[, .(combo, panel, thin, panel_unthinned, n_sites, underpowered, primary)], file.path(OUT, "panel_sizes.tsv"), sep = "\t")
fwrite(pairwise, file.path(OUT, "pairwise_metrics.tsv"), sep = "\t")
fwrite(assign_scores, file.path(OUT, "assignment_scores.tsv"), sep = "\t")
fwrite(verdict, file.path(OUT, "verdict.tsv"), sep = "\t")
fwrite(robust, file.path(OUT, "robustness.tsv"), sep = "\t")
fwrite(ibs0_by_chr, file.path(OUT, "ibs0_sites_per_chromosome.tsv"), sep = "\t")
if (nrow(knockin_tab)) fwrite(knockin_tab, file.path(OUT, "ibs0_near_knockin.tsv"), sep = "\t")
fwrite(G[key %chin% combo_sites[[PRIMARY]], c("CHROM", "POS", "REF", "ALT", "min_dp", "AF", SAMPLES), with = FALSE],
       file.path(OUT, "sites_used_primary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
git_hash <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(dirname(SCRIPT_RUN)), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)),
                     error = function(e) "unknown")
if (!length(git_hash)) git_hash <- "unknown"

# [CHANGE 20] RUN_LOG.md: changelog (copied from the header), assumptions, limitations, deviations
hdr_lines <- readLines(SCRIPT_RUN)
changelog <- hdr_lines[seq(grep("^# CHANGELOG", hdr_lines)[1], grep("^# END CHANGELOG", hdr_lines)[1])]
md_table <- function(x) { x <- as.data.table(x); if (!nrow(x)) return("(none)")
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))) }
fmt_num <- function(d) { d <- copy(d); for (cl in names(d)) if (is.double(d[[cl]])) set(d, j = cl, value = signif(d[[cl]], 4)); d }
log_entry <- c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_somalier", if (TEST_MODE) " (TEST_MODE, synthetic data)"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Tools:      ", somalier_version, "; ", paste(bcftools_version, collapse = " / ")),
  paste0("- Tool help:  ", if (length(help_files)) paste(help_files, collapse = ", ") else "(not saved: tools missing)"),
  paste0("- Script:     ", SCRIPT_RUN, "  git: ", git_hash),
  sprintf("- Parameters: MIN_QUAL=%s MIN_DP=%s bins hom_ref<=%s het %s-%s hom_alt>=%s; thinning %s (primary) + %s; MIN_PANEL_SITES=%s MIN_JOINT_SITES=%s; KNOCKIN=%s",
          MIN_QUAL, MIN_DP, HOM_REF_MAX, HET_MIN, HET_MAX, HOM_ALT_MIN, thin_label(THIN_PRIMARY), paste(thin_label(THIN_SENS), collapse = ", "),
          MIN_PANEL_SITES, MIN_JOINT_SITES, if (is.na(KNOCKIN_LOCUS$chrom)) "NA" else paste0(KNOCKIN_LOCUS$chrom, ":", format(KNOCKIN_LOCUS$pos, scientific = FALSE), " +-", format(KNOCKIN_FLANK, scientific = FALSE))),
  paste0("- Labelled:   ", paste(names(KI_NORMALS), KI_NORMALS, sep = "-", collapse = ", "), "; hypothesis tested: ", paste(names(HYPOTHESIS), HYPOTHESIS, sep = "-", collapse = ", ")),
  paste0("- Inputs:     ", paste(basename(fb_files), collapse = ", "), "; FASTA ", FASTA),
  paste0("- Outputs:    ", OUT),
  "", "### Changelog of the script", "", sub("^# ?", "", changelog),
  "", "### Filter counts", "", md_table(FILTER_COUNTS),
  "", "### AD check (AD vs RO,AO)", "", md_table(ad_check),
  "", "### Masked genotypes per sample (final sites)", "", md_table(masked_summary),
  "", "### Panel sizes", "", md_table(combos[, .(combo, panel_unthinned, n_sites, underpowered, primary)]),
  "", paste0("### Verdict (primary: panel A, thinned ", thin_label(THIN_PRIMARY), ", custom IBS0)"), "", md_table(fmt_num(verdict)),
  "", "### One-to-one assignments, primary, custom IBS0", "",
  md_table(fmt_num(assign_scores[combo == PRIMARY & method == "custom", .(rank, assignment, total_score, gap_to_next, clearly_separated, is_labelled, is_hypothesis)])),
  "", "### Robustness (best normal; ? = ambiguous)", "", md_table(dcast(robust, method + tumour ~ combo, value.var = "call")),
  "", "### Denominators and definitions",
  "- custom: n_joint = sites where both samples have an unmasked genotype; IBS0 = 0/0 vs 1/1; IBS2 = identical genotype; homozygous concordance = identical among sites where both are homozygous. 95% CIs: exact binomial (binom.test).",
  "- somalier: ibs0/n and ibs2/n with somalier's n = \"the number of sites where both samples have a known genotype\" (its genotypes, called from AD; v0.3.5 source: N counts sites known in both). relatedness = 2*(shared_hets - 2*ibs0)/het_ab (v0.3.5 source); `concordance` is somalier's adjusted homozygous concordance (stretched, with a low hom-alt penalty), not our homozygous concordance.",
  paste0("- somalier genotype calls (v0.3.5 source, relate.nim `alts`): unknown if REF+ALT < -d; hom-ref if ALT = 0 or AB < cutoff; hom-alt if REF = 0 or AB > 1 - cutoff; het if min-ab <= AB <= 1 - min-ab; else unknown. Run with -d ", MIN_DP, ", --min-ab ", HET_MIN, ", SOMALIER_AB_HOM_CUTOFF=", AB_HOM_CUTOFF, " so they equal our bins; masked genotypes have AD '.', read as depth 0 -> unknown. -u (unknown -> hom-ref) is never used."),
  "- gtcheck (1.4.1, no -g: cross-check of all samples): -G 1 uses GTs only; a pair is compared at sites where both GTs are non-missing; ERR = fraction of those sites with different genotypes.",
  "", "### Assumptions",
  "- Each FreeBayes file was called on its own labelled pair; a site missing from a file is unknown, never hom-ref.",
  "- VAF bins give an inferred genotype under a diploid model; FreeBayes' own GT is not used.",
  "- INFO/AF in the sites and masked files is computed in-sample from the 3 KI normals; somalier uses it only for its per-sample hom-ref/hom-alt counts and contamination estimate (per its help text), not for the pair statistics used here.",
  "- KI normals define the sites; tumours are only required to be callable (depth) and are never used to choose sites.",
  "", "### Interpretation and limitations",
  "- There are no human relatedness thresholds for inbred, genetically modified mice; genetic concordance shows which samples share inherited genotypes, it is not definitive proof of individual identity.",
  "- somalier, gtcheck and the custom IBS0 all use the same FreeBayes calls and the same masked file: they are complementary, not independent evidence.",
  "- Ascertainment: the per-pair BCFs contain only sites variant in that pair, so requiring all six samples keeps sites variant in all three pairs; this removes sites private to one animal and favours shared variation.",
  "- Tumours: purity, LOH and copy number make tumour hets unreliable. IBS0 is robust to LOH (het -> hom gives IBS1, not IBS0), IBS2 is not. Ambiguous tumour genotypes are masked.",
  "- Panels B and C are defined from differences between the KI normals, so normal-normal metrics on them are inflated by construction; only tumour-normal comparisons are interpretable there.",
  "- Residual inherited haplotype blocks (including donor-strain DNA flanking the knock-in) can concentrate IBS0 in few regions; thinning reduces but does not remove this dependence.",
  "- If results are ambiguous or panels underpowered: genotype the panel A sites from the BAMs in all six samples (e.g. bcftools mpileup -T sites) once BAMs are available.",
  "", "### Deviations", "", if (length(DEVIATIONS)) paste0("- DEVIATION: ", DEVIATIONS) else "- none", "")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(log_entry, old_log), log_file)

# [CHANGE 1] TEST_MODE checks
if (TEST_MODE) {
  truth <- c(HCC_1 = "KI2", HCC_2 = "KI1", HCC_3 = "KI3")
  got <- setNames(verdict$best_normal, verdict$tumour)[names(truth)]
  checks <- c(
    swap_detected            = all(got == truth),
    masked_never_hom_ref     = all(body[, unlist(lapply(.SD, function(x) !startsWith(x, "0/0:."))), .SDcols = SAMPLES]),
    ambiguous_tumours_masked = geno[sample %chin% TUMOUR_IDS & status == "ambiguous", all(masked)],
    six_somalier_files       = if (HAVE_SOMALIER) all(lengths(som_files) == 6) else NA,
    outputs_written          = all(file.exists(file.path(OUT, c("filter_counts.tsv", "panel_sizes.tsv", "pairwise_metrics.tsv",
                                                                "assignment_scores.tsv", "verdict.tsv", "robustness.tsv", "RUN_LOG.md",
                                                                "heatmap_custom_ibs0_rate.png", "ibs0_per_tumour.png", "robustness.png",
                                                                "chromosomal_distribution.png")))),
    run_log_sections         = all(vapply(c("### Changelog", "### Filter counts", "### Verdict", "### Interpretation and limitations",
                                            "### Deviations", "### Denominators"), function(s) any(startsWith(readLines(log_file), s)), TRUE)))
  print(checks)
  if (!all(checks, na.rm = TRUE)) stop("TEST_MODE checks failed: ", paste(names(checks)[!is.na(checks) & !checks], collapse = ", "))
  message("TEST_MODE: all checks passed", if (anyNA(checks)) paste0(" (not run: ", paste(names(checks)[is.na(checks)], collapse = ", "), ")") else "")
}
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Verdict: ", file.path(OUT, "verdict.tsv"),
        "; log: ", log_file)
