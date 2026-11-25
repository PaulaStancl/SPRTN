#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_freebayes_pass_check.R - do the FreeBayes PASS records satisfy the provider's prefilter?
#
#   INFO/MQM >= 40 && GT[normal]="0/0" && GT[tumor]="0/1" && FMT/DP[normal] >= 30 && FMT/DP[tumor] >= 30
#   && QUAL/INFO/AO >= 5 && SAF > 0 && SAR > 0
#
# The expression is evaluated here in R from the raw fields of each FreeBayes file (bcftools 1.4.1
# cannot address samples by name in -i), then compared with FILTER in both directions:
#   PASS but rule FALSE  -> PASS records that break the prefilter (which condition)
#   not PASS but rule TRUE -> records the prefilter would pass but FILTER does not
# GT is compared as a string, as in the expression ("0|1" would not count as "0/1"; counted separately).
# Multiallelic records (should be none, files split with norm -m -any): Number=A fields split with ALT.
#
# Blocks: 0 EDIT | 1 setup | 2 files and samples | 3 read and evaluate | 4 compare with FILTER | 5 outputs
# Env on ssi: r-mutation. Output: results/mouse_wgs/mouse_wgs_freebayes_pass_check/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
DIR_A    <- file.path(REPO, "data/mouse/wgs/tumor_vs_normal/snv_indel")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_freebayes_pass_check")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_freebayes_pass_check.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
TUMOURS  <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1",
              HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
              HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
bcf <- function(args) {
  out <- suppressWarnings(system2(BCFTOOLS, args, stdout = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) stop("bcftools failed: ", paste(args, collapse = " "))
  out
}

# ---- 2. FreeBayes files (single-caller, not the merged file) and their samples ----
files <- data.table(path = list.files(DIR_A, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf|vcf\\.gz)$"))
files <- files[grepl("freebayes", path, ignore.case = TRUE)]
files[, rel := substring(path, nchar(DIR_A) + 2)]
files[, has_CALLER := sapply(path, function(f) any(startsWith(bcf(c("view", "-h", shQuote(f))), "##INFO=<ID=CALLER,")))]
files <- files[has_CALLER == FALSE]
files[, tumour := NA_character_]
for (k in names(TUMOURS)) files[grepl(TUMOURS[[k]], basename(rel), fixed = TRUE), tumour := k]
files[, samples := sapply(path, function(f) paste(bcf(c("query", "-l", shQuote(f))), collapse = ","))]
print(files[, .(tumour, rel, samples)])
stopifnot(nrow(files) > 0, !anyNA(files$tumour))

# ---- 3. Read fields and evaluate the expression --------------------------------
read_fb <- function(i) {
  f <- files[i]
  smp <- strsplit(f$samples, ",", fixed = TRUE)[[1]]
  it  <- grep(TUMOURS[[f$tumour]], smp, fixed = TRUE)
  if (length(smp) != 2 || length(it) != 1) stop(f$rel, ": expected 2 samples, one of them the tumour; got ", f$samples)
  fmt <- "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%QUAL\\t%INFO/MQM\\t%INFO/AO\\t%INFO/SAF\\t%INFO/SAR[\\t%GT\\t%DP]\\n"
  cn  <- c("CHROM", "POS", "REF", "ALT", "FILTER", "QUAL", "MQM", "AO", "SAF", "SAR",
           paste0(c("GT_", "DP_"), rep(c("s1", "s2"), each = 2)))
  x <- fread(cmd = paste(shQuote(BCFTOOLS), "query -f", shQuote(fmt), shQuote(f$path)), sep = "\t", header = FALSE,
             col.names = cn, colClasses = "character", quote = "", na.strings = ".")
  tum <- paste0("s", it); nor <- paste0("s", 3 - it)
  setnames(x, paste0(c("GT_", "DP_", "GT_", "DP_"), c(tum, tum, nor, nor)), c("GT_tumour", "DP_tumour", "GT_normal", "DP_normal"))
  # Number=A fields with ALT (multiallelic records only)
  multi <- x[ALT %like% ",", .N]
  if (multi) {
    x[, rid := .I]
    x <- x[, .(ALT = strsplit(ALT, ",")[[1]], AO = strsplit(AO, ",")[[1]], MQM = strsplit(MQM, ",")[[1]],
               SAF = strsplit(SAF, ",")[[1]], SAR = strsplit(SAR, ",")[[1]]),
           by = c("rid", setdiff(names(x), c("rid", "ALT", "AO", "MQM", "SAF", "SAR")))][, rid := NULL]
  }
  num <- c("QUAL", "MQM", "AO", "SAF", "SAR", "DP_tumour", "DP_normal")
  x[, (num) := lapply(.SD, as.numeric), .SDcols = num]
  x[, `:=`(tumour = f$tumour, rel = f$rel, multiallelic_records = multi)]
  x
}
fb <- rbindlist(lapply(seq_len(nrow(files)), read_fb))

# Each condition; NA (missing field) counts as failing, as in bcftools
cond <- c("MQM_ge_40", "GT_normal_00", "GT_tumour_01", "DP_normal_ge_30", "DP_tumour_ge_30", "QUAL_AO_ge_5", "SAF_gt_0", "SAR_gt_0")
fb[, `:=`(MQM_ge_40       = MQM >= 40,
          GT_normal_00    = GT_normal == "0/0",
          GT_tumour_01    = GT_tumour == "0/1",
          DP_normal_ge_30 = DP_normal >= 30,
          DP_tumour_ge_30 = DP_tumour >= 30,
          QUAL_AO_ge_5    = QUAL / AO >= 5,
          SAF_gt_0        = SAF > 0,
          SAR_gt_0        = SAR > 0)]
for (cl in cond) set(fb, which(is.na(fb[[cl]])), cl, FALSE)
fb[, rule := Reduce(`&`, .SD), .SDcols = cond]
fb[, PASS := FILTER %chin% "PASS"]

# ---- 4. Compare with FILTER ----------------------------------------------------
print(fb[, .N, by = .(tumour, FILTER)][order(tumour, -N)])                      # FILTER values present
agree <- fb[, .(records = .N, PASS = sum(PASS), rule_TRUE = sum(rule),
                PASS_and_rule = sum(PASS & rule), PASS_but_rule_FALSE = sum(PASS & !rule),
                notPASS_but_rule_TRUE = sum(!PASS & rule)), by = .(tumour, rel)]
print(agree)

# Which conditions fail among PASS records that break the rule (one record can fail several)
fail_by_cond <- melt(fb[PASS & !rule, c("tumour", cond), with = FALSE], id.vars = "tumour",
                     variable.name = "condition", value.name = "ok")[ok == FALSE, .N, by = .(tumour, condition)]
print(fail_by_cond[order(tumour, -N)])
# Phased or other genotypes among PASS records (e.g. 0|1 would fail the string comparison)
print(fb[(PASS), .N, by = .(tumour, GT_normal, GT_tumour)][order(tumour, -N)])

# ---- 5. Outputs --------------------------------------------------------------
fwrite(agree, file.path(OUT, "freebayes_pass_vs_rule.tsv"), sep = "\t")
fwrite(fail_by_cond, file.path(OUT, "freebayes_pass_failing_conditions.tsv"), sep = "\t")
fwrite(fb[PASS != rule, c("tumour", "CHROM", "POS", "REF", "ALT", "FILTER", "QUAL", "MQM", "AO", "SAF", "SAR",
                          "GT_normal", "DP_normal", "GT_tumour", "DP_tumour", "rule", cond), with = FALSE],
       file.path(OUT, "freebayes_pass_rule_disagreements.tsv"), sep = "\t")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
writeLines(c(paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_freebayes_pass_check"),
             paste0("- Script: ", SCRIPT, "  git: ", git_hash, "  (", R.version.string, ", ", bcf("--version")[1], ")"),
             paste0("- Inputs: ", paste(files$rel, collapse = ", "), " in ", DIR_A),
             paste0("- Records ", nrow(fb), "; PASS ", sum(fb$PASS), "; PASS but rule FALSE ", sum(fb$PASS & !fb$rule),
                    "; not PASS but rule TRUE ", sum(!fb$PASS & fb$rule)),
             paste0("- Outputs: ", OUT), "", old_log), log_file)
message("Done. Outputs: ", OUT)
