#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_identity_crosscheck.R - same data, different method: bcftools gtcheck and KING-robust kinship
# (plink2 --make-king-table) on the SAME masked genotypes and site lists that mouse_wgs_somalier.R made.
# This checks the arithmetic of the somalier / custom IBS0 result; it is NOT independent evidence.
#
# Inputs (from mouse_wgs_somalier.R, run it first):
#   results/mouse_wgs/mouse_wgs_somalier/masked/masked_6samples.bcf   GT:AD, ambiguous tumour genotypes ./.
#   results/mouse_wgs/mouse_wgs_somalier/sites/targets_<combo>.tsv    CHROM POS per panel x thinning
#
# Methods:
#   gtcheck  bcftools gtcheck without -g = all sample pairs; GT only. bcftools >= 1.11: -u GT (+ -E 0 if
#            listed: discordance = number of differing genotypes); 1.4.1: -G 1 (ERR lines). Output parsed
#            from the column names in the version's own header line. Discordance counts ANY genotype
#            difference (het vs hom too), so tumour LOH raises it; compare normals relative to each other.
#   KING     plink2 --make-king-table: KING-robust kinship per pair (NSNP, HETHET, IBS0, KINSHIP).
#            Same individual ~0.5, first degree ~0.25, unrelated ~0. With the few hets of inbred mice the
#            value is unstable (can go strongly negative), like somalier relatedness: use the ranking.
#   peddy    not run: human-only (hg19/hg38 sites for sex, het and ancestry checks).
#
# Per tumour x combo x method: best-matching KI normal (gtcheck: lowest discordance rate; KING: highest
# kinship, tie -> lowest IBS0), next to the custom IBS0 best normal from mouse_wgs_somalier.R.
#
# Install plink2 once (separate env):
#   micromamba create -y -p /common/WORK/pstancl/envs/plink2 -c bioconda -c conda-forge plink2
# Env on ssi: r-mutation (data.table, ggplot2, bcftools). Run: Rscript scripts/mouse_wgs/mouse_wgs_identity_crosscheck.R
# Output: results/mouse_wgs/mouse_wgs_identity_crosscheck/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
SOM_OUT  <- file.path(REPO, "results/mouse_wgs/mouse_wgs_somalier")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_identity_crosscheck")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_identity_crosscheck.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
PLINK2   <- Sys.which("plink2");   if (!nzchar(PLINK2))   PLINK2   <- "/common/WORK/pstancl/envs/plink2/bin/plink2"

TUMOUR_IDS <- c("HCC_1", "HCC_2", "HCC_3")
NORMAL_IDS <- c("KI1", "KI2", "KI3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
started <- Sys.time()
MASKED <- file.path(SOM_OUT, "masked", "masked_6samples.bcf")
if (!file.exists(MASKED)) stop("masked file not found: ", MASKED, " - run mouse_wgs_somalier.R first")
targets <- list.files(file.path(SOM_OUT, "sites"), pattern = "^targets_.*\\.tsv$", full.names = TRUE)
if (!length(targets)) stop("no targets_<combo>.tsv in ", file.path(SOM_OUT, "sites"))
COMBOS <- sub("^targets_(.*)\\.tsv$", "\\1", basename(targets))
dir.create(file.path(OUT, "raw"), recursive = TRUE, showWarnings = FALSE)
NOTES <- character()
note <- function(...) { m <- paste0(...); NOTES <<- c(NOTES, m); message("NOTE: ", m) }
run <- function(cmd) { message("$ ", cmd); st <- system(cmd); if (st != 0) stop("failed (exit ", st, "): ", cmd) }
capture <- function(bin, args = character()) suppressWarnings(system2(bin, args, stdout = TRUE, stderr = TRUE))
HAVE_BCF <- file.exists(BCFTOOLS); HAVE_PLINK <- file.exists(PLINK2)
if (!HAVE_BCF) note("bcftools not found at ", BCFTOOLS, ": gtcheck skipped")
if (!HAVE_PLINK) note("plink2 not found at ", PLINK2, ": KING skipped (install line in the header)")
bcf_version   <- if (HAVE_BCF) capture(BCFTOOLS, "--version")[1] else "bcftools not available"
plink_version <- if (HAVE_PLINK) capture(PLINK2, "--version")[1] else "plink2 not available"
message(bcf_version, " | ", plink_version)
canon <- function(a, b) { s <- c(TUMOUR_IDS, NORMAL_IDS); sw <- match(a, s) > match(b, s); list(ifelse(sw, b, a), ifelse(sw, a, b)) }

# ---- 2. bcftools gtcheck -------------------------------------------------------
gt_usage <- if (HAVE_BCF) capture(BCFTOOLS, c("gtcheck", "-h")) else character()   # -h: never reads stdin
if (HAVE_BCF) writeLines(gt_usage, file.path(OUT, "raw", "bcftools_gtcheck_usage.txt"))
has_opt <- function(o) any(grepl(o, gt_usage, fixed = TRUE))
GT_ARGS <- if (has_opt("-u, --use")) {
  c("-u GT", if (has_opt("-E, --error-probability")) "-E 0" else { note("gtcheck has no -E: discordance is its default error-model score"); NULL },
    if (has_opt("--no-HWE-prob")) "--no-HWE-prob")
} else if (has_opt("-G, --GTs-only")) "-G 1" else { if (HAVE_BCF) note("gtcheck lists neither -u nor -G: skipped"); NULL }
if (!has_opt("-T, --targets-file")) GT_ARGS <- NULL
GT_ARGS

gtcheck_one <- function(cb, tg) {
  raw <- file.path(OUT, "raw", paste0("gtcheck_", cb, ".txt"))
  run(paste(shQuote(BCFTOOLS), "gtcheck", paste(GT_ARGS, collapse = " "), "-T", shQuote(tg), shQuote(MASKED),
            ">", shQuote(raw), "2>", shQuote(paste0(raw, ".stderr"))))
  l <- readLines(raw)
  # newer: "# DCv2 ... [2]query sample [3]genotyped sample [4]discordance ... [6]number of sites compared ..."
  hdr <- grep("^# (DCv2|DC)\\b.*\\[2\\]", l, value = TRUE, perl = TRUE)[1]
  if (!is.na(hdr)) {
    tag <- sub("^# (DCv2|DC)\\b.*", "\\1", hdr, perl = TRUE)
    m <- regmatches(hdr, gregexpr("\\[([0-9]+)\\]([^\\t\\[]+)", hdr, perl = TRUE))[[1]]   # perl: \\t is a tab inside [ ]
    idx <- setNames(as.integer(sub("^\\[([0-9]+)\\].*", "\\1", m)), trimws(sub("[,;:]\\s*$", "", trimws(sub("^\\[[0-9]+\\]", "", m)))))
    need <- c(q = "query sample", g = "genotyped sample", d = "discordance", n = "number of sites compared")
    if (!all(need %in% names(idx))) { note("gtcheck ", cb, ": header lacks ", paste(setdiff(need, names(idx)), collapse = ", "), " - raw kept: ", raw); return(NULL) }
    rows <- grep(paste0("^", tag, "\t"), l, value = TRUE)
    x <- as.data.table(do.call(rbind, strsplit(rows, "\t", fixed = TRUE)))
    d <- data.table(s1 = x[[idx[["query sample"]]]], s2 = x[[idx[["genotyped sample"]]]],
                    discordance = as.numeric(x[[idx[["discordance"]]]]), n_sites = as.numeric(x[[idx[["number of sites compared"]]]]))
  } else if (any(grepl("^# ERR", l))) {                  # 1.4.1: ERR error_rate n_sites sample_i sample_j
    e <- tstrsplit(grep("^ERR\t", l, value = TRUE), "\t", fixed = TRUE)
    d <- data.table(s1 = e[[4]], s2 = e[[5]], n_sites = as.numeric(e[[3]]))
    d[, discordance := round(as.numeric(e[[2]]) * n_sites)]
  } else { note("gtcheck ", cb, ": unknown output format - raw kept: ", raw); return(NULL) }
  d <- d[s1 != s2]
  cs <- canon(d$s1, d$s2); d[, `:=`(s1 = cs[[1]], s2 = cs[[2]])]
  d <- unique(d, by = c("s1", "s2"))
  d[, `:=`(combo = cb, method = "gtcheck", rate = discordance / n_sites)]
}
gtc <- if (length(GT_ARGS) && HAVE_BCF) rbindlist(Map(gtcheck_one, COMBOS, targets)) else data.table()

# ---- 3. KING-robust kinship (plink2) ------------------------------------------------
king_one <- function(cb, tg) {
  pre <- file.path(OUT, "raw", paste0("king_", cb))
  rng <- paste0(pre, ".range.txt")
  t <- fread(tg, header = FALSE, col.names = c("CHROM", "POS"), colClasses = list(character = 1))
  fwrite(t[, .(CHROM, POS, POS, paste0("s", .I))], rng, sep = "\t", col.names = FALSE)
  run(paste(shQuote(PLINK2), "--bcf", shQuote(MASKED), "--allow-extra-chr --chr-set 19 no-xy",
            "--extract bed1", shQuote(rng), "--set-all-var-ids @:#:\\$r:\\$a --make-king-table --out", shQuote(pre),
            ">", shQuote(paste0(pre, ".stdout")), "2>&1"))
  k <- fread(paste0(pre, ".kin0"))
  setnames(k, sub("^#", "", names(k)))
  need <- c("IID1", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP")
  if (!all(need %in% names(k))) { note("KING ", cb, ": .kin0 lacks ", paste(setdiff(need, names(k)), collapse = ", ")); return(NULL) }
  cs <- canon(k$IID1, k$IID2)
  data.table(combo = cb, method = "KING", s1 = cs[[1]], s2 = cs[[2]], n_sites = k$NSNP, hethet = k$HETHET,
             ibs0 = k$IBS0, kinship = k$KINSHIP, rate = k$IBS0 / k$NSNP)
}
king <- if (HAVE_PLINK) rbindlist(Map(king_one, COMBOS, targets)) else data.table()

# ---- 4. Best normal per tumour, next to the custom IBS0 result -----------------------
pairs <- rbindlist(list(gtc, king), fill = TRUE)
tn <- pairs[(s1 %chin% TUMOUR_IDS) != (s2 %chin% TUMOUR_IDS)]
tn[, `:=`(tumour = fifelse(s1 %chin% TUMOUR_IDS, s1, s2), normal = fifelse(s1 %chin% TUMOUR_IDS, s2, s1))]
# gtcheck: lowest discordance rate; KING: highest kinship, then lowest IBS0 rate
tn[, sort1 := fifelse(method == "KING", -kinship, rate)]
tn[, sort2 := fifelse(method == "KING", rate, 0)]
setorder(tn, combo, method, tumour, sort1, sort2)
best <- tn[, .(best_normal = normal[1], second_normal = normal[2],
               best_value = fifelse(method[1] == "KING", kinship[1], rate[1]),
               second_value = fifelse(method[1] == "KING", kinship[2], rate[2]), n_sites_best = n_sites[1]),
           by = .(combo, method, tumour)]
best[, labelled_normal := KI_NORMALS[tumour]]
# custom IBS0 best normal from mouse_wgs_somalier.R (robustness.tsv)
rob_f <- file.path(SOM_OUT, "robustness.tsv")
if (file.exists(rob_f)) {
  rob <- fread(rob_f)[method == "custom", .(combo, tumour, custom_ibs0_best = best_normal, custom_ambiguous = ambiguous)]
  best <- rob[best, on = .(combo, tumour)]
  best[, agrees_with_custom := best_normal == custom_ibs0_best]
} else note("robustness.tsv not found in ", SOM_OUT, ": no comparison with the custom IBS0")
best[, combo := factor(combo, unique(c(grep("^A_", COMBOS, value = TRUE), grep("^B_", COMBOS, value = TRUE), grep("^C_", COMBOS, value = TRUE))))]
setorder(best, method, tumour, combo)
print(dcast(best, method + tumour ~ combo, value.var = "best_normal"))
if ("agrees_with_custom" %in% names(best)) print(best[, .(combos = .N, agree = sum(agrees_with_custom)), by = method])

# ---- 5. Plot, outputs, RUN_LOG.md ----------------------------------------------
if (nrow(best)) {
  p <- ggplot(best, aes(combo, tumour, fill = best_normal)) + geom_tile(colour = "white") +
    geom_text(aes(label = best_normal), size = 3) + facet_wrap(~ method, ncol = 1) +
    labs(x = "panel_thinning", y = NULL, fill = "best normal",
         title = "Best-matching KI normal per tumour: gtcheck and KING on the masked genotypes",
         subtitle = "Same genotypes and sites as mouse_wgs_somalier.R - a check of the arithmetic, not independent evidence") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT, "best_normal.png"), p, width = 12, height = 2 + 1.8 * uniqueN(best$method), dpi = 150)
  ggsave(file.path(OUT, "best_normal.pdf"), p, width = 12, height = 2 + 1.8 * uniqueN(best$method))
}
fwrite(pairs, file.path(OUT, "pairwise_gtcheck_king.tsv"), sep = "\t")
fwrite(best, file.path(OUT, "best_normal.tsv"), sep = "\t")
git_hash <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)),
                     error = function(e) "unknown")
if (!length(git_hash)) git_hash <- "unknown"
md_table <- function(x) { x <- as.data.table(x); if (!nrow(x)) return("(none)")
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))) }
prim <- best[combo == "A_20kb", .(method, tumour, best_normal, best_value = signif(best_value, 4), second_normal,
                                   second_value = signif(second_value, 4), n_sites_best,
                                   custom_ibs0_best = if ("custom_ibs0_best" %in% names(best)) custom_ibs0_best else NA)]
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_identity_crosscheck"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Tools:      ", bcf_version, "; ", plink_version),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Inputs:     ", MASKED, "; ", length(targets), " site lists in ", file.path(SOM_OUT, "sites")),
  paste0("- gtcheck:    ", if (length(GT_ARGS)) paste("bcftools gtcheck", paste(GT_ARGS, collapse = " "), "-T targets_<combo>.tsv masked") else "not run"),
  paste0("- KING:       ", if (HAVE_PLINK) "plink2 --bcf masked --extract range <combo sites> --make-king-table" else "not run"),
  paste0("- Outputs:    ", OUT),
  "", "### Primary combination A_20kb", "", md_table(prim),
  "", "### Notes",
  "- Same masked genotypes and sites as mouse_wgs_somalier.R: a check of the arithmetic, not independent evidence.",
  "- gtcheck discordance counts any genotype difference (het vs hom included), so tumour LOH raises it; KING kinship is unstable with few hets (inbred mice): use rankings.",
  "- peddy not run: human-only site sets.",
  if (length(NOTES)) paste0("- NOTE: ", NOTES) else "- no notes",
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
