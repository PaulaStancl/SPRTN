# ---------------------------------------------------------------------------
# 04_summary_table.R - the "results compared" table: sarek vs oncoanalyser vs shared.
# Run after 03_compare_callers.R (SNVs / indels come from its comparison/ tables).
# Run:  Rscript 04_summary_table.R   (or line by line in R)
#
# Rows and their sources:
#   PASS SNVs / indels   comparison/overlap/mutation_callers.csv (03): MNVs split, indels
#                        normalised if 06 ran. sarek = Mutect2, Strelka2 (+ MuSE, SNVs, if run);
#                        oncoanalyser = SAGE (PURPLE's final VCF).
#                        sarek consensus: SNVs >= 2 of Mutect2, Strelka2, MuSE (indels: Mutect2 AND
#                        Strelka2 - MuSE calls no indels). Shared = sarek consensus AND SAGE.
#   Structural variants  PASS records of Manta (somatic) and ESVEE; shared = Manta records with a
#                        breakpoint within SV_TOL bp of an ESVEE breakpoint (same chromosome).
#   CN segments          ASCAT metrics n_segs / PURPLE qc CopyNumberSegments (segmentation
#                        granularity differs - not a measure of agreement)
#   Purity / ploidy, WGD, LOH   ASCAT purityploidy + metrics / PURPLE purity.tsv + qc
#   Driver events        PURPLE somatic driver catalogue (sarek calls no drivers); shared = PURPLE's
#                        reported driver mutations that are also PASS in Mutect2 or Strelka2
# Output: <OUT>/summary/summary_table.csv, summary_table.md (with footnotes), drivers.csv
# ---------------------------------------------------------------------------
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
od  <- file.path(OUT, "summary"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
cmp <- file.path(OUT, "comparison")
`%||%` <- function(a, b) if (is.null(a)) b else a   # (base R has it only from 4.4)
SV_TOL <- 100L                                        # bp, breakpoint matching tolerance
fmt <- function(x) if (length(x) == 0 || is.na(x)) "NA" else format(x, big.mark = ",", trim = TRUE)
pct <- function(a, b) if (is.na(b) || b == 0) "NA" else sprintf("%.0f%%", 100 * a / b)
kv  <- function(f) {                                  # two-column key/value file -> named vector
  if (is.na(f) || !file.exists(f)) return(NULL)
  d <- fread(f, header = FALSE, sep = "\t", fill = TRUE); setNames(as.character(d[[2]]), d[[1]])
}
tr1 <- function(f) {                                  # one-row table with header -> named vector
  if (is.na(f) || !file.exists(f)) return(NULL)
  d <- fread(f); setNames(as.character(unlist(d[1])), names(d))
}
rows <- list()

# ---- 1. SNVs and indels (from 03) ------------------------------------------------------
mc <- file.path(cmp, "overlap", "mutation_callers.csv")
if (!file.exists(mc)) stop("missing ", mc, " - run 03_compare_callers.R first", call. = FALSE)
mut <- fread(mc)
has <- function(cl) grepl(paste0("(^|\\+)", cl, "(\\+|$)"), mut$CALLERS)
mut[, `:=`(m2 = has("mutect2"), st = has("strelka"), mu = has("muse"), sg = has("sage"))]
# sarek consensus: SNVs - at least 2 of Mutect2, Strelka2, MuSE (if MuSE ran; else Mutect2 AND
# Strelka2); indels - Mutect2 AND Strelka2 (MuSE calls no indels). Shared = sarek consensus AND SAGE.
have_muse <- any(mut$mu)
for (ty in c("SNV", "INDEL")) {
  x <- mut[TYPE == ty]
  use_muse <- ty == "SNV" && have_muse
  x[, cons := if (use_muse) (m2 + st + mu) >= 2 else m2 & st]
  cons_lab <- if (use_muse) ">= 2 of 3" else "Mutect2+Strelka2"
  n_cons <- x[(cons), .N]; sage <- x[(sg), .N]; shared <- x[cons & sg, .N]
  sarek <- paste0("Mutect2 ", fmt(x[(m2), .N]), " · Strelka2 ", fmt(x[(st), .N]),
                  if (use_muse) paste0(" · MuSE ", fmt(x[(mu), .N])),
                  " · ", cons_lab, " ", fmt(n_cons),
                  if (use_muse) paste0(" (all 3: ", fmt(x[m2 & st & mu, .N]), ")"))
  rows[[length(rows) + 1]] <- data.table(
    Result = paste("PASS", if (ty == "SNV") "SNVs" else "indels"), sarek = sarek,
    oncoanalyser = paste("SAGE", fmt(sage)),
    Shared = sprintf("%s (%s of SAGE, %s of sarek %s)", fmt(shared), pct(shared, sage), pct(shared, n_cons), cons_lab))
}

# ---- 2. structural variants: Manta vs ESVEE ----------------------------------------------
manta_f <- find_one(file.path(SAREK, "variant_calling", "manta", PAIR), "\\.manta\\.somatic_sv\\.vcf\\.gz$", required = FALSE)
esvee_f <- find_one(file.path(ONCO, "esvee"), "\\.esvee\\.somatic\\.vcf\\.gz$", required = FALSE)
breakends <- function(sv) {                           # POS, plus END for intrachromosomal SVs (REC = record number)
  b <- sv[, .(REC, CHROM, BP = POS)]
  if ("info_END" %in% names(sv)) {
    e <- sv[!SVTYPE %chin% c("BND", "INS") & !is.na(info_END), .(REC, CHROM, BP = as.integer(info_END))]
    b <- rbind(b, e)
  }
  b[CHROM %chin% STD_CHR]
}
if (!is.na(manta_f) && !is.na(esvee_f)) {
  ma <- read_sv_vcf(manta_f)[FILTER == "PASS" & CHROM %chin% STD_CHR][, REC := .I]
  es <- read_sv_vcf(esvee_f)[FILTER == "PASS" & CHROM %chin% STD_CHR][, REC := .I]
  bm <- breakends(ma); be <- breakends(es)[, E := BP]
  setkey(be, CHROM, BP)
  hit <- be[bm, on = .(CHROM, BP), roll = "nearest"][!is.na(E) & abs(E - BP) <= SV_TOL, unique(i.REC)]
  rows[[length(rows) + 1]] <- data.table(
    Result = "Structural variants (PASS records)", sarek = paste("Manta", fmt(nrow(ma))),
    oncoanalyser = paste("ESVEE", fmt(nrow(es))),
    Shared = sprintf("%s Manta records (%s) with a breakpoint within %d bp of ESVEE", fmt(length(hit)),
                     pct(length(hit), nrow(ma)), SV_TOL))
} else message("SV row skipped - Manta or ESVEE VCF not found")

# ---- 3. copy number, purity / ploidy, WGD, LOH --------------------------------------------
asc <- file.path(SAREK, "variant_calling", "ascat", PAIR)
am  <- tr1(find_one(asc, "\\.metrics\\.txt$", required = FALSE))
app <- tr1(find_one(asc, "\\.purityploidy\\.txt$", required = FALSE))
pq  <- kv(find_one(file.path(ONCO, "purple"), "\\.purple\\.qc$", required = FALSE))
pp  <- tr1(find_one(file.path(ONCO, "purple"), "\\.purple\\.purity\\.tsv$", required = FALSE))
num <- function(x, d = 2) if (is.null(x) || is.na(x)) "NA" else formatC(as.numeric(x), format = "f", digits = d)
rows[[length(rows) + 1]] <- data.table(Result = "Copy-number segments",
  sarek = paste("ASCAT", am[["n_segs"]] %||% "NA"), oncoanalyser = paste("PURPLE", pq[["CopyNumberSegments"]] %||% "NA"), Shared = "—")
rows[[length(rows) + 1]] <- data.table(Result = "Purity / ploidy",
  sarek = sprintf("%s%% / %s", num(100 * as.numeric(app[["AberrantCellFraction"]]), 0), num(app[["Ploidy"]])),
  oncoanalyser = sprintf("%s%% / %s", num(100 * as.numeric(pp[["purity"]]), 0), num(pp[["ploidy"]])), Shared = "—")
wgd_a <- if (!is.null(am[["WGD"]])) c("0" = "no", "1" = "yes")[am[["WGD"]]] else "NA"
wgd_p <- if (!is.null(pp[["wholeGenomeDuplication"]])) c("false" = "no", "true" = "yes")[tolower(pp[["wholeGenomeDuplication"]])] else "NA"
rows[[length(rows) + 1]] <- data.table(Result = "Whole-genome doubling", sarek = wgd_a, oncoanalyser = wgd_p,
  Shared = if (identical(unname(wgd_a), unname(wgd_p))) "agree" else "differ")
rows[[length(rows) + 1]] <- data.table(Result = "LOH (% of genome)",
  sarek = paste0(num(100 * as.numeric(am[["LOH"]]), 1), "%"),
  oncoanalyser = paste0(num(100 * as.numeric(pq[["LohPercent"]]), 1), "%"), Shared = "—")

# ---- 4. driver events (PURPLE) ----------------------------------------------------------------
dc_f <- find_one(file.path(ONCO, "purple"), "\\.purple\\.driver\\.catalog\\.somatic\\.tsv$", required = FALSE)
if (!is.na(dc_f)) {
  dc <- fread(dc_f)
  fwrite(dc, file.path(od, "drivers.csv"))
  drv <- if (nrow(dc)) dc[, paste0(gene, " (", tolower(driver), ")")] else character(0)
  # PURPLE's reported driver mutations (INFO/REPORTED in its VCF) also PASS in a sarek caller
  pa <- file.path(cmp, "pass_calls", "snv_indel_pass_atomized.csv")
  shared_drv <- "NA"
  if (file.exists(pa)) {
    pr <- fread(pa, select = intersect(c("CALLER", "CHROM", "POS", "REF", "ALT1", "info_REPORTED"), names(fread(pa, nrows = 0))))
    if ("info_REPORTED" %in% names(pr)) {
      rep <- unique(pr[CALLER == "sage" & as.character(info_REPORTED) %chin% c("TRUE", "true", "1"),
                       .(KEY = paste(CHROM, POS, REF, ALT1, sep = ":"))])
      ok  <- rep[KEY %chin% mut[m2 | st, KEY], .N]
      shared_drv <- sprintf("%d of %d driver mutations also PASS in Mutect2 or Strelka2", ok, nrow(rep))
    }
  }
  rows[[length(rows) + 1]] <- data.table(Result = "Driver events",
    sarek = "— (no driver calling)",
    oncoanalyser = sprintf("%d: %s", nrow(dc), if (length(drv)) paste(drv, collapse = ", ") else "none"),
    Shared = shared_drv)
} else message("driver row skipped - PURPLE driver catalogue not found")

# ---- write ---------------------------------------------------------------------------------
tab <- rbindlist(rows)
print(tab); fwrite(tab, file.path(od, "summary_table.csv"))
md <- c("| Result | sarek | oncoanalyser | Shared |", "|---|---|---|---|",
        tab[, sprintf("| %s | %s | %s | %s |", Result, sarek, oncoanalyser, Shared)], "",
        "*sarek consensus: SNVs PASS in at least 2 of Mutect2, Strelka2 and MuSE (MuSE: PASS + Tier1-5); indels PASS in Mutect2 and Strelka2 (MuSE calls no indels). Shared: sarek consensus and PASS in SAGE. Matched on CHROM:POS:REF:ALT, MNVs split into SNVs, indels normalised with bcftools norm (06_normalize_vcfs.sh).*",
        sprintf("*Shared SVs: Manta PASS records with a breakpoint within %d bp of an ESVEE PASS breakpoint.*", SV_TOL),
        "*Copy-number segment counts reflect each tool's segmentation, not agreement.*",
        "*Driver events: PURPLE somatic driver catalogue; sarek does no driver calling.*")
writeLines(md, file.path(od, "summary_table.md"))
message("done: ", od)
