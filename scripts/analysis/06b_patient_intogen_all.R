#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), human WGS RJALS | Created: 2026-11-25
#
# 06b_patient_intogen_all.R - harmful patient mutations in ALL IntOGen driver genes (any tumour type) and
# the evidence to judge each one.
#
# Logic (2026-11-25):
#   hits       tumourevo mutations of the tumour (RJALS_RJALS_Tm in the driver.rds; PASS only if a FILTER column
#              exists), VEP IMPACT HIGH or MODERATE, gene in IntOGen Compendium_Cancer_Genes.tsv (any TUMOUR_TYPE).
#              Gene columns: ROLE, in_HCC, n_tumour_types, tumour_types, in_CGC; mutation: is_driver (tumourevo,
#              HCC genes only), driver_label.
#   check 1    known / hotspot: VEP CLIN_SIG and Existing_variation (COSV / rs IDs) if present; cancerhotspots v2
#              match on gene + amino-acid position if HOTSPOTS is set; links to cancerhotspots, OncoKB,
#              ClinVar, COSMIC for the manual check
#   check 2    VAF / clonality: VAF, depth and any CCF / cluster / copy-number columns tumourevo has
#   check 3    germline: same chrom/pos/ref/alt in the normal's table (RJALS_RJALS_N) of the rds; gnomAD / MAX_AF
#              > GERMLINE_AF -> germline_like
#
# Env on ssi: r-mutation. Output: results/wgs/analysis/drivers/ (next to 06_drivers.R)

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
DRV_RDS  <- file.path(REPO, "results/wgs/tumourevo/RJALS/driver_annotation/annotate_driver/RJALS/RJALS/RJALS_RJALS_Tm/RJALS_RJALS_RJALS_RJALS_Tm_driver.rds")
TUMOUR   <- "RJALS_RJALS_Tm"
NORMAL   <- "RJALS_RJALS_N"
INTOGEN  <- "/common/WORK/pstancl/references/intogen/Compendium_Cancer_Genes.tsv"
# cancerhotspots v2 (optional): download https://www.cancerhotspots.org/files/hotspots_v2.xls to this path; NA = skip
HOTSPOTS <- NA
GERMLINE_AF <- 0.001
OUT      <- file.path(REPO, "results/wgs/analysis/drivers")
SCRIPT   <- file.path(REPO, "scripts/analysis/06b_patient_intogen_all.R")

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
for (f in c(DRV_RDS, INTOGEN)) if (!file.exists(f)) stop("Not found: ", f)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
pick1 <- function(cands, nm) intersect(cands, nm)[1]

# ---- 2. Patient mutations and IntOGen genes --------------------------------------
x  <- readRDS(DRV_RDS)
pm <- as.data.table(x[[TUMOUR]]$mutations)
if ("FILTER" %in% names(pm)) pm <- pm[FILTER == "PASS"]
ig <- fread(INTOGEN)
tt <- pick1(c("TUMOUR_TYPE", "CANCER_TYPE"), names(ig))
igg <- ig[, .(ROLE = paste(sort(unique(ROLE)), collapse = ";"), in_HCC = any(get(tt) == "HCC"),
              n_tumour_types = uniqueN(get(tt)), tumour_types = paste(sort(unique(get(tt))), collapse = ";"),
              in_CGC = any(as.character(CGC_GENE) %chin% c("True", "TRUE"))), by = .(SYMBOL)]
message(nrow(igg), " IntOGen genes (all tumour types)")

hits <- pm[IMPACT %in% c("HIGH", "MODERATE") & SYMBOL %chin% igg$SYMBOL]
hits <- igg[hits, on = "SYMBOL"]
c_chr <- pick1(c("chr", "CHROM", "seqnames"), names(hits)); c_pos <- pick1(c("from", "POS", "start"), names(hits))
c_ref <- pick1(c("ref", "REF"), names(hits)); c_alt <- pick1(c("alt", "ALT"), names(hits))
hits[, aa := suppressWarnings(as.integer(sub("^.*p\\.[A-Za-z*]+?([0-9]+).*$", "\\1", as.character(HGVSp))))]

# ---- 3. Check 1: known / hotspot ---------------------------------------------------
known_cols <- grep("^CLIN_SIG$|^Existing_variation$|^SOMATIC$|^PHENO$|^SIFT$|^PolyPhen$", names(hits), value = TRUE)
hits[, hotspot := NA]
if (!is.na(HOTSPOTS) && file.exists(HOTSPOTS)) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl needed for the hotspots .xls")
  hs <- rbindlist(lapply(readxl::excel_sheets(HOTSPOTS), function(s) {
    d <- as.data.table(readxl::read_excel(HOTSPOTS, sheet = s))
    g <- pick1(c("Hugo_Symbol", "Hugo Symbol", "Gene"), names(d)); p <- pick1(c("Amino_Acid_Position", "Amino Acid Position", "Residue"), names(d))
    if (is.na(g) || is.na(p)) return(NULL)
    d[, .(SYMBOL = as.character(get(g)), aa = suppressWarnings(as.integer(gsub("[^0-9].*$", "", gsub("^[A-Za-z]", "", as.character(get(p)))))))]
  }), fill = TRUE)
  hits[, hotspot := paste(SYMBOL, aa) %chin% hs[, paste(SYMBOL, aa)]]
} else message("HOTSPOTS not set: hotspot column NA (use the links)")
hits[, `:=`(p_short = sub("^.*:", "", as.character(HGVSp)))]
hits[, `:=`(link_cancerhotspots = paste0("https://www.cancerhotspots.org/#/gene/", SYMBOL),
            link_OncoKB  = paste0("https://www.oncokb.org/gene/", SYMBOL, "/", sub("^p\\.", "", p_short)),
            link_ClinVar = paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=", SYMBOL, "%5Bgene%5D+", sub("^p\\.", "", p_short)),
            link_COSMIC  = paste0("https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=", SYMBOL))]

# ---- 4. Check 2: VAF / clonality ---------------------------------------------------
clonal_cols <- grep("VAF|^DP$|^NV$|^NR$|CCF|cluster|karyotype|multiplicity|^CN|minor|Major", names(hits), value = TRUE)

# ---- 5. Check 3: germline -------------------------------------------------------
af_cols <- grep("^MAX_AF$|gnomAD.*AF$|^AF$", names(hits), value = TRUE)
hits[, max_pop_AF := if (length(af_cols)) do.call(pmax, c(lapply(.SD, function(v) suppressWarnings(as.numeric(as.character(v)))), na.rm = TRUE)) else NA_real_,
     .SDcols = af_cols]
hits[, germline_like_popAF := !is.na(max_pop_AF) & max_pop_AF > GERMLINE_AF]
hits[, in_normal_table := NA]
if (NORMAL %in% names(x) && !is.null(x[[NORMAL]]$mutations) && !anyNA(c(c_chr, c_pos, c_ref, c_alt))) {
  pn <- as.data.table(x[[NORMAL]]$mutations)
  if (all(c(c_chr, c_pos, c_ref, c_alt) %in% names(pn)))
    hits[, in_normal_table := paste(get(c_chr), get(c_pos), get(c_ref), get(c_alt)) %chin%
           pn[, paste(get(c_chr), get(c_pos), get(c_ref), get(c_alt))]]
}

# ---- 6. Output ------------------------------------------------------------------
main <- unique(c("SYMBOL", c_chr, c_pos, c_ref, c_alt, "in_HCC", "ROLE", "tumour_types", "n_tumour_types", "in_CGC",
                 "IMPACT", "Consequence", "p_short", "aa", "is_driver", "driver_label", "hotspot", known_cols,
                 clonal_cols, "max_pop_AF", "germline_like_popAF", "in_normal_table",
                 "link_cancerhotspots", "link_OncoKB", "link_ClinVar", "link_COSMIC"))
main <- intersect(main[!is.na(main)], names(hits))
res <- hits[, ..main]
setorder(res, -in_HCC, SYMBOL)
for (i in seq_len(nrow(res))) { cat("\n====", res$SYMBOL[i], res$p_short[i], "====\n"); print(t(res[i, !grepl("^link_", names(res)), with = FALSE]), quote = FALSE) }
fwrite(res, file.path(OUT, "patient_intogen_all_harmful_evidence.tsv"), sep = "\t")

log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
writeLines(c(paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - 06b_patient_intogen_all"),
             paste0("- Script: ", SCRIPT, "  git: ", git_hash, "  (", R.version.string, ")"),
             paste0("- Inputs: ", DRV_RDS, " (", TUMOUR, ", normal ", NORMAL, "); ", INTOGEN, "; hotspots ", HOTSPOTS),
             paste0("- Parameters: VEP HIGH/MODERATE; germline_like if pop AF > ", GERMLINE_AF),
             paste0("- Result: ", nrow(res), " harmful mutations in ", uniqueN(res$SYMBOL), " IntOGen genes (", sum(res$in_HCC), " in HCC genes)"),
             paste0("- Output: ", file.path(OUT, "patient_intogen_all_harmful_evidence.tsv")), "", old_log), log_file)
message("Done. Output: ", file.path(OUT, "patient_intogen_all_harmful_evidence.tsv"))
