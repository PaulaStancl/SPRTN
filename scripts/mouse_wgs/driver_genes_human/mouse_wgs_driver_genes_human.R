#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_driver_genes_human.R - harmful mutations in mouse orthologs of human cancer driver genes.
#
# Logic (agreed 2026-11-25):
#   gene lists    two separate analyses (own tables and plots, suffix _IntOGen_HCC / _CGC):
#                   IntOGen_HCC  Compendium_Cancer_Genes.tsv, rows with TUMOUR_TYPE == "HCC" (106 genes; ROLE,
#                                IntOGen's CGC_GENE flag as in_CGC, number of HCC cohorts kept)
#                   CGC          COSMIC Cancer Gene Census, all genes (Tier 1 + 2; TIER and ROLE kept);
#                                old cancer_gene_census.csv or new Cosmic_CancerGeneCensus_*.tsv
#   patient_driver  TRUE if the human gene has is_driver == TRUE in the RJALS tumourevo driver.rds
#   orthologs     Ensembl Compara human -> mouse table from cross_species/h2m 00_resources.R
#                 (ensembl<R>/orthologs_hsapiens_mmusculus.tsv); all ortholog types, homology type reported;
#                 matched on the mouse Ensembl gene ID (VEP CSQ field Gene)
#   harmful       a VEP CSQ entry on that mouse gene with IMPACT HIGH or MODERATE (any transcript); per
#                 variant x gene the most severe entry is reported (canonical transcript preferred); SIFT shown
#   variants      SNVs (MNVs split) and indels PASS in any of the 9 sets of mouse_wgs_snv_indel_settings.R
#                 (its settings.rds); one TRUE/FALSE column per set. CSQ is read from the same files that
#                 set used (single-caller files of setting A, the B and C files); an MNV's CSQ is given to
#                 each of its SNVs (FROM_MNV = TRUE).
#
# Needs: mouse_wgs_snv_indel_settings.R run first (results/mouse_wgs/mouse_wgs_snv_indel_settings/settings.rds).
#
# Blocks: 0 EDIT | 1 setup | 2 gene lists, patient drivers, orthologs | 3 PASS variants per set
#         | 4 CSQ of the orthologs | 5 harmful hits (all genes) | 6 per list: tables + plot | 7 RUN_LOG.md
# Env on ssi: r-mutation. Output: results/mouse_wgs/driver_genes_human/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO      <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
SETTINGS  <- file.path(REPO, "results/mouse_wgs/mouse_wgs_snv_indel_settings/settings.rds")
OUT       <- file.path(REPO, "results/mouse_wgs/driver_genes_human")
SCRIPT    <- file.path(REPO, "scripts/mouse_wgs/driver_genes_human/mouse_wgs_driver_genes_human.R")
# tumourevo's IntOGen compendium; find it with:
#   find /common/WORK/pstancl ~/.nextflow -name Compendium_Cancer_Genes.tsv 2>/dev/null
INTOGEN   <- "/common/WORK/pstancl/references/intogen/Compendium_Cancer_Genes.tsv"
CANCER    <- "HCC"
# COSMIC Cancer Gene Census (download needs a COSMIC login): old cancer_gene_census.csv or new
# Cosmic_CancerGeneCensus_v*_GRCh38.tsv(.gz). NA or file missing = the CGC analysis is skipped (2026-11-25: not on ssi yet)
CGC       <- NA
# Patient (human RJALS) tumourevo driver annotation: is_driver per mutation
PATIENT_RDS <- file.path(REPO, "results/wgs/tumourevo/RJALS/driver_annotation/annotate_driver/RJALS/RJALS/RJALS_RJALS_Tm/RJALS_RJALS_RJALS_RJALS_Tm_driver.rds")
ORTHOLOGS <- "/common/WORK/pstancl/references/cross_species/ensembl110/orthologs_hsapiens_mmusculus.tsv"
IMPACTS   <- c("HIGH", "MODERATE")
ENV_BIN   <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS  <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
HAVE_CGC <- !is.na(CGC) && file.exists(CGC)
if (!HAVE_CGC) message("CGC file not set or not found (", CGC, ") - CGC analysis skipped, IntOGen HCC only")
for (f in c(SETTINGS, INTOGEN, ORTHOLOGS, PATIENT_RDS)) if (!file.exists(f)) stop("Not found: ", f,
  "\n  settings.rds: run mouse_wgs_snv_indel_settings.R; orthologs: run scripts/cross_species/h2m 00_resources.R;",
  "\n  IntOGen / CGC / patient rds: set the path in block 0")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
save_plot <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150, limitsize = FALSE)
}

# ---- 2. Gene lists, patient drivers, orthologs -----------------------------------
ig <- fread(INTOGEN)
tt <- intersect(c("TUMOUR_TYPE", "CANCER_TYPE"), names(ig))[1]          # tumourevo's copy has TUMOUR_TYPE
if (is.na(tt) || !"SYMBOL" %in% names(ig)) stop("IntOGen: need SYMBOL and TUMOUR_TYPE (or CANCER_TYPE) in ", INTOGEN)
intogen <- ig[get(tt) == CANCER, .(ROLE = if ("ROLE" %in% names(ig)) paste(sort(unique(na.omit(ROLE))), collapse = ";") else NA_character_,
                                   TIER = NA_character_,
                                   in_CGC = if ("CGC_GENE" %in% names(ig)) any(as.character(CGC_GENE) %chin% c("True", "TRUE")) else NA,
                                   n_cohorts = uniqueN(COHORT)), by = .(symbol_h = SYMBOL)]
LISTS <- list(IntOGen_HCC = intogen)
if (HAVE_CGC) {
  cg <- fread(CGC)
  pick <- function(cands) { k <- intersect(cands, names(cg)); if (!length(k)) stop("CGC: none of ", paste(cands, collapse = "/"), " in ", CGC); k[1] }
  c_sym <- pick(c("GENE_SYMBOL", "Gene Symbol")); c_tier <- pick(c("TIER", "Tier")); c_role <- pick(c("ROLE_IN_CANCER", "Role in Cancer"))
  LISTS$CGC <- unique(cg[, .(symbol_h = get(c_sym), ROLE = as.character(get(c_role)), TIER = as.character(get(c_tier)))], by = "symbol_h")
}
for (L in names(LISTS)) message(L, ": ", nrow(LISTS[[L]]), " genes")

# Patient: human genes with is_driver == TRUE in the tumourevo driver.rds (mutations table of the sample)
px <- readRDS(PATIENT_RDS)
pm <- NULL
for (el in px) if (is.list(el) && !is.null(el$mutations)) pm <- as.data.table(el$mutations)
if (is.null(pm) || !all(c("is_driver", "SYMBOL") %in% names(pm))) stop("PATIENT_RDS: no mutations table with is_driver and SYMBOL")
patient_genes <- unique(as.character(pm[is_driver %in% TRUE, SYMBOL]))
message("patient driver genes (is_driver): ", paste(patient_genes, collapse = ", "))

# Orthologs of every gene in either list (CSQ is read once for all; lists are separated in block 6)
all_genes <- unique(unlist(lapply(LISTS, `[[`, "symbol_h")))
orth <- fread(ORTHOLOGS, colClasses = "character", na.strings = c("", "NA"))
d_orth <- orth[!is.na(gene_id_m) & symbol_h %chin% all_genes, .(symbol_h, gene_id_h, gene_id_m, symbol_m, homology_type, perc_id_h_to_m)]
print(d_orth[, .N, by = homology_type])

# ---- 3. PASS variants per set (settings.rds) -----------------------------------
st <- readRDS(SETTINGS)
pass <- as.data.table(st$pass)[, .(tumour, SET = as.character(SET), CLASS, CHROM, POS, REF, ALT)]
SET_LEVELS <- st$params$SET_LEVELS
files <- as.data.table(st$files)                       # tumour, src, path of the files the sets came from
print(pass[, .N, by = .(tumour, SET)][order(tumour, match(SET, SET_LEVELS))])

# ---- 4. VEP CSQ of the driver orthologs, from the source files ------------------
ids_file <- tempfile(fileext = ".txt"); writeLines(unique(d_orth$gene_id_m), ids_file)
read_csq <- function(i) {
  f <- files[i]
  hdr <- system2(BCFTOOLS, c("view", "-h", shQuote(f$path)), stdout = TRUE)
  h <- grep("^##INFO=<ID=CSQ,", hdr, value = TRUE)
  if (!length(h)) { message(f$rel, ": no CSQ in header"); return(NULL) }
  fields <- strsplit(sub("\".*$", "", sub("^.*Format: ", "", h)), "|", fixed = TRUE)[[1]]
  # only records whose CSQ mentions a driver ortholog gene ID (grep -F keeps it fast)
  cmd <- paste(shQuote(BCFTOOLS), "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CSQ\\n"), shQuote(f$path),
               "| grep -F -f", shQuote(ids_file), "|| true")
  x <- tryCatch(fread(cmd = cmd, sep = "\t", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", "CSQ"),
                      colClasses = list(character = c("CHROM", "REF", "ALT", "CSQ"), integer = "POS"), quote = ""),
                error = function(e) data.table())
  if (!nrow(x)) { message(f$rel, ": no records in driver orthologs"); return(NULL) }
  x <- x[, .(entry = unlist(strsplit(CSQ, ",", fixed = TRUE))), by = .(CHROM, POS, REF, ALT)]
  m <- tstrsplit(x$entry, "|", fixed = TRUE, fill = NA_character_)
  length(m) <- length(fields)                                     # pad if trailing fields are empty everywhere
  names(m) <- fields
  keep <- intersect(c("Consequence", "IMPACT", "SYMBOL", "Gene", "Feature", "BIOTYPE", "HGVSc", "HGVSp",
                      "Protein_position", "Amino_acids", "SIFT", "CANONICAL", "EXON"), fields)
  x <- cbind(x[, .(CHROM, POS, REF, ALT)], as.data.table(lapply(m[keep], function(v) if (is.null(v)) rep(NA_character_, nrow(x)) else v)))
  for (k in setdiff(c("CANONICAL", "SIFT", "HGVSp", "HGVSc"), names(x))) x[, (k) := NA_character_]
  x[, `:=`(tumour = f$tumour, src = f$src)]
  x[Gene %chin% d_orth$gene_id_m]
}
csq <- rbindlist(lapply(seq_len(nrow(files)), read_csq), use.names = TRUE, fill = TRUE)
if (!nrow(csq)) stop("No CSQ entries on driver orthologs in any file - check that CSQ Gene holds mouse Ensembl gene IDs")
csq <- unique(csq)
print(csq[, .N, by = .(IMPACT)])

# MNVs: give the CSQ to each SNV of the MNV (as the sets split MNVs)
csq[, MNV := nchar(REF) == nchar(ALT) & nchar(REF) > 1L & REF %like% "^[ACGT]+$" & ALT %like% "^[ACGT]+$"]
mn <- csq[(MNV)]
if (nrow(mn)) {
  mn <- mn[rep(seq_len(.N), nchar(REF))][, k := sequence(nchar(mn$REF))]
  mn[, `:=`(POS = POS + k - 1L, REF = substring(REF, k, k), ALT = substring(ALT, k, k))]
  mn <- mn[REF != ALT][, k := NULL]
}
csq <- rbind(csq[!(MNV)][, FROM_MNV := FALSE], mn[, FROM_MNV := TRUE], fill = TRUE)[, MNV := NULL]

# ---- 5. Harmful hits in PASS variants -----------------------------------------
IMPACT_RANK <- c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)
h <- csq[IMPACT %chin% IMPACTS]
h[, `:=`(rank = IMPACT_RANK[IMPACT], canon = fifelse(CANONICAL %chin% "YES", 0L, 1L))]
setorder(h, tumour, CHROM, POS, REF, ALT, Gene, rank, canon)
h <- unique(h, by = c("tumour", "CHROM", "POS", "REF", "ALT", "Gene"))      # most severe, canonical first
h[, c("rank", "canon", "src") := NULL]

# PASS membership per set (wide), only variants with a harmful driver-ortholog entry
hits <- merge(pass, unique(h[, .(tumour, CHROM, POS, REF, ALT)]), by = c("tumour", "CHROM", "POS", "REF", "ALT"))
memb <- dcast(hits, tumour + CLASS + CHROM + POS + REF + ALT ~ SET, fun.aggregate = length, value.var = "SET")
for (s in SET_LEVELS) { if (!s %in% names(memb)) memb[, (s) := 0L]; set(memb, j = s, value = memb[[s]] > 0) }
memb[, n_sets := rowSums(.SD), .SDcols = SET_LEVELS]
setcolorder(memb, c("tumour", "CLASS", "CHROM", "POS", "REF", "ALT", "n_sets", SET_LEVELS))

driver_hits <- h[memb, on = .(tumour, CHROM, POS, REF, ALT)]
driver_hits <- d_orth[, .(gene_id_m, symbol_h, homology_type, perc_id_h_to_m)][driver_hits, on = c(gene_id_m = "Gene"), allow.cartesian = TRUE]
driver_hits[, patient_driver := symbol_h %chin% patient_genes]
setnames(driver_hits, c("SYMBOL"), c("symbol_m_vep"), skip_absent = TRUE)
setcolorder(driver_hits, intersect(c("tumour", "symbol_h", "symbol_m_vep", "gene_id_m", "patient_driver", "homology_type", "perc_id_h_to_m",
  "CLASS", "CHROM", "POS", "REF", "ALT", "FROM_MNV", "IMPACT", "Consequence", "HGVSp", "HGVSc", "SIFT", "CANONICAL",
  "Feature", "n_sets", SET_LEVELS), names(driver_hits)))
setorder(driver_hits, tumour, symbol_h, CHROM, POS)

# ---- 6. Per gene list: tables and plot ----------------------------------------
res_n <- list()
for (L in names(LISTS)) {
  gl <- LISTS[[L]]
  dh <- gl[driver_hits, on = "symbol_h", nomatch = NULL]                   # hits in this list's genes
  setcolorder(dh, intersect(c("tumour", "symbol_h", "symbol_m_vep", "ROLE", "TIER", "in_CGC", "n_cohorts", "patient_driver"), names(dh)))
  setorder(dh, tumour, symbol_h, CHROM, POS)
  message("== ", L, ": ", uniqueN(dh[, .(tumour, CHROM, POS, REF, ALT)]), " harmful PASS variants in ", uniqueN(dh$symbol_h), " genes")
  print(dh[, .(tumour, symbol_h, symbol_m_vep, patient_driver, CLASS, IMPACT, Consequence, HGVSp, n_sets)])
  no_orth <- gl[!symbol_h %chin% d_orth$symbol_h]
  res_n[[L]] <- c(genes = nrow(gl), without_ortholog = nrow(no_orth), variants = uniqueN(dh[, .(tumour, CHROM, POS, REF, ALT)]),
                  genes_hit = uniqueN(dh$symbol_h))

  # Gene x set per tumour: number of harmful PASS variants, any HIGH
  gs <- melt(dh, id.vars = c("tumour", "symbol_h", "patient_driver", "CHROM", "POS", "REF", "ALT", "IMPACT"),
             measure.vars = SET_LEVELS, variable.name = "SET", value.name = "PASS")[PASS == TRUE]
  gs <- gs[, .(n_variants = uniqueN(paste(CHROM, POS, REF, ALT)), n_HIGH = uniqueN(paste(CHROM, POS, REF, ALT)[IMPACT == "HIGH"])),
           by = .(tumour, symbol_h, patient_driver, SET)]
  print(dcast(gs, tumour + symbol_h ~ SET, value.var = "n_variants", fill = 0))

  if (nrow(gs)) {
    gs[, `:=`(SET = factor(SET, levels = SET_LEVELS), worst = fifelse(n_HIGH > 0, "HIGH", "MODERATE only"),
              gene = fifelse(patient_driver, paste0(symbol_h, " *"), symbol_h))]
    p <- ggplot(gs, aes(SET, gene, fill = worst)) +
      geom_tile(colour = "white") + geom_text(aes(label = n_variants), size = 3) +
      facet_wrap(~ tumour, nrow = 1) +
      scale_fill_manual(values = c(HIGH = "#D7301F", "MODERATE only" = "#FDBB84"), name = "VEP impact") +
      scale_y_discrete(limits = rev) +
      labs(x = NULL, y = paste0("human gene (", L, "), mouse ortholog mutated"),
           title = paste0("Harmful (VEP HIGH/MODERATE) PASS mutations in mouse orthologs of ", L, " genes"),
           subtitle = "number = PASS variants per gene and set; * = patient (RJALS) tumourevo driver gene") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, paste0("orthologs_harmful_by_set_", L), 16, 2 + 0.25 * uniqueN(gs$gene))
  }
  fwrite(gl, file.path(OUT, paste0("genes_", L, ".tsv")), sep = "\t")
  fwrite(no_orth, file.path(OUT, paste0("genes_without_mouse_ortholog_", L, ".tsv")), sep = "\t")
  fwrite(dh, file.path(OUT, paste0("orthologs_harmful_variants_", L, ".tsv")), sep = "\t")
  fwrite(gs, file.path(OUT, paste0("orthologs_harmful_by_gene_set_", L, ".tsv")), sep = "\t")
}

# ---- 7. Outputs --------------------------------------------------------------
fwrite(d_orth, file.path(OUT, "gene_orthologs.tsv"), sep = "\t")
fwrite(data.table(patient_driver_gene = patient_genes), file.path(OUT, "patient_driver_genes.tsv"), sep = "\t")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
writeLines(c(paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_driver_genes_human"),
             paste0("- Script: ", SCRIPT, "  git: ", git_hash, "  (", R.version.string, ")"),
             paste0("- Inputs: ", SETTINGS, "; ", INTOGEN, " (", CANCER, "); CGC ", if (HAVE_CGC) CGC else "not used", "; ", ORTHOLOGS, "; ", PATIENT_RDS),
             paste0("- Parameters: impacts ", paste(IMPACTS, collapse = "/"), "; sets ", paste(SET_LEVELS, collapse = ", ")),
             paste0("- ", names(res_n), ": ", sapply(res_n, function(v) paste(names(v), v, sep = " ", collapse = ", "))),
             paste0("- Patient driver genes: ", paste(patient_genes, collapse = ", ")),
             paste0("- Outputs: ", OUT), "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
