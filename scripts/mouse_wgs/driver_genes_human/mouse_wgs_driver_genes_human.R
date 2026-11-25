#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_driver_genes_human.R - harmful mutations in mouse orthologs of human HCC driver genes.
#
# Logic (agreed 2026-11-25):
#   driver genes  IntOGen Compendium_Cancer_Genes.tsv (tumourevo's), rows with CANCER_TYPE == "HCC"
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
# Blocks: 0 EDIT | 1 setup | 2 driver genes + orthologs | 3 PASS variants per set | 4 CSQ of driver orthologs
#         | 5 harmful hits | 6 plots | 7 outputs, RUN_LOG.md
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
ORTHOLOGS <- "/common/WORK/pstancl/references/cross_species/ensembl110/orthologs_hsapiens_mmusculus.tsv"
IMPACTS   <- c("HIGH", "MODERATE")
ENV_BIN   <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS  <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
for (f in c(SETTINGS, INTOGEN, ORTHOLOGS)) if (!file.exists(f)) stop("Not found: ", f,
  "\n  settings.rds: run mouse_wgs_snv_indel_settings.R; orthologs: run scripts/cross_species/h2m 00_resources.R;",
  "\n  IntOGen: set INTOGEN (find command in block 0)")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
save_plot <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150, limitsize = FALSE)
}

# ---- 2. Driver genes and their mouse orthologs ----------------------------------
ig <- fread(INTOGEN)
stopifnot(all(c("SYMBOL", "CANCER_TYPE") %in% names(ig)))
drivers <- ig[CANCER_TYPE == CANCER, .(ROLE = if ("ROLE" %in% names(ig)) paste(sort(unique(na.omit(ROLE))), collapse = ";") else NA_character_),
              by = .(symbol_h = SYMBOL)]
message(nrow(drivers), " IntOGen ", CANCER, " driver genes")

orth <- fread(ORTHOLOGS, colClasses = "character", na.strings = c("", "NA"))
d_orth <- orth[!is.na(gene_id_m), .(symbol_h, gene_id_h, gene_id_m, symbol_m, homology_type, perc_id_h_to_m)][drivers, on = "symbol_h"]
no_ortholog <- d_orth[is.na(gene_id_m), .(symbol_h, ROLE)]
d_orth <- d_orth[!is.na(gene_id_m)]
print(d_orth[, .N, by = homology_type])
message(nrow(no_ortholog), " driver genes without a mouse ortholog in the table (or symbol not found): ",
        paste(head(no_ortholog$symbol_h, 30), collapse = ", "))

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
driver_hits <- d_orth[, .(gene_id_m, symbol_h, ROLE, homology_type, perc_id_h_to_m)][driver_hits, on = c(gene_id_m = "Gene"), allow.cartesian = TRUE]
setnames(driver_hits, c("SYMBOL"), c("symbol_m_vep"), skip_absent = TRUE)
setcolorder(driver_hits, intersect(c("tumour", "symbol_h", "symbol_m_vep", "gene_id_m", "ROLE", "homology_type", "perc_id_h_to_m",
  "CLASS", "CHROM", "POS", "REF", "ALT", "FROM_MNV", "IMPACT", "Consequence", "HGVSp", "HGVSc", "SIFT", "CANONICAL",
  "Feature", "n_sets", SET_LEVELS), names(driver_hits)))
setorder(driver_hits, tumour, symbol_h, CHROM, POS)
print(driver_hits[, .(tumour, symbol_h, symbol_m_vep, CLASS, IMPACT, Consequence, HGVSp, n_sets)])

# Gene x set per tumour: number of harmful PASS variants
gene_set <- melt(driver_hits, id.vars = c("tumour", "symbol_h", "CHROM", "POS", "REF", "ALT", "IMPACT"),
                 measure.vars = SET_LEVELS, variable.name = "SET", value.name = "PASS")[PASS == TRUE]
gene_set <- gene_set[, .(n_variants = uniqueN(paste(CHROM, POS, REF, ALT)), n_HIGH = uniqueN(paste(CHROM, POS, REF, ALT)[IMPACT == "HIGH"])),
                     by = .(tumour, symbol_h, SET)]
print(dcast(gene_set, tumour + symbol_h ~ SET, value.var = "n_variants", fill = 0))

# ---- 6. Plot: gene x set per tumour, number of harmful PASS variants (label), any HIGH (fill) ----
if (nrow(gene_set)) {
  gene_set[, `:=`(SET = factor(SET, levels = SET_LEVELS), worst = fifelse(n_HIGH > 0, "HIGH", "MODERATE only"))]
  p <- ggplot(gene_set, aes(SET, symbol_h, fill = worst)) +
    geom_tile(colour = "white") + geom_text(aes(label = n_variants), size = 3) +
    facet_wrap(~ tumour, nrow = 1) +
    scale_fill_manual(values = c(HIGH = "#D7301F", "MODERATE only" = "#FDBB84"), name = "VEP impact") +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = "human HCC driver gene (IntOGen), mouse ortholog mutated",
         title = "Harmful (VEP HIGH/MODERATE) PASS mutations in mouse orthologs of human HCC driver genes",
         subtitle = "number = PASS variants per gene and set") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "driver_orthologs_harmful_by_set", 16, 2 + 0.25 * uniqueN(gene_set$symbol_h))
}

# ---- 7. Outputs --------------------------------------------------------------
fwrite(drivers, file.path(OUT, "intogen_drivers.tsv"), sep = "\t")
fwrite(d_orth, file.path(OUT, "driver_orthologs.tsv"), sep = "\t")
fwrite(no_ortholog, file.path(OUT, "drivers_without_mouse_ortholog.tsv"), sep = "\t")
fwrite(driver_hits, file.path(OUT, "driver_orthologs_harmful_variants.tsv"), sep = "\t")
fwrite(gene_set, file.path(OUT, "driver_orthologs_harmful_by_gene_set.tsv"), sep = "\t")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
writeLines(c(paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_driver_genes_human"),
             paste0("- Script: ", SCRIPT, "  git: ", git_hash, "  (", R.version.string, ")"),
             paste0("- Inputs: ", SETTINGS, "; ", INTOGEN, " (", CANCER, ", ", nrow(drivers), " genes); ", ORTHOLOGS),
             paste0("- Parameters: impacts ", paste(IMPACTS, collapse = "/"), "; sets ", paste(SET_LEVELS, collapse = ", ")),
             paste0("- Result: ", uniqueN(driver_hits[, .(tumour, CHROM, POS, REF, ALT)]), " harmful PASS variants in ",
                    uniqueN(driver_hits$symbol_h), " driver-gene orthologs"),
             paste0("- Outputs: ", OUT), "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
