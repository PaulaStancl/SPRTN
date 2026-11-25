#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_driver_lollipop.R - lollipop plots (trackViewer) of the patient's driver genes: the mouse
# orthologs with the mouse mutations per setting, and the human protein with the patient's mutations.
#
# Logic (agreed 2026-11-25):
#   genes      human genes with is_driver == TRUE in the RJALS tumourevo driver.rds (tumour element)
#              -> mouse orthologs (Ensembl Compara table, h2m resources)
#   mouse      every PASS variant (any of the 9 sets of mouse_wgs_snv_indel_settings.R, settings.rds) with a
#              protein position on the canonical transcript, any VEP impact (HIGH / MODERATE / LOW);
#              colour = impact, label = tumour + HGVSp; one track per set that has >= 1 variant (sets without
#              any are named in the title); MNVs split as in the sets (an MNV's CSQ goes to each SNV)
#   patient    coding mutations of the human gene in the tumourevo table (RJALS_RJALS_Tm), top track
#   axis       amino acids of the canonical transcript; protein length and Pfam domains from the Ensembl
#              REST API (rest.ensembl.org); without internet: no domains, axis to the last mutation
#   table      every coding variant in these genes x the 9 sets (TRUE/FALSE), mouse and patient
#
# Needs: settings.rds (mouse_wgs_snv_indel_settings.R). Runs in its OWN env r-lollipop (trackViewer cannot be
#   solved into r-mutation, 2026-11-25); bcftools is taken from r-mutation's bin/ (ENV_BIN). Create once:
#   micromamba create -p /common/WORK/pstancl/envs/r-lollipop -c conda-forge -c bioconda \
#     r-base=4.5 bioconductor-trackviewer bioconductor-genomicranges r-data.table r-jsonlite
# Blocks: 0 EDIT | 1 setup | 2 genes | 3 mouse coding variants per set | 4 patient variants
#         | 5 protein length + domains | 6 lollipop plots | 7 outputs, RUN_LOG.md
# Env on ssi: r-lollipop (not r-mutation). Output: results/mouse_wgs/driver_genes_human/lollipop/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO      <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
SETTINGS  <- file.path(REPO, "results/mouse_wgs/mouse_wgs_snv_indel_settings/settings.rds")
OUT       <- file.path(REPO, "results/mouse_wgs/driver_genes_human/lollipop")
SCRIPT    <- file.path(REPO, "scripts/mouse_wgs/driver_genes_human/mouse_wgs_driver_lollipop.R")
ORTHOLOGS <- "/common/WORK/pstancl/references/cross_species/ensembl110/orthologs_hsapiens_mmusculus.tsv"
PATIENT_SAMPLE <- "RJALS_RJALS_Tm"
PATIENT_RDS <- file.path(REPO, "results/wgs/tumourevo/RJALS/driver_annotation/annotate_driver/RJALS/RJALS/RJALS_RJALS_Tm/RJALS_RJALS_RJALS_RJALS_Tm_driver.rds")
REST      <- "https://rest.ensembl.org"
ENV_BIN   <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS  <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
IMPACT_COLS <- c(HIGH = "#D7301F", MODERATE = "#FC8D59", LOW = "#9ECAE1", MODIFIER = "grey70")

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
if (!requireNamespace("trackViewer", quietly = TRUE)) stop("trackViewer not installed - run this script in the r-lollipop env (header)")
if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS)
for (f in c(SETTINGS, ORTHOLOGS, PATIENT_RDS)) if (!file.exists(f)) stop("Not found: ", f)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
first_int <- function(x) suppressWarnings(as.integer(sub("^([0-9]+).*$", "\\1", x)))     # "123-124" -> 123, "?" -> NA

# ---- 2. Genes: patient drivers and their mouse orthologs ----------------------
px <- readRDS(PATIENT_RDS)
if (!PATIENT_SAMPLE %in% names(px)) stop("PATIENT_RDS: no element ", PATIENT_SAMPLE)
pm <- as.data.table(px[[PATIENT_SAMPLE]]$mutations)
genes_h <- unique(as.character(pm[is_driver %in% TRUE, SYMBOL]))
orth <- fread(ORTHOLOGS, colClasses = "character", na.strings = c("", "NA"))
go <- orth[symbol_h %chin% genes_h & !is.na(gene_id_m), .(symbol_h, gene_id_m, symbol_m, homology_type)]
print(go)
if (!nrow(go)) stop("No mouse orthologs for ", paste(genes_h, collapse = ", "))

# ---- 3. Mouse coding variants per set ---------------------------------------
st <- readRDS(SETTINGS)
pass <- as.data.table(st$pass)[, .(tumour, SET = as.character(SET), CHROM = as.character(CHROM), POS = as.integer(POS),
                                    REF = as.character(REF), ALT = as.character(ALT))]
SET_LEVELS <- st$params$SET_LEVELS
files <- as.data.table(st$files)
ids_file <- tempfile(fileext = ".txt"); writeLines(unique(go$gene_id_m), ids_file)
read_csq <- function(i) {
  f <- files[i]
  hdr <- system2(BCFTOOLS, c("view", "-h", shQuote(f$path)), stdout = TRUE)
  h <- grep("^##INFO=<ID=CSQ,", hdr, value = TRUE)
  if (!length(h)) return(NULL)
  fields <- strsplit(sub("\".*$", "", sub("^.*Format: ", "", h)), "|", fixed = TRUE)[[1]]
  cmd <- paste(shQuote(BCFTOOLS), "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CSQ\\n"), shQuote(f$path),
               "| grep -F -f", shQuote(ids_file), "|| true")
  x <- tryCatch(fread(cmd = cmd, sep = "\t", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", "CSQ"),
                      colClasses = list(character = c("CHROM", "REF", "ALT", "CSQ"), integer = "POS"), quote = ""),
                error = function(e) data.table())
  if (!nrow(x)) return(NULL)
  x <- x[, .(entry = unlist(strsplit(CSQ, ",", fixed = TRUE))), by = .(CHROM, POS, REF, ALT)]
  m <- tstrsplit(x$entry, "|", fixed = TRUE, fill = NA_character_); length(m) <- length(fields); names(m) <- fields
  keep <- c("Consequence", "IMPACT", "Gene", "Feature", "HGVSp", "Protein_position", "Amino_acids", "CANONICAL")
  x <- cbind(x[, .(CHROM, POS, REF, ALT)], as.data.table(lapply(setNames(keep, keep), function(k)
    if (k %in% fields && !is.null(m[[k]])) m[[k]] else rep(NA_character_, nrow(x)))))
  x[Gene %chin% go$gene_id_m]
}
csq <- unique(rbindlist(lapply(seq_len(nrow(files)), read_csq), fill = TRUE))
if (!nrow(csq)) stop("No CSQ entries on ", paste(go$symbol_m, collapse = ", "), " in any file")
# MNVs -> SNVs (as the sets)
csq[, MNV := nchar(REF) == nchar(ALT) & nchar(REF) > 1L & REF %like% "^[ACGT]+$" & ALT %like% "^[ACGT]+$"]
mn <- csq[(MNV)]
if (nrow(mn)) { mn <- mn[rep(seq_len(.N), nchar(REF))][, k := sequence(nchar(mn$REF))]
  mn[, `:=`(POS = POS + k - 1L, REF = substring(REF, k, k), ALT = substring(ALT, k, k))]; mn <- mn[REF != ALT][, k := NULL] }
csq <- rbind(csq[!(MNV)], mn)[, MNV := NULL]
csq[, `:=`(CHROM = as.character(CHROM), POS = as.integer(POS))]          # same types as pass (newer data.table reads 1..19 as int)

# Canonical transcript per mouse gene (most frequent CANONICAL == YES feature; else most frequent feature)
pick_tx <- function(d) {
  can <- d[CANONICAL %chin% "YES" & !is.na(Feature), .N, by = Feature][order(-N)]
  if (nrow(can)) can$Feature[1] else d[!is.na(Feature), .N, by = Feature][order(-N)]$Feature[1]
}
tx_m <- csq[, .(Feature = pick_tx(.SD)), by = Gene]
mc <- csq[tx_m, on = .(Gene, Feature)][!is.na(first_int(Protein_position))]
mc[, aa := first_int(Protein_position)]
mc <- go[, .(gene_id_m, symbol_h, symbol_m)][mc, on = c(gene_id_m = "Gene")]

# PASS membership per set
mv <- merge(pass, unique(mc[, .(CHROM, POS, REF, ALT)]), by = c("CHROM", "POS", "REF", "ALT"))
mouse_tab <- dcast(mv, tumour + CHROM + POS + REF + ALT ~ SET, fun.aggregate = length, value.var = "SET")
for (s in SET_LEVELS) { if (!s %in% names(mouse_tab)) mouse_tab[, (s) := 0L]; set(mouse_tab, j = s, value = mouse_tab[[s]] > 0) }
mouse_tab <- unique(mc[, .(symbol_h, symbol_m, Feature, CHROM, POS, REF, ALT, aa, IMPACT, Consequence, HGVSp)])[mouse_tab, on = .(CHROM, POS, REF, ALT)]
mouse_tab[, n_sets := rowSums(.SD), .SDcols = SET_LEVELS]
setcolorder(mouse_tab, c("symbol_h", "symbol_m", "tumour", "CHROM", "POS", "REF", "ALT", "aa", "IMPACT", "Consequence", "HGVSp", "n_sets", SET_LEVELS))
setorder(mouse_tab, symbol_h, aa, tumour)
print(mouse_tab[, .(symbol_m, tumour, aa, IMPACT, Consequence, HGVSp, n_sets)])

# ---- 4. Patient coding variants -------------------------------------------------
need <- c("SYMBOL", "IMPACT", "Protein_position", "Feature")
if (!all(need %in% names(pm))) stop("tumourevo table lacks: ", paste(setdiff(need, names(pm)), collapse = ", "))
pv <- pm[SYMBOL %chin% genes_h & !is.na(first_int(as.character(Protein_position)))]
pv <- pv[, .(symbol_h = as.character(SYMBOL), Feature = as.character(Feature), aa = first_int(as.character(Protein_position)),
             IMPACT = as.character(IMPACT), HGVSp = if ("HGVSp" %in% names(pv)) as.character(HGVSp) else NA_character_,
             Consequence = if ("Consequence" %in% names(pv)) as.character(Consequence) else NA_character_,
             is_driver = is_driver %in% TRUE, CANONICAL = if ("CANONICAL" %in% names(pv)) as.character(CANONICAL) else NA_character_)]
tx_h <- pv[, .(Feature = pick_tx(.SD)), by = symbol_h]
pv <- pv[tx_h, on = .(symbol_h, Feature)]
print(pv)

# ---- 5. Protein length and Pfam domains (Ensembl REST) ---------------------------
rest_json <- function(path) tryCatch(jsonlite::fromJSON(paste0(REST, path, if (grepl("\\?", path)) ";" else "?", "content-type=application/json")),
                                     error = function(e) { message("REST failed (", path, "): ", conditionMessage(e)); NULL })
protein_info <- function(tx) {
  tx <- sub("\\..*$", "", tx)
  lk <- rest_json(paste0("/lookup/id/", tx, "?expand=1"))
  if (is.null(lk) || is.null(lk$Translation)) return(list(len = NA_integer_, dom = data.table()))
  dom <- rest_json(paste0("/overlap/translation/", lk$Translation$id, "?feature=protein_feature;type=Pfam"))
  dom <- if (is.data.frame(dom) && nrow(dom)) as.data.table(dom)[, .(start, end, name = fifelse(nzchar(description), description, id))] else data.table()
  list(len = as.integer(lk$Translation$length), dom = unique(dom))
}

# ---- 6. Lollipop plots --------------------------------------------------------
make_gr <- function(d, seqname) {
  gr <- GRanges(seqname, IRanges(d$aa, width = 1))
  gr$color <- unname(IMPACT_COLS[d$IMPACT]); gr$color[is.na(gr$color)] <- "grey70"
  gr$border <- "grey30"
  names(gr) <- d$label
  gr
}
make_feat <- function(dom, seqname, len) {
  if (!nrow(dom)) return(GRanges(seqname, IRanges(1, len), fill = "grey90", height = 0.02))
  ft <- GRanges(seqname, IRanges(dom$start, dom$end))
  ft$fill <- grDevices::hcl.colors(max(3, nrow(dom)), "Set 2")[seq_len(nrow(dom))]
  ft$height <- 0.04
  names(ft) <- dom$name
  ft
}
legend <- list(labels = names(IMPACT_COLS), fill = unname(IMPACT_COLS))
plotted <- character()
for (g in unique(go$symbol_h)) {
  gm <- go[symbol_h == g]
  tab_g <- mouse_tab[symbol_h == g]
  sym_m <- gm$symbol_m[1]
  pi_m <- protein_info(tx_m[Gene %chin% gm$gene_id_m, Feature][1])
  pi_h <- if (nrow(tx_h[symbol_h == g])) protein_info(tx_h[symbol_h == g, Feature][1]) else list(len = NA_integer_, dom = data.table())
  len_m <- if (is.na(pi_m$len)) max(c(tab_g$aa, 100L), na.rm = TRUE) + 20L else pi_m$len
  len_h <- if (is.na(pi_h$len)) max(c(pv[symbol_h == g, aa], 100L), na.rm = TRUE) + 20L else pi_h$len

  snps <- list(); feats <- list(); rngs <- list()
  p_g <- pv[symbol_h == g]
  if (nrow(p_g)) {
    p_g[, label := paste0(fifelse(is_driver, "* ", ""), fifelse(is.na(HGVSp) | !nzchar(HGVSp), paste0("aa", aa), sub("^.*:", "", HGVSp)))]
    nm <- paste0("patient ", g)
    snps[[nm]] <- make_gr(p_g, g); feats[[nm]] <- make_feat(pi_h$dom, g, len_h); rngs[[nm]] <- GRanges(g, IRanges(1, len_h))
  }
  for (s in SET_LEVELS) {
    d <- tab_g[get(s) == TRUE]
    if (!nrow(d)) next
    d[, label := paste(tumour, fifelse(is.na(HGVSp) | !nzchar(HGVSp), paste0("aa", aa), sub("^.*:", "", HGVSp)))]
    snps[[s]] <- make_gr(d, sym_m); feats[[s]] <- make_feat(pi_m$dom, sym_m, len_m); rngs[[s]] <- GRanges(sym_m, IRanges(1, len_m))
  }
  empty <- setdiff(SET_LEVELS, names(snps))
  if (length(snps) == 0) { message(g, ": no coding variants to plot"); next }
  ttl <- paste0(g, " (patient, human) / ", sym_m, " (mouse) - PASS coding variants per setting",
                if (length(empty)) paste0("\nno variant in: ", paste(empty, collapse = ", ")))
  h <- 2 + 1.6 * length(snps)
  for (dev in c("pdf", "png")) {
    fn <- file.path(OUT, paste0("lollipop_", g, ".", dev))
    if (dev == "pdf") pdf(fn, width = 14, height = h) else png(fn, width = 14, height = h, units = "in", res = 150)
    ok <- tryCatch({
      trackViewer::lolliplot(snps, feats, ranges = rngs, legend = rep(list(legend), length(snps)), ylab = names(snps), cex = 0.7)
      grid::grid.text(ttl, x = 0.5, y = 0.985, just = "top", gp = grid::gpar(cex = 0.9, fontface = "bold"))
      TRUE }, error = function(e) { message("lolliplot ", g, " failed: ", conditionMessage(e)); FALSE })
    dev.off()
    if (ok && dev == "pdf") plotted <- c(plotted, g)
  }
}

# ---- 7. Outputs --------------------------------------------------------------
fwrite(mouse_tab, file.path(OUT, "mouse_coding_variants_by_set.tsv"), sep = "\t")
fwrite(pv, file.path(OUT, "patient_coding_variants.tsv"), sep = "\t")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
writeLines(c(paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_driver_lollipop"),
             paste0("- Script: ", SCRIPT, "  git: ", git_hash, "  (", R.version.string, ", trackViewer ", packageVersion("trackViewer"), ")"),
             paste0("- Inputs: ", SETTINGS, "; ", ORTHOLOGS, "; ", PATIENT_RDS, " (", PATIENT_SAMPLE, "); domains from ", REST),
             paste0("- Genes: ", paste(genes_h, collapse = ", "), " -> ", paste(go$symbol_m, collapse = ", ")),
             paste0("- Mouse coding variants: ", nrow(mouse_tab), "; patient coding variants: ", nrow(pv), "; plotted: ", paste(plotted, collapse = ", ")),
             paste0("- Outputs: ", OUT), "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
