#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_snv_indel_settings.R - compare the three calling settings per tumour.
#
# Sets (per tumour; PASS = FILTER == PASS in each file's own FILTER, never the merged file):
#   A_freebayes            A: tumor_vs_normal, FreeBayes PASS (FILTER as it is in the file)
#   A_freebayes_pass_expr  A: FreeBayes PASS AND the prefilter expression (below) on the split record
#   A_freebayes_expr       A: FreeBayes, prefilter expression only; the FILTER column is not used
#   A_mutect2_strelka    A: PASS first in Mutect2's and in Strelka's own file, then intersected
#                        ("all records" = found by both, any FILTER)
#   A_strelka            A: Strelka PASS
#   A_mutect2            A: Mutect2 PASS
#   B_mutect2_joint_wt   B: tumor_vs_joint_wt/  Mutect2 joint calls with the WT normals, PASS
#   C_mutect2_pon        C: tumor_only_pon/     Mutect2 tumour-only + WT panel of normals, PASS
#
# FreeBayes prefilter expression (provider), evaluated in R on each record of the split file:
#   INFO/MQM >= 40 && GT[normal]="0/0" && GT[tumor]="0/1" && FMT/DP[normal] >= 30 && FMT/DP[tumor] >= 30
#   && QUAL/INFO/AO >= 5 && SAF > 0 && SAR > 0      (missing value = fails)
# The file's FILTER does not follow it record by record (mouse_wgs_freebayes_pass_check.R), hence the 3 sets.
#
# Questions answered, per tumour:
#   1. SNV / indel counts per set (all records vs PASS)
#   2. proportion of PASS variants shared between the sets (pairwise, exclusive categories)
#   3. SBS96 profile of each set, cosine similarity between sets and between tumours
#
# Blocks:
#   0  EDIT THIS BLOCK
#   1  setup and helpers
#   2  files (settings A, B, C) and the samples in each
#   3  read records; multiallelics and MNVs split; the sets
#   4  counts per set (+ qcVCF plot_mutation_counts)
#   5  sharing between sets (+ qcVCF plot_pairwise_shared_mutations), exclusive categories
#   6  SBS96: profiles, plots (ggplot + qcVCF plot96_matrix), cosine similarity
#   7  outputs, summary.md, RUN_LOG.md
#
# Comparison unit: chrom:pos:ref:alt. MNVs are split into SNVs (decided 2026-11-25); doublet halves
# (adjacent SNVs within one set) are left out of SBS96 (one DBS event). Indels are not left-aligned.
#
# Env on ssi: r-mutation (data.table, ggplot2, qcVCF, bcftools, samtools).
# Output: results/mouse_wgs/mouse_wgs_snv_indel_settings/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
DATA     <- file.path(REPO, "data/mouse/wgs")
DIR_A    <- file.path(DATA, "tumor_vs_normal/snv_indel")
DIR_B    <- file.path(DATA, "tumor_vs_joint_wt")
DIR_C    <- file.path(DATA, "tumor_only_pon")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_snv_indel_settings")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_snv_indel_settings.R")

# GRCm38 FASTA (Ensembl 102 primary assembly, from scripts/cross_species/h2m 00_resources.R), + .fai
FASTA    <- "/common/WORK/pstancl/references/cross_species/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa"

# Tools: RStudio often does not see the env's PATH, so fall back to the env's bin/
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
SAMTOOLS <- Sys.which("samtools"); if (!nzchar(SAMTOOLS)) SAMTOOLS <- file.path(ENV_BIN, "samtools")

TUMOURS <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1",
             HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
             HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")
SET_LABELS <- c(A_freebayes        = "A: tumor_vs_normal, FreeBayes PASS",
                A_freebayes_pass_expr = "A: tumor_vs_normal, FreeBayes PASS & prefilter expression",
                A_freebayes_expr   = "A: tumor_vs_normal, FreeBayes prefilter expression only (FILTER not used)",
                A_mutect2_strelka  = "A: tumor_vs_normal, Mutect2 PASS & Strelka PASS (intersect)",
                A_strelka          = "A: tumor_vs_normal, Strelka PASS",
                A_mutect2          = "A: tumor_vs_normal, Mutect2 PASS",
                B_mutect2_joint_wt = "B: tumor_vs_joint_wt, Mutect2 PASS",
                C_mutect2_pon      = "C: tumor_only_pon, Mutect2 PASS")
SET_LEVELS <- names(SET_LABELS)


# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# Console as wide as the terminal (falls back to 80 when there is no terminal, e.g. batch Rscript)
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
HAVE_QCVCF <- requireNamespace("qcVCF", quietly = TRUE)
if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS, " - install it in r-mutation or set BCFTOOLS")
HAVE_FASTA <- file.exists(SAMTOOLS) && file.exists(FASTA) && file.exists(paste0(FASTA, ".fai"))
if (!HAVE_FASTA) message("No samtools or FASTA(.fai): the SBS96 part (block 6) will be skipped")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()

# bcftools: stdout lines; on failure shows bcftools' own message
bcf <- function(args) {
  err <- tempfile()
  out <- tryCatch(suppressWarnings(system2(BCFTOOLS, args, stdout = TRUE, stderr = err)),
                  error = function(e) structure(character(), status = -1L, msg = conditionMessage(e)))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0)
    stop("bcftools failed: bcftools ", paste(args, collapse = " "), "\n",
         paste(c(attr(out, "msg"), readLines(err, warn = FALSE)), collapse = "\n"))
  out
}
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150, limitsize = FALSE)
  invisible(p)
}
# One row per ALT of multiallelic records
split_alt <- function(x) {
  m <- x[ALT %like% ","]
  if (!nrow(m)) return(x)
  s <- strsplit(m$ALT, ",", fixed = TRUE)
  rbind(x[!ALT %like% ","], m[rep(seq_len(.N), lengths(s))][, ALT := unlist(s)])
}
# Split MNVs into their SNVs (positions where REF and ALT differ)
split_mnv <- function(x) {
  m <- x[TYPE == "MNV"]
  if (!nrow(m)) return(x[, FROM_MNV := FALSE][])
  n <- nchar(m$REF)
  e <- m[rep(seq_len(.N), n)]
  e[, k := sequence(n)]
  e[, `:=`(POS = POS + k - 1L, REF = substring(REF, k, k), ALT = substring(ALT, k, k))]
  e <- e[REF != ALT][, `:=`(TYPE = "SNV", FROM_MNV = TRUE, k = NULL)]
  rbind(x[TYPE != "MNV"][, FROM_MNV := FALSE], e)
}

# ---- 2. Files ----------------------------------------------------------------
list_vcf <- function(d, setting) {
  p <- list.files(d, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf|vcf\\.gz)$")
  if (!length(p)) stop("No VCF/BCF files in ", d)
  data.table(setting = setting, path = p)
}
files <- rbind(list_vcf(DIR_A, "A"), list_vcf(DIR_B, "B"), list_vcf(DIR_C, "C"))
files[, rel := substring(path, nchar(DATA) + 2)]
files[, tumour := NA_character_]
for (k in names(TUMOURS)) files[grepl(TUMOURS[[k]], basename(rel), fixed = TRUE), tumour := k]
hdr <- setNames(lapply(files$path, function(f) bcf(c("view", "-h", shQuote(f)))), files$rel)
files[, has_CALLER := sapply(rel, function(r) any(startsWith(hdr[[r]], "##INFO=<ID=CALLER,")))]
files[, caller := fcase(grepl("mutect2", rel, ignore.case = TRUE), "mutect2",
                        grepl("strelka", rel, ignore.case = TRUE), "strelka",
                        grepl("freebayes", rel, ignore.case = TRUE), "freebayes", default = NA_character_)]
# Setting A: the single-caller files (not the merged file)
files[, src := fcase(setting == "A" & !has_CALLER & caller == "freebayes", "A_freebayes",
                     setting == "A" & !has_CALLER & caller == "mutect2", "A_mutect2",
                     setting == "A" & !has_CALLER & caller == "strelka", "A_strelka",
                     setting == "B", "B_mutect2_joint_wt",
                     setting == "C", "C_mutect2_pon", default = NA_character_)]
files[, .(setting, tumour, src, has_CALLER, file = rel)]
files <- files[!is.na(src)]
stopifnot(!anyNA(files$tumour))
print(dcast(files[, .N, by = .(tumour, src)], tumour ~ src, value.var = "N", fill = 0))   # expect >= 1 everywhere
miss <- CJ(tumour = names(TUMOURS), src = c("A_freebayes", "A_mutect2", "A_strelka", "B_mutect2_joint_wt", "C_mutect2_pon"))[!files, on = .(tumour, src)]
if (nrow(miss)) stop("Missing files:\n", paste(capture.output(print(miss)), collapse = "\n"))

# Samples in each file
files[, samples := sapply(path, function(f) paste(bcf(c("query", "-l", shQuote(f))), collapse = ","))]
print(files[, .(tumour, src, samples)])

# ---- 3. Read records and build the sets ----------------------------------------
read_rec <- function(i) {
  f <- files[i]
  is_fb <- f$src == "A_freebayes"
  fmt <- paste0("%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER",
                if (is_fb) "\\t%QUAL\\t%INFO/MQM\\t%INFO/AO\\t%INFO/SAF\\t%INFO/SAR[\\t%GT\\t%DP]", "\\n")
  cn  <- c("CHROM", "POS", "REF", "ALT", "FILTER",
           if (is_fb) c("QUAL", "MQM", "AO", "SAF", "SAR", "GT_s1", "DP_s1", "GT_s2", "DP_s2"))
  cmd <- paste(shQuote(BCFTOOLS), "query", "-f", shQuote(fmt), shQuote(f$path))
  x <- tryCatch(fread(cmd = cmd, sep = "\t", header = FALSE, col.names = cn, quote = "",
                      colClasses = list(character = setdiff(cn, "POS"), integer = "POS")),
                error = function(e) { message(f$rel, ": no records read (", conditionMessage(e), ")"); data.table() })
  if (!nrow(x)) return(NULL)
  if (is_fb) {
    # FreeBayes prefilter expression per record (see header); tumour/normal column from the sample names
    smp <- strsplit(f$samples, ",", fixed = TRUE)[[1]]
    it  <- grep(TUMOURS[[f$tumour]], smp, fixed = TRUE)
    if (length(smp) != 2 || length(it) != 1) stop(f$rel, ": expected 2 samples, one of them the tumour; got ", f$samples)
    g <- function(v, s) x[[paste0(v, "_s", s)]]
    n <- function(v) suppressWarnings(as.numeric(v))
    x[, RULE := n(MQM) >= 40 & g("GT", 3 - it) == "0/0" & g("GT", it) == "0/1" &
                n(g("DP", 3 - it)) >= 30 & n(g("DP", it)) >= 30 & n(QUAL) / n(AO) >= 5 & n(SAF) > 0 & n(SAR) > 0]
    x[is.na(RULE), RULE := FALSE]
    x[, c("QUAL", "MQM", "AO", "SAF", "SAR", "GT_s1", "DP_s1", "GT_s2", "DP_s2") := NULL]
  }
  x[, `:=`(tumour = f$tumour, src = f$src, rel = f$rel)]
}
raw <- rbindlist(lapply(seq_len(nrow(files)), read_rec), use.names = TRUE, fill = TRUE)
raw[, PASS := FILTER == "PASS"]
# FreeBayes expression as its own temporary source (PASS = expression TRUE); the two extra FreeBayes sets
# are built from it at variant level after the MNV split (block 3, below)
raw <- rbind(raw, raw[src == "A_freebayes"][, `:=`(src = "fb_expression", PASS = RULE)])
records_per_file <- raw[, .(records = .N, PASS = sum(PASS)), by = .(tumour, src, rel)][order(tumour, src)]
print(records_per_file)
filter_values <- raw[, .N, by = .(tumour, src, FILTER)][order(tumour, src, -N)]
print(filter_values[, head(.SD, 8), by = .(tumour, src)])          # top FILTER values per file

# Multiallelics split; only plain ACGT alleles (drops "*" and symbolic ALTs); MNVs split into SNVs.
# One row per tumour x source x variant; PASS if PASS in any of the source's records (e.g. snvs + indels files).
rec <- split_alt(raw[, .(tumour, src, CHROM, POS, REF, ALT, PASS)])
dropped_alleles <- rec[!(REF %like% "^[ACGT]+$" & ALT %like% "^[ACGT]+$"), .N, by = .(tumour, src)]
print(dropped_alleles)
rec <- rec[REF %like% "^[ACGT]+$" & ALT %like% "^[ACGT]+$"]
rec[, TYPE := fifelse(nchar(REF) == nchar(ALT), fifelse(nchar(REF) == 1L, "SNV", "MNV"), "INDEL")]
rec <- split_mnv(rec)[, .(PASS = any(PASS)), by = .(tumour, src, CLASS = fifelse(TYPE == "SNV", "SNV", "INDEL"), CHROM, POS, REF, ALT)]

# A_mutect2_strelka: PASS first in each caller's own file, then intersected. "all records" = found by both
# (own files, any FILTER); PASS = PASS in both own files
key <- c("tumour", "CLASS", "CHROM", "POS", "REF", "ALT")
m2s <- merge(rec[src == "A_mutect2", c(key, "PASS"), with = FALSE], rec[src == "A_strelka", c(key, "PASS"), with = FALSE],
             by = key, suffixes = c("_m2", "_st"))[, .(PASS = PASS_m2 & PASS_st), by = key]
one <- function(s) rec[src == s, c(key, "PASS"), with = FALSE][, SET := s]
# FreeBayes, per variant (MNVs split): file PASS x expression
#   A_freebayes_pass_expr = file PASS AND expression TRUE
#   A_freebayes_expr      = expression TRUE, FILTER column not used (expression only)
fbv <- merge(rec[src == "A_freebayes", c(key, "PASS"), with = FALSE], rec[src == "fb_expression", c(key, "PASS"), with = FALSE],
             by = key, suffixes = c("_file", "_expr"))
print(fbv[, .(variants = .N, PASS_only = sum(PASS_file & !PASS_expr), expr_only = sum(!PASS_file & PASS_expr),
              PASS_and_expr = sum(PASS_file & PASS_expr), neither = sum(!PASS_file & !PASS_expr)), by = .(tumour, CLASS)])
allset <- rbind(one("A_freebayes"),
                fbv[, c(key, "PASS_file", "PASS_expr"), with = FALSE][, .(PASS = PASS_file & PASS_expr), by = key][, SET := "A_freebayes_pass_expr"],
                fbv[, c(key, "PASS_file", "PASS_expr"), with = FALSE][, .(PASS = PASS_expr), by = key][, SET := "A_freebayes_expr"],
                m2s[, SET := "A_mutect2_strelka"], one("A_strelka"), one("A_mutect2"),
                one("B_mutect2_joint_wt"), one("C_mutect2_pon"), use.names = TRUE)
rm(m2s, fbv)
allset[, ID := paste(CHROM, POS, REF, ALT, sep = ":")]
allset[, SET := factor(SET, levels = SET_LEVELS)]
pass <- allset[PASS == TRUE]
pass[, .N, by = .(tumour, SET, CLASS)][order(tumour, SET, CLASS)]

# ---- 4. Counts per set ------------------------------------------------------------
counts_filter <- rbind(allset[, .N, by = .(tumour, SET, CLASS)][, filter := "all"],
                       pass[, .N, by = .(tumour, SET, CLASS)][, filter := "PASS"])
counts_filter[, filter := factor(filter, levels = c("all", "PASS"))]
setorder(counts_filter, tumour, filter, SET, CLASS)
counts <- dcast(counts_filter, tumour + SET + CLASS ~ filter, value.var = "N", fill = 0)
counts[, pct_PASS := round(100 * PASS / all, 1)]
print(counts)

p_counts <- ggplot(counts_filter, aes(SET, N, fill = CLASS)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = N), position = position_dodge(width = 0.9), vjust = -0.3, size = 2.8) +
  facet_grid(filter ~ tumour, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "variants (MNVs split)", fill = NULL,
       title = "SNVs and indels per setting, all records vs PASS",
       subtitle = "A_mutect2_strelka: all = found by both; PASS = PASS in both own files") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.subtitle = element_text(size = 8))
save_plot(p_counts, "counts_per_setting", 16, 8)

# PASS only, one panel per tumour (the main count plot)
p_counts_pass <- ggplot(counts_filter[filter == "PASS"], aes(SET, N, fill = CLASS)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = N), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  facet_wrap(~ tumour, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "PASS variants (MNVs split)", fill = NULL, title = "PASS SNVs and indels per setting") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.subtitle = element_text(size = 8))
save_plot(p_counts_pass, "counts_per_setting_PASS", 16, 5.5)

# qcVCF: PASS counts per set (classes SNV / INDEL after the MNV split)
if (HAVE_QCVCF) {
  for (t in names(TUMOURS)) {
    q <- pass[tumour == t, .(CHROM, POS, REF, ALT, subjectID = tumour, tool = as.character(SET), mutType = CLASS)]
    if (!nrow(q)) next
    pq <- tryCatch(qcVCF::plot_mutation_counts(q, grp_col = "tool", cohort = t),
                   error = function(e) { message("qcVCF plot_mutation_counts ", t, " skipped: ", conditionMessage(e)); NULL })
    if (!is.null(pq)) save_plot(pq, paste0("qcVCF_mutation_counts_", t), 8, 6)
  }
}

# ---- 5. Sharing between sets (PASS) ----------------------------------------------
# 5a. Pairwise: shared, % of the row set, Jaccard
pairwise <- pass[, {
  s <- split(ID, SET, drop = TRUE)
  pr <- CJ(set_1 = names(s), set_2 = names(s))
  pr[, `:=`(n_1 = lengths(s[set_1]), n_2 = lengths(s[set_2]),
            shared = mapply(function(a, b) length(intersect(s[[a]], s[[b]])), set_1, set_2))]
  pr[, `:=`(pct_of_1 = round(100 * shared / n_1, 1), jaccard = round(shared / (n_1 + n_2 - shared), 3))]
  pr
}, by = .(tumour, CLASS)]
pairwise[, `:=`(set_1 = factor(set_1, levels = SET_LEVELS), set_2 = factor(set_2, levels = SET_LEVELS))]
setorder(pairwise, tumour, CLASS, set_1, set_2)
print(pairwise[set_1 != set_2])

p_pair <- ggplot(pairwise, aes(set_2, set_1, fill = pct_of_1)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f", pct_of_1)), size = 3) +
  facet_grid(CLASS ~ tumour) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100)) +
  scale_y_discrete(limits = rev) +
  labs(x = "found also in", y = "set (row = 100%)", fill = "% of row",
       title = "Shared PASS variants between settings") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.subtitle = element_text(size = 8))
save_plot(p_pair, "pairwise_shared", 19, 9)

# qcVCF pairwise plots (SNV, INDEL): every tumour x set against every other, so sharing between tumours shows too.
# qcVCF::plot_pairwise_shared_mutations() merges the table with itself without allow.cartesian: a variant in
# k columns gives k^2 rows, which data.table refuses for the SNVs here (18 columns, many shared variants).
# Fallback: the same numbers (100 * shared / N of the row column) from a sparse variant x column matrix,
# drawn like qcVCF's heatmap; columns ordered tumour, then set.
pairwise_shared_sparse <- function(q, cohort) {
  q <- unique(q[, .(ID = paste(CHROM, POS, REF, ALT, sep = ":"), sample_tool)])
  cols <- unique(q$sample_tool)
  cols <- cols[order(sub("^(HCC_[0-9]+)_.*$", "\\1", cols), match(sub("^HCC_[0-9]+_", "", cols), SET_LEVELS))]
  ids <- unique(q$ID)
  M <- Matrix::sparseMatrix(i = match(q$ID, ids), j = match(q$sample_tool, cols), x = 1, dims = c(length(ids), length(cols)),
                            dimnames = list(NULL, cols))
  S <- as.matrix(Matrix::crossprod(M))                      # S[x, y] = variants in both x and y
  d <- data.table(sampleID.x = rownames(S)[row(S)], sampleID.y = colnames(S)[col(S)], N.y = as.vector(S))
  d[, N.x := diag(S)[sampleID.x]]
  d[, prop_shared := 100 * N.y / N.x]
  d[, `:=`(sampleID.x = factor(sampleID.x, levels = cols), sampleID.y = factor(sampleID.y, levels = cols))]
  p <- ggplot(d, aes(x = sampleID.y, y = sampleID.x, fill = prop_shared)) +
    geom_tile(color = "black") +
    geom_text(aes(label = sprintf("%.1f", prop_shared)), size = 3.5) +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "% of shared\nmutations") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.title = element_text(size = rel(0.9)), legend.text = element_text(size = rel(0.9))) +
    labs(title = paste0("Pairwise comparison of shared mutations between samples (", cohort, ")"),
         subtitle = "only PASS variants before QC filtering", x = NULL, y = NULL,
         caption = "Proportion of mutations from sampleID.x shared with sampleID.y")
  list(data = d, plot = p)
}
qc_pairwise <- data.table()
if (HAVE_QCVCF) {
  qc_pairwise <- rbindlist(lapply(c("SNV", "INDEL"), function(cl) {
    q <- pass[CLASS == cl, .(CHROM, POS, REF, ALT, subjectID = tumour, tool = as.character(SET))]
    q[, sample_tool := paste(subjectID, tool, sep = "_")]
    pw <- tryCatch(qcVCF::plot_pairwise_shared_mutations(q, sample_col = "sample_tool", cohort = cl),
                   error = function(e) { message("qcVCF pairwise ", cl, " failed (", conditionMessage(e), ") - same plot from pairwise_shared_sparse()"); NULL })
    if (is.null(pw)) pw <- pairwise_shared_sparse(q, cohort = cl)
    save_plot(pw$plot, paste0("qcVCF_pairwise_shared_", cl), 17, 14)
    as.data.table(pw$data)[, CLASS := cl]
  }), fill = TRUE)
  print(qc_pairwise)
}

# 5b. Exclusive categories: in which sets each PASS variant is
memb <- dcast(pass, tumour + CLASS + ID ~ SET, fun.aggregate = length, value.var = "ID")
for (s in SET_LEVELS) if (!s %in% names(memb)) memb[, (s) := 0L]
memb[, N_SETS := rowSums(.SD > 0), .SDcols = SET_LEVELS]
memb[, SETS := apply(.SD > 0, 1, function(v) paste(SET_LEVELS[v], collapse = "+")), .SDcols = SET_LEVELS]
exclusive <- memb[, .N, by = .(tumour, CLASS, SETS, N_SETS)]
exclusive[, pct := round(100 * N / sum(N), 1), by = .(tumour, CLASS)]
setorder(exclusive, tumour, CLASS, -N)
print(exclusive)

lvl <- unlist(lapply(seq_along(SET_LEVELS), function(k) combn(SET_LEVELS, k, paste, collapse = "+")))
p_excl <- ggplot(exclusive, aes(factor(SETS, levels = lvl), N, fill = factor(N_SETS))) +
  geom_col() + geom_text(aes(label = N), vjust = -0.3, size = 2.6) +
  facet_grid(CLASS ~ tumour, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "found in (exclusive)", y = "PASS variants", fill = "sets",
       title = "Overlap between settings, exclusive categories (UpSet-style)",
       subtitle = "A_mutect2_strelka is a subset of A_mutect2 and A_strelka by definition") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), plot.subtitle = element_text(size = 8))
save_plot(p_excl, "overlap_exclusive", 22, 9)

# 5c. Proportion of each set's PASS variants also found in each other set, as stacked bars:
#     for each set, its variants split by how many OTHER sets also have them
share_n <- melt(memb, id.vars = c("tumour", "CLASS", "ID", "N_SETS"), measure.vars = SET_LEVELS,
                variable.name = "SET", value.name = "in_set")[in_set > 0]
share_n <- share_n[, .N, by = .(tumour, CLASS, SET, other_sets = N_SETS - 1L)]
share_n[, pct := 100 * N / sum(N), by = .(tumour, CLASS, SET)]
p_share <- ggplot(share_n, aes(SET, pct, fill = factor(other_sets))) +
  geom_col() +
  facet_grid(CLASS ~ tumour) +
  scale_fill_brewer(palette = "Blues", name = "also in\n# other sets") +
  labs(x = NULL, y = "% of the set's PASS variants", title = "Proportion of PASS variants shared with the other settings") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.subtitle = element_text(size = 8))
save_plot(p_share, "proportion_shared", 16, 7)

# ---- 6. SBS96 context ----------------------------------------------------------
SBS96_TYPES <- sort(CJ(b5 = c("A", "C", "G", "T"), sub = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                       b3 = c("A", "C", "G", "T"))[, paste0(b5, "[", sub, "]", b3)], method = "radix")
# + strand trinucleotide from the FASTA (samtools faidx); contig names matched to the .fai
trinuc_context <- function(d) {
  fai <- fread(paste0(FASTA, ".fai"), header = FALSE, select = 1)[[1]]
  chrom <- fifelse(d$CHROM %chin% fai, d$CHROM, sub("^chr", "", d$CHROM))
  reg <- paste0(chrom, ":", d$POS - 1L, "-", d$POS + 1L)
  rf <- tempfile(fileext = ".txt"); on.exit(unlink(rf))
  writeLines(unique(reg), rf)
  fa <- system2(SAMTOOLS, c("faidx", "-r", shQuote(rf), shQuote(FASTA)), stdout = TRUE)
  h <- startsWith(fa, ">")
  seqs <- data.table(region = sub("^>", "", fa[h]), NC_3 = toupper(fa[which(h) + 1L]))
  seqs[data.table(region = reg), on = "region"]$NC_3
}
# SBS96 class on the pyrimidine strand, e.g. A[C>T]G (G/A references reverse-complemented)
sbs96 <- function(ref, alt, nc3) {
  rc <- function(x) chartr("ACGT", "TGCA", x)
  pur <- ref %chin% c("A", "G")
  b5 <- fifelse(pur, rc(substr(nc3, 3, 3)), substr(nc3, 1, 1))
  b3 <- fifelse(pur, rc(substr(nc3, 1, 1)), substr(nc3, 3, 3))
  out <- paste0(b5, "[", fifelse(pur, rc(ref), ref), ">", fifelse(pur, rc(alt), alt), "]", b3)
  out[is.na(nc3)] <- NA_character_
  out
}
# Doublet half: another SNV of the same set (and same FILTER group) at the neighbouring position
in_doublet <- function(d) {
  here <- paste(d$tumour, d$SET, d$CHROM, d$POS)
  paste(d$tumour, d$SET, d$CHROM, d$POS - 1L) %chin% here | paste(d$tumour, d$SET, d$CHROM, d$POS + 1L) %chin% here
}

m96 <- data.table(); cosine <- data.table(); cosine_samples <- data.table(); sbs96_filter_n <- data.table()
if (HAVE_FASTA) {
  # All SNV records per set (PASS flag as defined per set); context looked up once per position
  a <- allset[CLASS == "SNV", .(tumour, SET, CHROM, POS, REF, ALT, PASS)]
  u <- unique(a[, .(CHROM, POS)])
  u[, NC_3 := trinuc_context(u)]
  a <- u[a, on = .(CHROM, POS)]
  rm(u)
  ref_mismatch <- a[!is.na(NC_3) & substr(NC_3, 2, 2) != REF, .N]
  if (ref_mismatch) warning(ref_mismatch, " SNVs whose REF is not the FASTA base - is this the genome the callers used?")
  # doublet halves: within the PASS set for the PASS profile, within all records for the qcVCF all/not PASS rows
  a[, IN_DBS_all := in_doublet(a)]
  a[(PASS), IN_DBS_pass := in_doublet(.SD)]
  a <- a[!is.na(NC_3) & substr(NC_3, 2, 2) == REF]
  a[, SBS96 := sbs96(REF, ALT, NC_3)]
  a[, `:=`(NC_3 = paste0(substr(SBS96, 1, 1), substr(SBS96, 3, 3), substr(SBS96, 7, 7)),   # pyrimidine strand
           filter = fifelse(PASS, "PASS", "not PASS"))]
  sbs96_filter_n <- a[, .(SNVs = .N, PASS = sum(PASS), PASS_doublet_halves = sum(PASS & IN_DBS_pass, na.rm = TRUE),
                          not_PASS = sum(!PASS)), by = .(tumour, SET)][order(tumour, SET)]
  print(sbs96_filter_n)
  p96set <- a[PASS & !IN_DBS_pass]

  # Count matrix: rows = 96 types, columns = tumour|set (SigProfiler format), PASS SNVs without doublet halves
  m96 <- dcast(p96set[, .N, by = .(SBS96, col = paste(tumour, SET, sep = "|"))], SBS96 ~ col, value.var = "N", fill = 0)
  m96 <- m96[data.table(SBS96 = SBS96_TYPES), on = "SBS96"]
  for (cl in setdiff(names(m96), "SBS96")) set(m96, which(is.na(m96[[cl]])), cl, 0L)
  setnames(m96, "SBS96", "MutationType")

  # 6b. ggplot per tumour: rows = set, PASS SNVs; % of each row's SNVs
  prof <- p96set[, .N, by = .(tumour, SET, SBS96)]
  prof <- prof[CJ(tumour = unique(prof$tumour), SET = factor(SET_LEVELS, levels = SET_LEVELS),
                  SBS96 = SBS96_TYPES, unique = TRUE), on = .(tumour, SET, SBS96)]
  prof[is.na(N), N := 0L][, SUB := substr(SBS96, 3, 5)]
  prof[, `:=`(n_row = sum(N), pct = if (sum(N)) 100 * N / sum(N) else 0), by = .(tumour, SET)]
  prof <- prof[n_row > 0]
  setorder(prof, tumour, SET)
  prof[, row := paste0(SET, "\nPASS (n=", n_row, ")")]
  prof[, row := factor(row, levels = unique(row))]
  prof[, SBS96 := factor(SBS96, levels = SBS96_TYPES[order(substr(SBS96_TYPES, 3, 5), SBS96_TYPES)])]
  sub_cols <- c("C>A" = "#1EBFF0", "C>G" = "#050708", "C>T" = "#E62725", "T>A" = "#CBCACB", "T>C" = "#A1CF64", "T>G" = "#EDC8C5")
  for (t in unique(prof$tumour)) {
    p96 <- ggplot(prof[tumour == t], aes(SBS96, pct, fill = SUB)) +
      geom_col(width = 0.8) +
      facet_grid(row ~ SUB, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = sub_cols, guide = "none") +
      labs(x = NULL, y = "% of SNVs", title = paste(t, "- SBS96 per setting, PASS SNVs (GRCm38 context)")) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5, family = "mono"),
                         strip.text.y = element_text(angle = 0), panel.spacing.x = unit(0.5, "mm"))
    save_plot(p96, paste0("sbs96_", t), 16, 2 + 1.4 * uniqueN(prof[tumour == t, row]))
  }

  # 6c. qcVCF 96-context plot per tumour: the sets stacked, rows all SNVs / PASS / not PASS (doublet halves
  #     of all records out)
  if (HAVE_QCVCF && all(sapply(c("cowplot", "stringr", "stringi", "ggtext"), requireNamespace, quietly = TRUE))) {
    for (t in names(TUMOURS)) {
      q96 <- a[tumour == t & !IN_DBS_all][order(SET), .(CHROM, POS, REF, ALT, NC_3, tool = as.character(SET), filter)]
      if (!nrow(q96)) next
      fig <- tryCatch(qcVCF::plot96_matrix(q96, rowsplit = "filter", plotsplitrow = "tool", orderplots = "no", showperc = "yes",
                                           dropempty = "yes", dontshowall = "no", roworder = "all SNVs,PASS,not PASS"),
                      error = function(e) { message("qcVCF plot96_matrix ", t, " skipped: ", conditionMessage(e)); NULL })
      rows <- q96[, uniqueN(filter) + 1L, by = tool]$V1
      for (k in seq_along(fig)) save_plot(fig[[k]], paste0("qcVCF_sbs96_", t, if (k > 1) paste0("_", k)),
                                          22 / 2.54, (sum(rows) * 4 + length(rows) * 3) / 2.54)
    }
  }
  rm(a)

  # 6d. Cosine similarity of the SBS96 profiles, per tumour, between sets
  cosine <- rbindlist(lapply(names(TUMOURS), function(t) {
    cols <- grep(paste0("^", t, "\\|"), names(m96), value = TRUE)
    if (length(cols) < 2) return(NULL)
    M <- as.matrix(m96[, ..cols]); colnames(M) <- sub("^[^|]*\\|", "", cols)
    nrm <- sqrt(colSums(M^2)); C <- crossprod(M) / outer(nrm, nrm)
    data.table(tumour = t, set_1 = rownames(C)[row(C)], set_2 = colnames(C)[col(C)], cosine = round(as.vector(C), 3))
  }))
  if (nrow(cosine)) {
    cosine[, `:=`(set_1 = factor(set_1, levels = SET_LEVELS), set_2 = factor(set_2, levels = SET_LEVELS))]
    print(dcast(cosine, tumour + set_1 ~ set_2, value.var = "cosine"))
    p_cos <- ggplot(cosine, aes(set_2, set_1, fill = cosine)) + geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.2f", cosine)), size = 3) + facet_wrap(~ tumour) +
      scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +
      scale_y_discrete(limits = rev) +
      labs(x = NULL, y = NULL, title = "Cosine similarity of SBS96 profiles between settings (PASS SNVs)") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.subtitle = element_text(size = 8))
    save_plot(p_cos, "sbs96_cosine", 19, 7)
  }

  # 6e. Cosine similarity per set, between tumours
  cosine_samples <- rbindlist(lapply(SET_LEVELS, function(s) {
    cols <- names(m96)[endsWith(names(m96), paste0("|", s))]
    if (length(cols) < 2) return(NULL)
    M <- as.matrix(m96[, ..cols]); colnames(M) <- sub("\\|.*$", "", cols)
    M <- M[, colSums(M) > 0, drop = FALSE]
    if (ncol(M) < 2) return(NULL)
    nrm <- sqrt(colSums(M^2)); C <- crossprod(M) / outer(nrm, nrm)
    data.table(set = s, tumour_1 = rownames(C)[row(C)], tumour_2 = colnames(C)[col(C)], cosine = round(as.vector(C), 3))
  }))
  if (nrow(cosine_samples)) {
    cosine_samples[, set := factor(set, levels = SET_LEVELS)]
    print(dcast(cosine_samples, set + tumour_1 ~ tumour_2, value.var = "cosine"))
    off <- cosine_samples[tumour_1 != tumour_2]
    cos_lo <- if (nrow(off)) min(0.95, floor(min(off$cosine) * 20) / 20) else 0
    cosine_samples[, fill := fifelse(tumour_1 == tumour_2, NA_real_, cosine)]
    cosine_samples[, txt_col := fifelse(!is.na(fill) & fill > cos_lo + 0.6 * (1 - cos_lo), "white", "black")]
    p_cos_s <- ggplot(cosine_samples, aes(tumour_2, tumour_1, fill = fill)) + geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.2f", cosine), colour = txt_col), size = 3) + scale_colour_identity() +
      facet_wrap(~ set, nrow = 1) +
      scale_fill_gradientn(colours = c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"), limits = c(cos_lo, 1),
                           oob = scales::squish, na.value = "grey85", name = "cosine") +
      labs(x = NULL, y = NULL, title = "Cosine similarity of SBS96 profiles between tumours, per setting (PASS SNVs)") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p_cos_s, "sbs96_cosine_between_tumours", 22, 4.5)
    cosine_samples[, c("fill", "txt_col") := NULL]
  }
}

# ---- 7. Outputs --------------------------------------------------------------
w <- function(x, name) if (nrow(x)) fwrite(x, file.path(OUT, paste0(name, ".tsv")), sep = "\t")
w(files[, !"path"], "files"); w(records_per_file, "records_per_file"); w(filter_values, "filter_values")
w(dropped_alleles, "dropped_non_ACGT_alleles")
w(counts, "counts_per_setting"); w(pairwise, "pairwise_shared"); w(qc_pairwise, "qcVCF_pairwise_shared")
w(exclusive, "overlap_exclusive"); w(share_n, "proportion_shared")
w(m96, "sbs96_matrix"); w(cosine, "sbs96_cosine"); w(cosine_samples, "sbs96_cosine_between_tumours")
w(sbs96_filter_n, "sbs96_snvs_by_filter")
fwrite(pass[order(tumour, SET, CHROM, POS), .(tumour, SET, CLASS, CHROM, POS, REF, ALT)],
       file.path(OUT, "pass_variants_per_setting.tsv.gz"), sep = "\t")

saveRDS(list(files = files, records_per_file = records_per_file, pass = pass, counts = counts, pairwise = pairwise,
             memb = memb, exclusive = exclusive, share_n = share_n, m96 = m96, cosine = cosine,
             cosine_samples = cosine_samples, params = list(FASTA = FASTA, SET_LEVELS = SET_LEVELS)),
        file.path(OUT, "settings.rds"))

fwrite(data.table(parameter = c("DIR_A", "DIR_B", "DIR_C", "OUT", "FASTA", "BCFTOOLS", "SAMTOOLS",
                                paste0("TUMOUR_", names(TUMOURS)), paste0("SET_", names(SET_LABELS))),
                  value = c(DIR_A, DIR_B, DIR_C, OUT, FASTA, BCFTOOLS, SAMTOOLS, TUMOURS, SET_LABELS)),
       file.path(OUT, "parameters.tsv"), sep = "\t")
fwrite(data.table(file = files$rel, md5 = unname(tools::md5sum(files$path))), file.path(OUT, "input_md5.tsv"), sep = "\t")
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
writeLines(c(R.version.string, bcf("--version")[1],
             if (file.exists(SAMTOOLS)) system2(SAMTOOLS, "--version", stdout = TRUE)[1] else "samtools not found",
             paste("data.table", packageVersion("data.table")), paste("ggplot2", packageVersion("ggplot2")),
             if (HAVE_QCVCF) paste("qcVCF", packageVersion("qcVCF")) else "qcVCF not installed",
             paste("FASTA", FASTA), paste("git", git_hash)), file.path(OUT, "versions.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

md_table <- function(x) {
  x <- as.data.table(x)
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
writeLines(c(
  "# mouse_wgs_snv_indel_settings summary", "",
  paste0("Run ", format(started, "%Y-%m-%d %H:%M"), ". PASS = FILTER of each file's own FILTER; MNVs split into SNVs."), "",
  "Sets:", "", paste0("- `", SET_LEVELS, "`: ", SET_LABELS[SET_LEVELS]), "",
  "## Samples per file", "", md_table(files[, .(tumour, src, samples)]), "",
  "## Counts per setting", "", md_table(counts), "",
  "## Pairwise shared PASS (row set = 100%)", "", md_table(pairwise[set_1 != set_2]), "",
  "## Exclusive categories", "", md_table(exclusive), "",
  if (nrow(cosine)) c("## SBS96 cosine similarity between settings", "",
                      md_table(dcast(cosine, tumour + set_1 ~ set_2, value.var = "cosine")), "",
                      if (nrow(cosine_samples)) c("## SBS96 cosine similarity between tumours, per setting", "",
                        md_table(dcast(cosine_samples, set + tumour_1 ~ tumour_2, value.var = "cosine")), "")) else
    "SBS96 skipped (no FASTA/samtools).",
  "## Notes", "",
  "- B and C are Mutect2 only; A_mutect2 is the same-caller comparison.",
  "- C has no matched normal: germline variants not in the PoN can remain.",
  "- Doublet halves are left out of SBS96 (counted as one DBS event).",
  "- Counts, not mutations per Mb (no BAMs, no callable territory)."),
  file.path(OUT, "summary.md"))

log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_snv_indel_settings"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Parameters: sets ", paste(SET_LEVELS, collapse = ", "), "; PASS from each file's own FILTER"),
  paste0("- Inputs:     ", DIR_A, " (FreeBayes, Mutect2, Strelka files), ", DIR_B, ", ", DIR_C, " (", nrow(files), " files); FASTA ", FASTA),
  paste0("- Outputs:    ", OUT),
  paste0("- Counts/QC:  records read ", nrow(raw), ", PASS variants (all sets) ", nrow(pass)),
  "", old_log), log_file)

message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
