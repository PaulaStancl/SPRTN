#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24
#
# mouse_wgs_snv_indel.R - Task 1: SNV/indel characterization of the three normal-reference settings.
#
#   A  tumour vs matched non-adjacent KI normal  (Strelka2 + Mutect2 + FreeBayes, merged, INFO/CALLER)
#   B  tumour vs joint WT normals                (Mutect2; exploratory)
#   C  tumour-only + WT panel of normals         (Mutect2)
#
# Blocks:
#   0  EDIT THIS BLOCK
#   1  setup, helpers
#   2  find the files (same layout as 00_input_inspection.R)
#   3  read records + tumour FORMAT (AD, DP, AF) + SnpEff ANN
#   4  counts per tumour x setting (PASS / non-PASS, SNV / MNV / INS / DEL, Ts/Tv, ins/del)
#   5  FILTER reasons
#   6  distributions: VAF, depth, alt depth, indel length, per chromosome (autosomes vs chrX/chrY)
#   7  setting comparison per tumour (SNVs and indels separately): exclusive categories,
#      UpSet-style plot, Jaccard, % of A recovered in B and C, PASS in one / filtered in another
#   8  qcVCF: mutation counts and pairwise sharing between settings
#   9  annotation: consequence classes, Sprtn, LoF, missense, cancer genes (if a list is given)
#   10 outputs, summary.md, RUN_LOG.md
#
# Rules kept here (see the project brief):
#   - Counts only, not mutations per Mb: there are no BAMs, so no callable territory.
#   - Agreement between settings shows robustness to the normal, not somatic origin.
#     A variant in all three settings is NOT thereby a real driver, and a real driver need
#     not be in all three (A drops it if the KI normal has reads for it, C if the PoN does).
#   - n = 3 mice: descriptive, per tumour.
#   - Indels were split (bcftools norm -m -any) but NOT left-aligned (no mm10 FASTA), so the
#     same indel may be written differently between callers/settings and miss an overlap.
#
# The setting comparison uses Mutect2 in all three settings by default (COMPARE_A_CALLER), so
# the differences come from the normal, not the caller. Setting A's other callers are counted
# in blocks 4-6.
#
# qcVCF functions used (as in scripts/analysis/01_sarek.R, same arguments):
#   plot_mutation_counts(), plot_pairwise_shared_mutations(), overlap_and_annotate_shared_mutations()
# Everything else is data.table / ggplot2.
#
# Env on ssi: r-mutation (data.table, ggplot2, qcVCF, bcftools on the PATH).
# Output: results/mouse_wgs/mouse_wgs_snv_indel/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
ROOT     <- file.path(REPO, "data/mouse/wgs")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_snv_indel")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_snv_indel.R")
# bcftools: RStudio often does not see the micromamba env's PATH, so Sys.which() can be "".
# Then the env's own bin/bcftools is used. Check with: system2(BCFTOOLS, "--version")
BCFTOOLS <- Sys.which("bcftools")
if (!nzchar(BCFTOOLS)) BCFTOOLS <- "/common/WORK/pstancl/envs/r-mutation/bin/bcftools"

TUMOURS <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1",
             HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
             HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")

TEST_CHROM       <- TRUE        # TRUE: chr19 only (quick test); FALSE: genome-wide
COMPARE_A_CALLER <- "mutect2"   # setting A records used in the comparison (block 7): "mutect2" or "any"
ENGINEERED_GENE  <- "Sprtn"     # mm10 chr8
CANCER_GENES     <- NA          # optional: path to a one-column file of mouse gene symbols (orthologs of
                                # cancer genes, e.g. COSMIC CGC 1:1 orthologs); NA = skip that part

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
HAVE_QCVCF <- requireNamespace("qcVCF", quietly = TRUE)
if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS, ": install it (micromamba install -c bioconda bcftools) or set BCFTOOLS in block 0")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()

# bcftools helper: returns stdout lines; on failure shows bcftools' own error message
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
# Save a ggplot as PDF + PNG
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  invisible(p)
}
SETTING_LAB <- c(A = "A: vs KI normal", B = "B: vs joint WT", C = "C: tumour-only + WT PoN")

# ---- 2. Find the files -------------------------------------------------------
files <- data.table(path = list.files(ROOT, recursive = TRUE, full.names = TRUE,
                                      pattern = "\\.(bcf|vcf|vcf\\.gz)$"))
files[, rel := substring(path, nchar(ROOT) + 2)]
files[, setting := fcase(startsWith(rel, "tumor_vs_normal/snv_indel/"), "A",
                         startsWith(rel, "tumor_vs_joint_wt/"),         "B",
                         startsWith(rel, "tumor_only_pon/"),            "C",
                         default = NA_character_)]
files <- files[!is.na(setting)]
files[, tumour := NA_character_]
for (k in names(TUMOURS)) files[grepl(TUMOURS[[k]], basename(path), fixed = TRUE), tumour := k]
files[, .(setting, tumour, file = basename(path))]
stopifnot(!anyNA(files$tumour))                       # every file must name its tumour
files[, .N, by = setting]                             # expect 3 per setting

# Header: which INFO / FORMAT fields and samples each file has
hdr <- setNames(lapply(files$path, function(f) bcf(c("view", "-h", shQuote(f)))), files$rel)
has_tag <- function(h, kind, id) any(startsWith(h, paste0("##", kind, "=<ID=", id, ",")))
files[, `:=`(
  has_CALLER = sapply(rel, function(r) has_tag(hdr[[r]], "INFO", "CALLER")),
  has_ANN    = sapply(rel, function(r) has_tag(hdr[[r]], "INFO", "ANN")),
  has_AD     = sapply(rel, function(r) has_tag(hdr[[r]], "FORMAT", "AD")),
  has_DP     = sapply(rel, function(r) has_tag(hdr[[r]], "FORMAT", "DP")),
  has_AF     = sapply(rel, function(r) has_tag(hdr[[r]], "FORMAT", "AF")))]
# Sample names from the #CHROM header line (columns 10 onwards)
files[, samples := sapply(rel, function(r) {
  s <- strsplit(tail(grep("^#CHROM", hdr[[r]], value = TRUE), 1), "\t", fixed = TRUE)[[1]]
  paste(s[-(1:9)], collapse = ",")
})]
# Tumour column: the sample named like the tumour; else Strelka's placeholder "TUMOR"
files[, tumour_col := mapply(function(s, t) {
  s <- strsplit(s, ",", fixed = TRUE)[[1]]
  hit <- s[s == TUMOURS[[t]]]
  if (!length(hit)) hit <- s[toupper(s) %chin% c("TUMOR", "TUMOUR")]
  if (length(hit)) hit[1] else NA_character_
}, samples, tumour)]
files[, .(setting, tumour, samples, tumour_col, has_CALLER, has_ANN, has_AD, has_DP, has_AF)]
stopifnot(!anyNA(files$tumour_col))

# chr19 is "19" (Ensembl/iGenomes) or "chr19" (UCSC)
TEST_CONTIG <- if (any(grepl("^##contig=<ID=chr19,", unlist(hdr)))) "chr19" else "19"
TEST_CONTIG

# ---- 3. Read records ---------------------------------------------------------
# One row per record: position, FILTER, CALLER, tumour AD/DP/AF, and the raw ANN string.
read_calls <- function(i) {
  f <- files[i]
  fmt <- paste0("%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER",
                if (f$has_CALLER) "\\t%INFO/CALLER",
                if (f$has_ANN) "\\t%INFO/ANN",
                "[", if (f$has_AD) "\\t%AD", if (f$has_DP) "\\t%DP", if (f$has_AF) "\\t%AF", "]\\n")
  cn <- c("CHROM", "POS", "REF", "ALT", "FILTER", if (f$has_CALLER) "CALLER", if (f$has_ANN) "ANN",
          if (f$has_AD) "AD", if (f$has_DP) "DP", if (f$has_AF) "AF")
  cmd <- paste(shQuote(BCFTOOLS), "query", if (TEST_CHROM) paste("-t", TEST_CONTIG),
               "-s", shQuote(f$tumour_col), "-f", shQuote(fmt), shQuote(f$path))
  x <- tryCatch(fread(cmd = cmd, sep = "\t", header = FALSE, col.names = cn, quote = "",
                      colClasses = list(character = setdiff(cn, "POS"), integer = "POS")),
                error = function(e) { message(f$rel, ": ", conditionMessage(e)); data.table() })
  if (!nrow(x)) return(NULL)
  if (!f$has_CALLER) x[, CALLER := "mutect2"]
  for (cl in c("ANN", "AD", "DP", "AF")) if (!cl %in% names(x)) x[, (cl) := NA_character_]
  x[, `:=`(setting = f$setting, tumour = f$tumour)]
  x
}
calls <- rbindlist(lapply(seq_len(nrow(files)), read_calls), use.names = TRUE)
calls[, .N, by = .(setting, tumour)]                  # row counts straight from the files

# Variant key, type, PASS
calls[, ID := paste(CHROM, POS, REF, ALT, sep = ":")]
calls[, TYPE := fcase(nchar(REF) == 1 & nchar(ALT) == 1, "SNV",
                      nchar(REF) == nchar(ALT),          "MNV",
                      nchar(REF) <  nchar(ALT),          "INS",
                      default =                          "DEL")]
calls[, CLASS := fifelse(TYPE %chin% c("INS", "DEL"), "indel", TYPE)]   # SNV / MNV / indel
calls[, PASS := FILTER == "PASS"]
calls[, CHROM_CLASS := fcase(sub("^chr", "", CHROM) %chin% as.character(1:19), "autosome",
                             sub("^chr", "", CHROM) == "X", "chrX",
                             sub("^chr", "", CHROM) == "Y", "chrY",
                             default = "other")]
# Not left-aligned / complex records: REF and ALT longer than 1 and different first base
calls[, COMPLEX := TYPE %chin% c("INS", "DEL") & substr(REF, 1, 1) != substr(ALT, 1, 1)]

# Tumour depth, alt reads, VAF. AD after the split is "ref,alt". VAF from AD when possible
# (same definition for all callers), else FORMAT AF. Strelka writes neither (NA here).
num1 <- function(x, k) suppressWarnings(as.numeric(vapply(strsplit(x, ",", fixed = TRUE),
                                                          function(v) if (length(v) >= k) v[k] else NA_character_, "")))
calls[, `:=`(REF_READS = num1(AD, 1), ALT_READS = num1(AD, 2), DEPTH = suppressWarnings(as.numeric(DP)))]
calls[, VAF := fifelse(!is.na(ALT_READS) & (REF_READS + ALT_READS) > 0,
                       ALT_READS / (REF_READS + ALT_READS), num1(AF, 1))]
calls[, .(records = .N, with_VAF = sum(!is.na(VAF)), with_DP = sum(!is.na(DEPTH))), by = .(setting, CALLER)]

# Setting A: one record may carry several callers ("mutect2,strelka"). Long form, one row per caller:
calls_long <- calls[, .(CALLER1 = unlist(strsplit(CALLER, ",", fixed = TRUE))), by = .(setting, tumour, ID)]
calls_long[, .N, by = .(setting, CALLER1)]

# ---- 4. Counts per tumour x setting -------------------------------------------
# Ts/Tv on PASS SNVs: transitions are A<>G and C<>T
calls[TYPE == "SNV", TS := paste0(REF, ALT) %chin% c("AG", "GA", "CT", "TC")]
counts <- calls[, .(records = .N, PASS = sum(PASS), nonPASS = sum(!PASS),
                    SNV = sum(PASS & TYPE == "SNV"), MNV = sum(PASS & TYPE == "MNV"),
                    INS = sum(PASS & TYPE == "INS"), DEL = sum(PASS & TYPE == "DEL"),
                    complex_indel = sum(PASS & COMPLEX),
                    TsTv = round(sum(PASS & TYPE == "SNV" & TS, na.rm = TRUE) /
                                 sum(PASS & TYPE == "SNV" & !TS, na.rm = TRUE), 2),
                    ins_del = round(sum(PASS & TYPE == "INS") / sum(PASS & TYPE == "DEL"), 2),
                    PASS_autosome = sum(PASS & CHROM_CLASS == "autosome"),
                    PASS_chrX = sum(PASS & CHROM_CLASS == "chrX"), PASS_chrY = sum(PASS & CHROM_CLASS == "chrY")),
                by = .(setting, tumour)]
setorder(counts, tumour, setting)
counts

# Setting A per caller (a record counts for each caller that made it)
counts_A_caller <- calls_long[setting == "A"][calls[setting == "A", .(ID, tumour, PASS, TYPE)], on = .(tumour, ID), allow.cartesian = TRUE][
  , .(records = .N, PASS = sum(PASS), SNV = sum(PASS & TYPE == "SNV"), indel = sum(PASS & TYPE %chin% c("INS", "DEL"))),
  by = .(tumour, CALLER1)]
setorder(counts_A_caller, tumour, CALLER1)
counts_A_caller

p_counts <- ggplot(melt(counts, id.vars = c("setting", "tumour"), measure.vars = c("SNV", "MNV", "INS", "DEL"),
                        variable.name = "type", value.name = "n"),
                   aes(setting, n, fill = type)) +
  geom_col() + facet_wrap(~ tumour) +
  labs(x = NULL, y = "PASS records", fill = NULL,
       title = paste("PASS SNV / MNV / indel counts", if (TEST_CHROM) paste0("(", TEST_CONTIG, ")") else ""),
       caption = paste(paste(names(SETTING_LAB), SETTING_LAB, sep = " = "), collapse = "; ")) +
  theme_bw()
save_plot(p_counts, "counts_per_setting", 9, 5)

# ---- 5. FILTER reasons --------------------------------------------------------
filter_reasons <- calls[, .(reason = unlist(strsplit(FILTER, ";", fixed = TRUE))), by = .(setting, tumour, CALLER)][
  , .N, by = .(setting, tumour, CALLER, reason)]
setorder(filter_reasons, setting, tumour, -N)
filter_reasons

p_filter <- ggplot(filter_reasons[reason != "PASS"], aes(N, reorder(reason, N), fill = setting)) +
  geom_col(position = "dodge") + facet_wrap(~ tumour, scales = "free_x") +
  labs(x = "records", y = NULL, title = "FILTER reasons (non-PASS records)") + theme_bw()
save_plot(p_filter, "filter_reasons", 11, 7)

# ---- 6. Distributions --------------------------------------------------------
pass <- calls[PASS == TRUE]

# 6a. VAF, autosomes and chrX separately (males: clonal chrX variants sit near purity, not purity/2)
p_vaf <- ggplot(pass[!is.na(VAF) & CHROM_CLASS %chin% c("autosome", "chrX")],
                aes(VAF, colour = setting)) +
  geom_density(adjust = 0.7) + facet_grid(CHROM_CLASS ~ tumour, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Tumour VAF of PASS calls", subtitle = "from FORMAT AD (else AF); Strelka-only records have none") +
  theme_bw()
save_plot(p_vaf, "vaf_density", 11, 6)
vaf_summary <- pass[!is.na(VAF), .(n = .N, median_VAF = round(median(VAF), 3),
                                   q25 = round(quantile(VAF, 0.25), 3), q75 = round(quantile(VAF, 0.75), 3)),
                    by = .(setting, tumour, CHROM_CLASS, CLASS)]
setorder(vaf_summary, tumour, setting, CHROM_CLASS, CLASS)
vaf_summary

# 6b. Depth and alt reads
p_depth <- ggplot(melt(pass[!is.na(DEPTH)], id.vars = c("setting", "tumour"),
                       measure.vars = c("DEPTH", "ALT_READS"), variable.name = "metric"),
                  aes(setting, value, fill = setting)) +
  geom_boxplot(outlier.size = 0.3) + facet_grid(metric ~ tumour, scales = "free_y") +
  scale_y_log10() + labs(x = NULL, y = "reads (log10)", title = "Tumour depth and alt reads, PASS calls") +
  theme_bw() + theme(legend.position = "none")
save_plot(p_depth, "depth_altreads", 9, 6)
depth_summary <- pass[, .(median_DP = median(DEPTH, na.rm = TRUE), median_alt = median(ALT_READS, na.rm = TRUE)),
                      by = .(setting, tumour)]
depth_summary

# 6c. Indel length (ALT minus REF length: + insertion, - deletion)
pass[CLASS == "indel", INDEL_LEN := nchar(ALT) - nchar(REF)]
p_indel <- ggplot(pass[CLASS == "indel"], aes(pmax(pmin(INDEL_LEN, 20), -20), fill = setting)) +
  geom_bar(position = "dodge") + facet_wrap(~ tumour, ncol = 1, scales = "free_y") +
  labs(x = "indel length (bp; capped at +/-20)", y = "PASS indels", title = "Indel length") + theme_bw()
save_plot(p_indel, "indel_length", 9, 7)

# 6d. Per chromosome
chrom_counts <- pass[, .N, by = .(setting, tumour, CHROM, CLASS)]
chrom_counts[, CHROM := factor(CHROM, levels = unique(CHROM[order(suppressWarnings(as.integer(sub("^chr", "", CHROM))), CHROM)]))]
p_chrom <- ggplot(chrom_counts, aes(CHROM, N, fill = setting)) +
  geom_col(position = "dodge") + facet_grid(CLASS ~ tumour, scales = "free_y") +
  labs(x = NULL, y = "PASS records", title = "PASS calls per chromosome") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
save_plot(p_chrom, "per_chromosome", 14, 7)

# ---- 7. Setting comparison ---------------------------------------------------
# Unit: chrom:pos:ref:alt, PASS only, per tumour; SNVs and indels separately (MNVs reported, not compared).
# Setting A: Mutect2 records only (COMPARE_A_CALLER = "mutect2"), so all three settings share the caller.
cmp_calls <- calls[CLASS %chin% c("SNV", "indel")]
if (COMPARE_A_CALLER != "any")
  cmp_calls <- cmp_calls[setting != "A" | grepl(COMPARE_A_CALLER, CALLER, fixed = TRUE)]
cmp <- cmp_calls[PASS == TRUE, .(ID, setting, tumour, CLASS)]
cmp <- unique(cmp)
cmp[, .N, by = .(tumour, CLASS, setting)][order(tumour, CLASS, setting)]

# 7a. Exclusive categories: A only, B only, C only, A+B, A+C, B+C, A+B+C
memb <- dcast(cmp, tumour + CLASS + ID ~ setting, fun.aggregate = length, value.var = "ID")
for (s in c("A", "B", "C")) if (!s %in% names(memb)) memb[, (s) := 0L]
memb[, category := paste0(fifelse(A > 0, "A", ""), fifelse(B > 0, "B", ""), fifelse(C > 0, "C", ""))]
memb[, category := gsub("(?<=.)(?=.)", "+", category, perl = TRUE)]
overlap <- memb[, .N, by = .(tumour, CLASS, category)]
overlap[, pct := round(100 * N / sum(N), 1), by = .(tumour, CLASS)]
setorder(overlap, tumour, CLASS, -N)
overlap

# UpSet-style plot: one bar per exclusive category
lvl <- c("A", "B", "C", "A+B", "A+C", "B+C", "A+B+C")
p_upset <- ggplot(overlap, aes(factor(category, levels = lvl), N)) +
  geom_col(fill = "grey30") + geom_text(aes(label = N), vjust = -0.3, size = 3) +
  facet_grid(CLASS ~ tumour, scales = "free_y") +
  labs(x = "found in setting(s), exclusive", y = "PASS variants",
       title = "Overlap between settings (exclusive categories)",
       subtitle = paste0("A = ", COMPARE_A_CALLER, " records of setting A. Agreement = robust to the normal, not proof of somatic origin.")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_upset, "overlap_upset", 11, 6)

# 7b. Jaccard between each pair of settings, and % of A recovered in B and C
pairs <- CJ(s1 = c("A", "B", "C"), s2 = c("A", "B", "C"))[s1 < s2]
jaccard <- memb[, rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  a <- .SD[[pairs$s1[i]]] > 0; b <- .SD[[pairs$s2[i]]] > 0
  data.table(s1 = pairs$s1[i], s2 = pairs$s2[i], n1 = sum(a), n2 = sum(b), shared = sum(a & b),
             jaccard = round(sum(a & b) / sum(a | b), 3))
})), by = .(tumour, CLASS), .SDcols = c("A", "B", "C")]
jaccard
recovery <- memb[A > 0, .(n_A = .N, pct_A_in_B = round(100 * mean(B > 0), 1),
                          pct_A_in_C = round(100 * mean(C > 0), 1),
                          pct_A_in_B_and_C = round(100 * mean(B > 0 & C > 0), 1)), by = .(tumour, CLASS)]
recovery

p_jac <- ggplot(jaccard, aes(s1, s2, fill = jaccard)) + geom_tile() +
  geom_text(aes(label = sprintf("%.2f\n(%d)", jaccard, shared)), size = 3) +
  facet_grid(CLASS ~ tumour) + scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
  labs(x = NULL, y = NULL, title = "Jaccard index between settings (shared PASS variants)") + theme_bw()
save_plot(p_jac, "jaccard", 9, 6)

# 7c. PASS in one setting but present and filtered in another (needs non-PASS records in the files)
all_filters <- unique(cmp_calls[, .(tumour, CLASS, ID, setting, FILTER)])
all_filters <- all_filters[, .(FILTER = paste(unique(FILTER), collapse = "|")), by = .(tumour, CLASS, ID, setting)]
fw <- dcast(all_filters, tumour + CLASS + ID ~ setting, value.var = "FILTER")
for (s in c("A", "B", "C")) if (!s %in% names(fw)) fw[, (s) := NA_character_]
state <- function(x) fcase(is.na(x), "absent", grepl("(^|\\|)PASS($|\\|)", x), "PASS", default = "filtered")
fw[, `:=`(A_state = state(A), B_state = state(B), C_state = state(C))]
pass_vs_filtered <- fw[(A_state == "PASS" | B_state == "PASS" | C_state == "PASS") &
                         (A_state == "filtered" | B_state == "filtered" | C_state == "filtered"),
                       .N, by = .(tumour, CLASS, A_state, B_state, C_state)]
setorder(pass_vs_filtered, tumour, CLASS, -N)
pass_vs_filtered
# Which FILTER removed a call that another setting kept?
why_filtered <- melt(fw[(A_state == "PASS" | B_state == "PASS" | C_state == "PASS")],
                     id.vars = c("tumour", "CLASS", "ID"), measure.vars = c("A", "B", "C"),
                     variable.name = "setting", value.name = "FILTER")[!is.na(FILTER) & !grepl("PASS", FILTER)]
why_filtered <- why_filtered[, .(reason = unlist(strsplit(FILTER, "[;|]"))), by = .(tumour, CLASS, setting)][, .N, by = .(tumour, CLASS, setting, reason)]
setorder(why_filtered, tumour, CLASS, setting, -N)
why_filtered

# ---- 8. qcVCF ----------------------------------------------------------------
# Settings play the role of "tool" in qcVCF (one tumour at a time). mutType: SNV, MNV, INS, DEL.
qc_out <- list()
if (HAVE_QCVCF) {
  library(qcVCF)
  to_qc <- function(x) x[, .(CHROM, POS, REF, ALT, subjectID = tumour, tool = setting, mutType = TYPE)]
  qc_pass <- if (COMPARE_A_CALLER != "any") pass[setting != "A" | grepl(COMPARE_A_CALLER, CALLER, fixed = TRUE)] else pass
  for (t in names(TUMOURS)) {
    q <- to_qc(qc_pass[tumour == t])
    if (!nrow(q)) next
    save_plot(plot_mutation_counts(q, grp_col = "tool", cohort = t), paste0("qcVCF_mutation_counts_", t), 8, 6)
    q[, varClass := fifelse(mutType %chin% c("INS", "DEL"), "INDEL", mutType)]
    for (cl in c("SNV", "INDEL")) {
      # qcVCF's by = get() fails on very small classes; skip those instead of stopping
      pw <- tryCatch(plot_pairwise_shared_mutations(q[varClass == cl], sample_col = "tool", cohort = paste(t, cl)),
                     error = function(e) { message("qcVCF pairwise ", t, " ", cl, " skipped: ", conditionMessage(e)); NULL })
      if (is.null(pw)) next
      save_plot(pw$plot, paste0("qcVCF_pairwise_shared_", t, "_", cl), 7, 6)
      qc_out[[paste(t, cl)]] <- pw$data[, `:=`(tumour = t, varClass = cl)]
    }
  }
  qc_pairwise <- rbindlist(qc_out, fill = TRUE)
  print(qc_pairwise)
} else message("qcVCF not installed in this env: block 8 skipped")

# ---- 9. Annotation -----------------------------------------------------------
# SnpEff ANN: "Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Biotype|Rank|HGVS.c|HGVS.p|..."
# After the split each record still lists every allele of the original record, so keep only the
# entries whose Allele equals this record's ALT.
# Transcript rule: SnpEff writes entries most severe first, so the first entry for the allele is the
# "worst consequence" call (gene and transcript taken from that entry). Stated, not tuned.
ann <- pass[!is.na(ANN) & ANN != ".", .(setting, tumour, ID, ALT, CALLER, TYPE, VAF, ANN)]
if (nrow(ann)) {
  ann_long <- ann[, .(entry = unlist(strsplit(ANN, ",", fixed = TRUE))), by = .(setting, tumour, ID, ALT, CALLER, TYPE, VAF)]
  ann_long[, f := strsplit(entry, "|", fixed = TRUE)]
  ann_long[, `:=`(Allele = vapply(f, `[`, "", 1), Consequence = vapply(f, `[`, "", 2), Impact = vapply(f, `[`, "", 3),
                  Gene = vapply(f, `[`, "", 4), Transcript = vapply(f, `[`, "", 7),
                  HGVSc = vapply(f, function(v) if (length(v) >= 10) v[10] else NA_character_, ""),
                  HGVSp = vapply(f, function(v) if (length(v) >= 11) v[11] else NA_character_, ""))]
  ann_long[, c("f", "entry") := NULL]
  ann_long <- ann_long[Allele == ALT]
  ann_top <- ann_long[, .SD[1], by = .(setting, tumour, ID)]          # first (most severe) entry per variant

  # Consequence classes (first term of SnpEff's "a&b" consequences)
  ann_top[, CONS := sub("&.*", "", Consequence)]
  ann_top[, CONS_CLASS := fcase(
    CONS %chin% c("stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
                  "start_lost", "stop_lost", "exon_loss_variant"), "LoF",
    CONS == "missense_variant", "missense",
    CONS %chin% c("conservative_inframe_insertion", "conservative_inframe_deletion",
                  "disruptive_inframe_insertion", "disruptive_inframe_deletion"), "inframe_indel",
    CONS %chin% c("synonymous_variant", "stop_retained_variant", "start_retained_variant"), "synonymous",
    CONS %chin% c("splice_region_variant"), "splice_region",
    CONS %chin% c("5_prime_UTR_variant", "3_prime_UTR_variant", "5_prime_UTR_premature_start_codon_gain_variant"), "UTR",
    CONS %chin% c("intron_variant"), "intronic",
    CONS %chin% c("upstream_gene_variant", "downstream_gene_variant", "intergenic_region"), "intergenic/flanking",
    default = "other non-coding")]
  cons_counts <- ann_top[, .N, by = .(tumour, setting, CONS_CLASS)]
  setorder(cons_counts, tumour, setting, -N)
  print(dcast(cons_counts, tumour + CONS_CLASS ~ setting, value.var = "N", fill = 0))

  # Cancer gene list (optional, mouse symbols)
  cancer_genes <- if (!is.na(CANCER_GENES) && file.exists(CANCER_GENES)) fread(CANCER_GENES, header = FALSE)[[1]] else character()

  # Separate groups. No damage scores in SnpEff ANN, so missense stays "missense (no damage score)".
  # None of these is called a driver here.
  ann_top[, GROUP := fcase(
    Gene == ENGINEERED_GENE,                                "engineered gene (Sprtn)",
    CONS_CLASS == "LoF",                                    "loss-of-function",
    CONS_CLASS == "missense",                               "missense (no damage score)",
    CONS_CLASS %chin% c("inframe_indel", "splice_region"),  "uncertain protein effect",
    default = "other")]
  ann_top[, CANCER_GENE := Gene %chin% cancer_genes]
  protein_altering <- ann_top[GROUP != "other"]
  setorder(protein_altering, tumour, GROUP, Gene)
  print(protein_altering[, .N, by = .(tumour, setting, GROUP, CANCER_GENE)])

  # Variants in the engineered gene, every setting (expect the construct to show up here)
  sprtn <- ann_long[Gene == ENGINEERED_GENE]
  print(sprtn)
} else {
  message("No SnpEff ANN in these files: block 9 skipped")
  ann_top <- protein_altering <- sprtn <- cons_counts <- data.table()
}

# ---- 10. Outputs -------------------------------------------------------------
tag <- if (TEST_CHROM) paste0("_", TEST_CONTIG) else ""
w <- function(x, name) fwrite(x, file.path(OUT, paste0(name, tag, ".tsv")), sep = "\t")
w(files[, !"path"], "files")
w(counts, "counts"); w(counts_A_caller, "counts_settingA_per_caller")
w(filter_reasons, "filter_reasons"); w(vaf_summary, "vaf_summary"); w(depth_summary, "depth_summary")
w(chrom_counts, "per_chromosome")
w(overlap, "overlap_exclusive_categories"); w(jaccard, "jaccard"); w(recovery, "recovery_of_A")
w(pass_vs_filtered, "pass_in_one_filtered_in_other"); w(why_filtered, "why_filtered")
w(memb, "variant_setting_membership")
if (exists("qc_pairwise")) w(qc_pairwise, "qcVCF_pairwise_shared")
if (nrow(ann_top)) { w(cons_counts, "consequence_counts"); w(protein_altering, "protein_altering"); w(sprtn, "engineered_gene_variants") }

saveRDS(list(files = files, calls = calls, calls_long = calls_long, counts = counts, counts_A_caller = counts_A_caller,
             filter_reasons = filter_reasons, vaf_summary = vaf_summary, depth_summary = depth_summary,
             chrom_counts = chrom_counts, memb = memb, overlap = overlap, jaccard = jaccard, recovery = recovery,
             fw = fw, pass_vs_filtered = pass_vs_filtered, why_filtered = why_filtered,
             ann_top = ann_top, protein_altering = protein_altering, sprtn = sprtn,
             params = list(TEST_CHROM = TEST_CHROM, TEST_CONTIG = TEST_CONTIG, COMPARE_A_CALLER = COMPARE_A_CALLER)),
        file.path(OUT, paste0("snv_indel", tag, ".rds")))

params <- data.table(parameter = c("ROOT", "OUT", "BCFTOOLS", "TEST_CHROM", "TEST_CONTIG", "COMPARE_A_CALLER",
                                   "ENGINEERED_GENE", "CANCER_GENES", paste0("TUMOUR_", names(TUMOURS))),
                     value = c(ROOT, OUT, BCFTOOLS, TEST_CHROM, TEST_CONTIG, COMPARE_A_CALLER,
                               ENGINEERED_GENE, CANCER_GENES, TUMOURS))
fwrite(params, file.path(OUT, "parameters.tsv"), sep = "\t")
fwrite(data.table(file = files$rel, md5 = unname(tools::md5sum(files$path))), file.path(OUT, "input_md5.tsv"), sep = "\t")
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
writeLines(c(R.version.string, bcf("--version")[1],
             paste("data.table", packageVersion("data.table")), paste("ggplot2", packageVersion("ggplot2")),
             if (HAVE_QCVCF) paste("qcVCF", packageVersion("qcVCF")) else "qcVCF not installed",
             paste("git", git_hash)), file.path(OUT, "versions.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

md_table <- function(x) {
  x <- as.data.table(x)
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
writeLines(c(
  "# mouse_wgs_snv_indel summary", "",
  paste0("Run ", format(started, "%Y-%m-%d %H:%M"), "; ", if (TEST_CHROM) paste0("**", TEST_CONTIG, " only (test run)**") else "genome-wide",
         "; setting A in the comparison = ", COMPARE_A_CALLER, " records."), "",
  "## Counts per tumour and setting (PASS unless stated)", "", md_table(counts), "",
  "## Setting A per caller", "", md_table(counts_A_caller), "",
  "## % of A recovered in B and C", "", md_table(recovery), "",
  "## Exclusive overlap categories", "", md_table(overlap), "",
  "## Notes", "",
  "- Counts, not mutations per Mb (no BAMs, no callable territory).",
  "- Overlap shows robustness to the choice of normal, not somatic origin; B is exploratory (multi-normal mode with unrelated mice).",
  "- Indels not left-aligned (no FASTA): some indel mismatches between settings may be representation only.",
  "- Strelka records carry no AD/AF, so their VAF is missing.",
  "- Annotation groups are not driver calls; missense has no damage score in SnpEff ANN."),
  file.path(OUT, "summary.md"))

log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_snv_indel"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Parameters: TEST_CHROM=", TEST_CHROM, " (", TEST_CONTIG, "), COMPARE_A_CALLER=", COMPARE_A_CALLER,
         ", CANCER_GENES=", CANCER_GENES),
  paste0("- Inputs:     ", ROOT, " (", nrow(files), " BCFs; split with bcftools norm -m -any by the provider, not left-aligned)"),
  paste0("- Outputs:    ", OUT),
  paste0("- Counts/QC:  records ", nrow(calls), ", PASS ", nrow(pass)),
  "", old_log), log_file)

message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
