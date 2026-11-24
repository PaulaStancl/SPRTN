#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24
#
# mouse_wgs_example_sites.R - diagnostics for two SNVs where the callers disagree and the merged
# (union) file's consensus misleads:
#   chr19:30419577 T>C   HCC_1: Mutect2 PASS, FreeBayes PASS, Strelka2 LowEVS -> merged LowEVS
#   chr10:3184122  G>T   HCC_1: Strelka2 PASS, FreeBayes "." -> merged PASS with "strelka,freebayes"
#
# Everything comes from the ORIGINAL single-caller files (not the merged FILTER/INFO):
#   setting A  tumor_vs_normal/snv_indel/  Mutect2, Strelka2, FreeBayes (tumour + matched KI normal)
#   setting B  tumor_vs_joint_wt/          Mutect2 joint calls: the only files with the 3 WT normals
#
# Blocks:
#   0  EDIT THIS BLOCK (sites, paths)
#   1  setup, helpers
#   2  files
#   3  read the records at the two sites from every file (all samples, all INFO/FORMAT)
#   4  per-sample metrics: GT, DP, REF/ALT reads, VAF, GQ, strand ALT reads
#   5  per-record metrics: FILTER, QUAL, MQ, BQ, strand, Mutect2 TLOD/NLOD, Strelka2 SomaticEVS, FreeBayes MQM ...
#   6  FreeBayes: Patrik's expression, each criterion separately
#   7  read counts in all nine samples, per caller, from VCF FORMAT (no BAMs - see note)
#   8  plots 1-8
#   9  outputs, RUN_LOG.md
#
# Notes
# - Strelka2 has no GT; its REF/ALT counts are the tier-1 {BASE}U counts. Its strand evidence is
#   SNVSB (a strand-bias score), not per-strand counts.
# - FreeBayes INFO (AO, SAF, SAR, MQM, QUAL) is pooled over tumour + normal, so QUAL/AO and SAF/SAR in
#   Patrik's expression are not tumour-only values.
# - Nine-sample read counts: there are no BAMs, so they are VCF FORMAT counts, shown PER CALLER:
#   tumour and KI normal from each setting-A caller file, WT normals (and tumours) from the setting-B
#   Mutect2 joint file. These are the callers' own (filtered, realigned) read counts, not a uniform
#   pileup. A caller without a record at the site is "absent", not 0 reads.
#
# Env on ssi: r-mutation (data.table, ggplot2, bcftools).
# Output: results/mouse_wgs/mouse_wgs_example_sites/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
ROOT     <- file.path(REPO, "data/mouse/wgs")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_example_sites")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_example_sites.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")

SITES <- data.table::data.table(
  variant = c("chr19:30419577 T>C", "chr10:3184122 G>T"),
  CHROM   = c("19", "10"),
  POS     = c(30419577L, 3184122L),
  REF     = c("T", "G"),
  ALT     = c("C", "T"))

TUMOURS    <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1", HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
                HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")
CALLERS    <- c("strelka", "mutect2", "freebayes")

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# Console as wide as the terminal (falls back to 80 when there is no terminal, e.g. batch Rscript)
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()

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
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  invisible(p)
}
# INFO string -> named list (flags = "TRUE"); the long VEP CSQ / SnpEff ANN are dropped
kv <- function(s) {
  p <- strsplit(s, ";", fixed = TRUE)[[1]]
  key <- sub("=.*$", "", p); val <- ifelse(grepl("=", p, fixed = TRUE), sub("^[^=]*=", "", p), "TRUE")
  keep <- !key %chin% c("CSQ", "ANN")
  as.list(setNames(val[keep], key[keep]))
}
# i-th comma-separated value as number (NA if missing)
pick <- function(s, i) {
  if (is.null(s) || is.na(s) || is.na(i)) return(NA_real_)
  suppressWarnings(as.numeric(strsplit(s, ",", fixed = TRUE)[[1]][i]))
}
# Short sample name and class
short_name <- function(s) sub("^(id_[0-9]+_)?(sfb_dikic_)?([0-9]{4}_[0-9]{2}_[0-9]{2}_)?", "", s)
sample_class <- function(s) fcase(grepl("KI_HCC_[0-9]|^TUMOU?R$", s, ignore.case = TRUE), "tumour",
                                  grepl("KI[0-9]$|^NORMAL$", s), "GM normal",
                                  default = "WT normal")

CALLER_COLS <- c(strelka = "#4C72B0", mutect2 = "#DD8452", freebayes = "#55A868")
STATUS_COLS <- c(PASS = "#2E7D32", filtered = "#C62828", unassessed = "#BDBDBD", absent = "white")
CLASS_LEVELS <- c("tumour", "GM normal", "WT normal")

# ---- 2. Files ----------------------------------------------------------------
files <- data.table(path = list.files(ROOT, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf|vcf\\.gz)$"))
files[, rel := substring(path, nchar(ROOT) + 2)]
files[, setting := fcase(startsWith(rel, "tumor_vs_normal/snv_indel/"), "A",
                         startsWith(rel, "tumor_vs_joint_wt/"), "B", default = NA_character_)]
files <- files[!is.na(setting) & !grepl("union|merged", rel, ignore.case = TRUE)]   # single-caller files only
files[, tumour := NA_character_]
for (k in names(TUMOURS)) files[grepl(TUMOURS[[k]], rel, fixed = TRUE), tumour := k]
files[, caller := NA_character_]
for (cl in CALLERS) files[is.na(caller) & grepl(cl, rel, ignore.case = TRUE), caller := cl]
files[setting == "B" & is.na(caller), caller := "mutect2"]
files[, .(setting, tumour, caller, rel)]
stopifnot(!anyNA(files$tumour), !anyNA(files$caller))

# ---- 3. Records at the two sites ---------------------------------------------
# One bcftools call per file for both sites (-t streams the file, no index needed).
REGIONS <- paste(SITES[, paste0(CHROM, ":", POS)], collapse = ",")
read_records <- function(i) {
  f <- files[i]
  v <- bcf(c("view", "-H", "-t", REGIONS, shQuote(f$path)))
  if (!length(v)) return(NULL)
  smp <- bcf(c("query", "-l", shQuote(f$path)))
  rbindlist(lapply(v, function(l) {
    x <- strsplit(l, "\t", fixed = TRUE)[[1]]
    s <- SITES[CHROM == x[1] & POS == as.integer(x[2]) & REF == x[4]]
    if (!nrow(s)) return(NULL)
    k <- match(s$ALT, strsplit(x[5], ",", fixed = TRUE)[[1]])
    if (is.na(k)) return(NULL)
    keys <- strsplit(x[9], ":", fixed = TRUE)[[1]]
    info <- kv(x[8])
    rbindlist(lapply(seq_along(smp), function(j) {
      vals <- strsplit(x[9 + j], ":", fixed = TRUE)[[1]]
      data.table(variant = s$variant, setting = f$setting, tumour = f$tumour, caller = f$caller, rel = f$rel,
                 REF = s$REF, ALT = s$ALT, k = k, QUAL = x[6], FILTER = x[7], sample_raw = smp[j],
                 info = list(info), fmt = list(as.list(setNames(vals[seq_along(keys)], keys))))
    }))
  }))
}
rec <- rbindlist(lapply(seq_len(nrow(files)), read_records))
# Short names: Strelka's NORMAL/TUMOR placeholders become the file's KI normal / tumour
rec[, class := sample_class(sample_raw)]
rec[, sample := fcase(sample_raw %chin% c("TUMOR", "TUMOUR"), tumour, sample_raw == "NORMAL", KI_NORMALS[tumour],
                      class == "tumour", sub(".*KI_(HCC_[0-9]).*", "\\1", sample_raw), default = short_name(sample_raw))]
rec[, .N, by = .(variant, setting, tumour, caller)][order(variant, setting, tumour, caller)]   # records found
rec[, .(variant, setting, tumour, caller, FILTER, sample_raw, sample, class)]

# ---- 4. Per-sample metrics ---------------------------------------------------
sample_metrics <- function(fm, cl, ref, alt, k) {
  g <- function(key) if (!is.null(fm[[key]])) fm[[key]] else NA_character_
  if (cl == "strelka") {                                 # tier-1 counts of the REF and ALT base
    ref_n <- pick(g(paste0(ref, "U")), 1); alt_n <- pick(g(paste0(alt, "U")), 1)
  } else {
    ref_n <- pick(g("AD"), 1); alt_n <- pick(g("AD"), k + 1)
  }
  # Mutect2 SB = ref_fwd, ref_rev, alt_fwd, alt_rev
  sb <- if (cl == "mutect2") c(pick(g("SB"), 3), pick(g("SB"), 4)) else c(NA_real_, NA_real_)
  list(GT = g("GT"), DP = pick(g("DP"), 1), REF_N = ref_n, ALT_N = alt_n,
       VAF = if (!is.na(ref_n + alt_n) && ref_n + alt_n > 0) alt_n / (ref_n + alt_n) else NA_real_,
       AF_reported = pick(g("AF"), k), GQ = pick(g("GQ"), 1), ALT_FWD = sb[1], ALT_REV = sb[2])
}
rec[, c("GT", "DP", "REF_N", "ALT_N", "VAF", "AF_reported", "GQ", "ALT_FWD", "ALT_REV") :=
      rbindlist(Map(sample_metrics, fmt, caller, REF, ALT, k))]
per_sample <- rec[, .(variant, setting, tumour, caller, FILTER, sample, class, GT, DP, REF_N, ALT_N,
                      VAF = round(VAF, 3), AF_reported, GQ, ALT_FWD, ALT_REV)]
per_sample[order(variant, setting, tumour, caller, class)]

# ---- 5. Per-record metrics (one row per site x file) --------------------------
record_metrics <- function(inf, cl, k, qual) {
  gi <- function(key) if (!is.null(inf[[key]])) inf[[key]] else NA_character_
  list(QUAL = suppressWarnings(as.numeric(qual)),
       # mapping quality of the ALT reads: Mutect2 MMQ (ref,alt), FreeBayes MQM; Strelka2 MQ (all reads)
       MQ_alt  = switch(cl, mutect2 = pick(gi("MMQ"), k + 1), freebayes = pick(gi("MQM"), k), strelka = pick(gi("MQ"), 1)),
       MQ_ref  = switch(cl, mutect2 = pick(gi("MMQ"), 1), freebayes = pick(gi("MQMR"), 1), strelka = NA_real_),
       # base quality of the ALT reads: Mutect2 MBQ (ref,alt); FreeBayes QA/AO (pooled); Strelka2 none
       BQ_alt  = switch(cl, mutect2 = pick(gi("MBQ"), k + 1),
                        freebayes = pick(gi("QA"), k) / pick(gi("AO"), k), strelka = NA_real_),
       TLOD = pick(gi("TLOD"), k), NLOD = pick(gi("NLOD"), k), NALOD = pick(gi("NALOD"), k),
       SomaticEVS = pick(gi("SomaticEVS"), 1), QSS = pick(gi("QSS"), 1), SNVSB = pick(gi("SNVSB"), 1),
       MQM = pick(gi("MQM"), k), AO_pooled = pick(gi("AO"), k),
       SAF = pick(gi("SAF"), k), SAR = pick(gi("SAR"), k), SRF = pick(gi("SRF"), 1), SRR = pick(gi("SRR"), 1))
}
recs <- unique(rec, by = c("variant", "setting", "tumour", "caller", "rel"))
recs[, c("QUAL_n", "MQ_alt", "MQ_ref", "BQ_alt", "TLOD", "NLOD", "NALOD", "SomaticEVS", "QSS", "SNVSB",
         "MQM", "AO_pooled", "SAF", "SAR", "SRF", "SRR") := rbindlist(Map(record_metrics, info, caller, k, QUAL))]
per_record <- recs[, .(variant, setting, tumour, caller, FILTER, QUAL = QUAL_n, MQ_alt, MQ_ref, BQ_alt, TLOD, NLOD, NALOD,
                       SomaticEVS, QSS, SNVSB, MQM, AO_pooled, SAF, SAR, SRF, SRR)]
per_record[order(variant, setting, tumour, caller)]

# Status grid for setting A: every variant x tumour x caller; absent = no record in that caller's file
status <- per_record[setting == "A", .(variant, tumour, caller, FILTER)][
  CJ(variant = SITES$variant, tumour = names(TUMOURS), caller = CALLERS), on = .(variant, tumour, caller)]
status[, status := fcase(is.na(FILTER), "absent", FILTER == "PASS", "PASS", default = "filtered")]
status[is.na(FILTER), FILTER := "absent"]
status

# ---- 6. FreeBayes: Patrik's expression, criterion by criterion ---------------
# 'INFO/MQM >= 40 && GT[normal]="0/0" && GT[tumor]="0/1" && FMT/DP[normal] >= 30 &&
#  FMT/DP[tumor] >= 30 && QUAL/INFO/AO >= 5 && SAF > 0 && SAR > 0'
fb <- rec[setting == "A" & caller == "freebayes"]
fb_crit <- fb[, {
  t <- .SD[class == "tumour"][1]; n <- .SD[class == "GM normal"][1]
  inf <- t$info[[1]]
  vals <- list(`MQM >= 40` = pick(inf$MQM, t$k), `GT normal = 0/0` = n$GT, `GT tumour = 0/1` = t$GT,
               `DP normal >= 30` = n$DP, `DP tumour >= 30` = t$DP,
               `QUAL/AO >= 5` = as.numeric(t$QUAL) / pick(inf$AO, t$k),
               `SAF > 0` = pick(inf$SAF, t$k), `SAR > 0` = pick(inf$SAR, t$k))
  ok <- c(vals[[1]] >= 40, vals[[2]] == "0/0", vals[[3]] == "0/1", vals[[4]] >= 30, vals[[5]] >= 30,
          vals[[6]] >= 5, vals[[7]] > 0, vals[[8]] > 0)
  .(criterion = names(vals), value = vapply(vals, function(v) if (is.numeric(v)) format(round(v, 2)) else as.character(v), ""),
    pass = ok, FILTER = FILTER[1])
}, by = .(variant, tumour)]
fb_crit <- fb_crit[CJ(variant = SITES$variant, tumour = names(TUMOURS), criterion = unique(fb_crit$criterion)),
                   on = .(variant, tumour, criterion)]
fb_crit[, state := fcase(is.na(pass), "unassessed", pass, "pass", default = "fail")]
fb_crit[, criterion := factor(criterion, levels = unique(fb_crit$criterion))]
fb_expr <- fb_crit[, .(expression_pass = if (any(state == "unassessed")) NA else all(state == "pass")), by = .(variant, tumour)]
fb_expr[status[caller == "freebayes"], on = .(variant, tumour), FB_FILTER := i.FILTER]
fb_expr        # expression_pass should match FB_FILTER == "PASS"

# ---- 7. Read counts in all nine samples, PER CALLER --------------------------
# Every caller's own counts, nothing picked: tumour and KI normal from the setting-A files of
# Strelka2, Mutect2, FreeBayes; WT normals (and tumours again) from the setting-B joint Mutect2 file.
# A tumour's KI normal is taken from that tumour's own files only. "absent" = that caller's file has
# no record at the site for this sample.
SOURCES <- c("strelka", "mutect2", "freebayes", "mutect2 joint (B)")
SOURCE_COLS <- c(CALLER_COLS, `mutect2 joint (B)` = "#8C564B")
pc <- per_sample[(setting == "A" & class != "WT normal") | setting == "B",
                 .(variant, sample, class, source = fifelse(setting == "B", "mutect2 joint (B)", caller),
                   FILTER, REF_N, ALT_N, DP, VAF)]
pc <- pc[order(variant, sample, source, -DP)][, .SD[1], by = .(variant, sample, source)]   # one record per source
all_samples <- unique(c(names(TUMOURS), KI_NORMALS, rec[class == "WT normal", unique(sample)]))
if (length(all_samples) != 9) message("Found ", length(all_samples), " samples, not 9 - check the WT names in the setting-B files")
grid <- CJ(variant = SITES$variant, sample = all_samples, source = SOURCES)
grid[, class := sample_class(fifelse(sample %chin% names(TUMOURS), paste0("KI_", sample), sample))]
grid <- grid[class != "WT normal" | source == "mutect2 joint (B)"]             # WT only exist in setting B
grid <- grid[!(class == "GM normal" & source == "mutect2 joint (B)")]          # KI normals are not in B
per_caller_counts <- pc[, !"class"][grid, on = .(variant, sample, source)]
per_caller_counts[, status := fcase(is.na(FILTER), "absent", FILTER == "PASS", "PASS", default = "filtered")]
per_caller_counts[, `:=`(class = factor(class, levels = CLASS_LEVELS), sample = factor(sample, levels = all_samples),
                         source = factor(source, levels = SOURCES))]
per_caller_counts[order(variant, class, sample, source), .(variant, sample, source, FILTER, REF_N, ALT_N, DP, VAF = round(VAF, 3))]

# ---- 8. Plots ------------------------------------------------------------------
# Tumour vs matched normal, per caller (setting A)
tn <- dcast(rec[setting == "A" & class %chin% c("tumour", "GM normal"),
                .(variant, tumour, caller, class = fifelse(class == "tumour", "T", "N"), VAF, DP, ALT_N)],
            variant + tumour + caller ~ class, value.var = c("VAF", "DP", "ALT_N"))
tn <- status[, .(variant, tumour, caller, status)][tn, on = .(variant, tumour, caller)]
tn[, lab := paste0(tumour, " vs ", KI_NORMALS[tumour])]
# Point labels: real sample names; ggrepel keeps them apart if installed
lab_layer <- function(...) if (requireNamespace("ggrepel", quietly = TRUE))
  ggrepel::geom_text_repel(aes(label = lab), size = 3, show.legend = FALSE, min.segment.length = 0, ...) else
  geom_text(aes(label = lab), size = 3, hjust = -0.15, show.legend = FALSE, ...)
caller_scale <- scale_colour_manual(values = CALLER_COLS)
status_shape <- scale_shape_manual(values = c(PASS = 16, filtered = 1), drop = FALSE)

# 1. Tumour VAF vs matched normal VAF
p1 <- ggplot(tn, aes(VAF_N, VAF_T, colour = caller, shape = status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, stroke = 1.2) + lab_layer() +
  facet_wrap(~ variant) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + caller_scale + status_shape +
  labs(x = "matched KI normal VAF", y = "tumour VAF", title = "1. Tumour vs matched normal VAF (each caller's own counts)") + theme_bw()
save_plot(p1, "01_vaf_tumour_vs_normal", 10, 5.5)

# 2. Tumour DP vs normal DP
p2 <- ggplot(tn, aes(DP_N, DP_T, colour = caller, shape = status)) +
  geom_hline(yintercept = 30, linetype = "dashed", colour = "grey50") + geom_vline(xintercept = 30, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, stroke = 1.2) + lab_layer() +
  facet_wrap(~ variant) + caller_scale + status_shape + expand_limits(x = 0, y = 0) +
  labs(x = "matched KI normal DP", y = "tumour DP", title = "2. Tumour vs normal depth (dashed: DP = 30, FreeBayes rule)") + theme_bw()
save_plot(p2, "02_dp_tumour_vs_normal", 10, 5)

# 3. Tumour VAF vs ALT read count
p3 <- ggplot(tn, aes(ALT_N_T, VAF_T, colour = caller, shape = status)) +
  geom_point(size = 3, stroke = 1.2) + lab_layer() +
  facet_wrap(~ variant, scales = "free_x") + caller_scale + status_shape + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "tumour ALT reads", y = "tumour VAF", title = "3. Tumour VAF vs ALT read count") + theme_bw()
save_plot(p3, "03_vaf_vs_altreads_tumour", 10, 5)

# 4. REF and ALT reads, all nine samples, one bar per caller (absent = no record)
facet_ns <- if (requireNamespace("ggh4x", quietly = TRUE)) ggh4x::facet_nested(variant ~ class + sample, scales = "free_x", space = "free_x") else
  facet_grid(variant ~ class + sample, scales = "free_x", space = "free_x")
pcl <- melt(per_caller_counts, id.vars = c("variant", "sample", "class", "source", "status", "VAF"),
            measure.vars = c("REF_N", "ALT_N"), variable.name = "allele", value.name = "reads")
pcl[, allele := factor(sub("_N", "", allele), levels = c("REF", "ALT"))]
lab4 <- per_caller_counts[, .(variant, sample, class, source, top = fcoalesce(REF_N + ALT_N, 0),
                              label = fifelse(status == "absent", "absent", sprintf("%.2f", VAF)))]
p4 <- ggplot(pcl[!is.na(reads)], aes(source, reads, fill = allele)) +
  geom_col(width = 0.75) +
  geom_text(data = lab4, aes(source, top, label = label), inherit.aes = FALSE, vjust = -0.4, size = 2.5) +
  facet_ns + scale_fill_manual(values = c(REF = "grey70", ALT = "#C62828")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "reads (each caller's own count)", fill = NULL,
       title = "4. REF and ALT reads per caller in all nine samples", subtitle = "number above bar = VAF; absent = no record in that caller's file") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), strip.text = element_text(size = 7))
save_plot(p4, "04_ref_alt_nine_samples_per_caller", 16, 7)

# 5. VAF per caller across the nine samples
p5 <- ggplot(per_caller_counts[!is.na(VAF)], aes(sample, VAF, colour = source, shape = status)) +
  geom_point(size = 3, stroke = 1.1, position = position_dodge(width = 0.6)) +
  facet_grid(variant ~ class, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = SOURCE_COLS) + scale_shape_manual(values = c(PASS = 16, filtered = 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "VAF (each caller's own counts)", colour = "caller", shape = "FILTER",
       title = "5. VAF per caller in tumours, GM (KI) normals and WT normals",
       subtitle = "Missing point = that caller has no record at the site for the sample (no data, not VAF 0)") + theme_bw()
save_plot(p5, "05_vaf_nine_samples_per_caller", 11, 6)

# 6. Caller FILTER / QC status per tumour (+ FreeBayes expression)
st6 <- rbind(status[, .(variant, tumour, row = caller, status, label = FILTER)],
             fb_expr[, .(variant, tumour, row = "freebayes expression",
                         status = fcase(is.na(expression_pass), "unassessed", expression_pass, "PASS", default = "filtered"),
                         label = fcase(is.na(expression_pass), "n/a", expression_pass, "pass", default = "fail"))])
st6[, row := factor(row, levels = rev(c(CALLERS, "freebayes expression")))]
p6 <- ggplot(st6, aes(tumour, row, fill = status)) +
  geom_tile(colour = "grey40") + geom_text(aes(label = label), size = 3) +
  facet_wrap(~ variant) + scale_fill_manual(values = STATUS_COLS, drop = FALSE) +
  labs(x = NULL, y = NULL, title = "6. FILTER in each caller's own file, per tumour") + theme_minimal()
save_plot(p6, "06_filter_status_heatmap", 10, 4)

# 7. FreeBayes genotype vs tumour VAF (FreeBayes' own counts)
fb7 <- rec[setting == "A" & caller == "freebayes" & class == "tumour", .(variant, tumour, GT, VAF, GQ)]
fb7 <- status[caller == "freebayes", .(variant, tumour, status)][fb7, on = .(variant, tumour)]
p7 <- ggplot(fb7, aes(GT, VAF, shape = status)) +
  geom_hline(yintercept = c(0.5, 1), linetype = "dotted", colour = "grey50") +
  geom_point(size = 3.5, stroke = 1.2, colour = CALLER_COLS[["freebayes"]]) +
  geom_text(aes(label = paste0(tumour, " (GQ ", GQ, ")")), nudge_x = 0.3, size = 3) +
  facet_wrap(~ variant) + scale_y_continuous(limits = c(0, 1)) + status_shape +
  labs(x = "FreeBayes tumour genotype", y = "tumour VAF (FreeBayes AD)",
       title = "7. FreeBayes genotype vs observed tumour VAF", subtitle = "Patrik's expression keeps tumour GT 0/1 only") + theme_bw()
save_plot(p7, "07_freebayes_gt_vs_vaf", 10, 5)

# 8. FreeBayes expression: each criterion
p8 <- ggplot(fb_crit, aes(criterion, tumour, fill = state)) +
  geom_tile(colour = "grey40") + geom_text(aes(label = fifelse(is.na(value), "", value)), size = 2.8) +
  facet_wrap(~ variant, ncol = 1) +
  scale_fill_manual(values = c(pass = "#A5D6A7", fail = "#EF9A9A", unassessed = "#E0E0E0")) +
  labs(x = NULL, y = NULL, fill = NULL, title = "8. Patrik's FreeBayes expression, criterion by criterion",
       subtitle = "MQM, QUAL/AO, SAF, SAR are pooled over tumour + normal (FreeBayes INFO)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p8, "08_freebayes_criteria_heatmap", 10, 6)

# ---- 9. Outputs --------------------------------------------------------------
fwrite(per_sample, file.path(OUT, "per_sample_metrics.tsv"), sep = "\t")
fwrite(per_record, file.path(OUT, "per_record_metrics.tsv"), sep = "\t")
fwrite(status, file.path(OUT, "filter_status.tsv"), sep = "\t")
fwrite(fb_crit, file.path(OUT, "freebayes_criteria.tsv"), sep = "\t")
fwrite(fb_expr, file.path(OUT, "freebayes_expression.tsv"), sep = "\t")
fwrite(per_caller_counts, file.path(OUT, "nine_sample_counts_per_caller.tsv"), sep = "\t")
saveRDS(list(sites = SITES, files = files, rec = rec, per_sample = per_sample, per_record = per_record,
             status = status, fb_crit = fb_crit, fb_expr = fb_expr, per_caller_counts = per_caller_counts),
        file.path(OUT, "example_sites.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_example_sites"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ", ", bcf("--version")[1], ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Sites:      ", paste(SITES$variant, collapse = "; ")),
  paste0("- Inputs:     ", nrow(files), " single-caller files under ", ROOT, " (settings A and B; merged/union files not used)"),
  "- Nine-sample counts: VCF FORMAT (no BAMs)",
  paste0("- Outputs:    ", OUT, " (plots 01-08 PDF+PNG, TSVs, example_sites.rds)"),
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "secs")), " s. Outputs: ", OUT)
