#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24
#
# mouse_wgs_sample_identity.R - are the tumour / matched KI normal pairs from the same mouse?
#
# Why: at chr10:3184122 G>T, HCC_1 is 1/1 while its "matched" KI1 is 0/0, and HCC_2 is 0/0 while KI2
# is 1/1 - swapped normals would make both pairs consistent. A germline genotype must agree between a
# tumour and the normal liver of the same mouse (except rare LOH / somatic sites), so genome-wide
# genotype concordance tells which normal belongs to which tumour.
#
# Data: the FreeBayes files of setting A (tumor_vs_normal/snv_indel/, one file per tumour/normal pair).
# FreeBayes is a germline-style caller: it genotypes BOTH samples at every variable site, germline
# sites included, whatever its FILTER. All records are used (not only PASS).
#
# Blocks:
#   0  EDIT THIS BLOCK
#   1  setup
#   2  FreeBayes files
#   3  read SNV genotypes of both samples per file (QUAL >= MIN_QUAL, autosomes)
#   4  one table: site x sample (dosage 0/1/2, DP, VAF); informative sites
#   5  pairwise concordance for all sample pairs
#   6  verdict: best-matching normal per tumour vs the labelled one
#   7  plots: concordance heatmap, VAF scatter tumour x normal
#   8  outputs, RUN_LOG.md
#
# Definitions
#   dosage      0/0 = 0, 0/1 = 1, 1/1 = 2 (FreeBayes GT; phased "|" treated as "/")
#   usable      the sample's genotype is called and its DP >= MIN_DP
#   informative a site where >= 2 samples are usable and they do not all have the same dosage
#               (sites identical in all mice, e.g. fixed differences to the reference, say nothing
#               about identity in an inbred colony)
#   concordance fraction of informative sites, usable in both samples, with the same dosage
#   opposite_hom sites where one sample is 0/0 and the other 1/1 (strong evidence of different animals)
# Same mouse: concordance close to 1 (> 0.95 expected; tumour LOH and somatic calls lower it a little).
#
# Limitations
# - A site is in a pair's file only if FreeBayes found variation in that pair; when a file has no
#   record the two samples are "unknown" (not assumed 0/0). Cross-file pairs are therefore compared
#   only at sites present in both files.
# - Inbred C57BL/6J mice share most of the genome: there may be few informative sites. Their number
#   is reported for every pair; a verdict from a handful of sites is weak.
# - WT normals are not included: they exist only in the Mutect2 joint files (setting B), and Mutect2
#   does not genotype germline sites in the same way.
#
# Env on ssi: r-mutation (data.table, ggplot2, bcftools).
# Output: results/mouse_wgs/mouse_wgs_sample_identity/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
DIR_A    <- file.path(REPO, "data/mouse/wgs/tumor_vs_normal/snv_indel")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_sample_identity")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_sample_identity.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")

TUMOURS    <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1", HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
                HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing

MIN_QUAL <- 30      # FreeBayes site QUAL
MIN_DP   <- 20      # per-sample depth for a genotype to count

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
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  invisible(p)
}
short_name <- function(s) fcase(grepl("KI_HCC_[0-9]", s), sub(".*KI_(HCC_[0-9]).*", "\\1", s),
                                default = sub("^(id_[0-9]+_)?(sfb_dikic_)?([0-9]{4}_[0-9]{2}_[0-9]{2}_)?", "", s))

# ---- 2. FreeBayes files --------------------------------------------------------
fb_files <- data.table(path = list.files(DIR_A, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf|vcf\\.gz)$"))
fb_files <- fb_files[grepl("freebayes", path, ignore.case = TRUE) & !grepl("union|merged", path, ignore.case = TRUE)]
fb_files[, tumour := NA_character_]
for (k in names(TUMOURS)) fb_files[grepl(TUMOURS[[k]], basename(path), fixed = TRUE), tumour := k]
fb_files[, file := basename(path)]
fb_files[, .(tumour, file)]
stopifnot(nrow(fb_files) == length(TUMOURS), !anyNA(fb_files$tumour))

# ---- 3. Read genotypes ---------------------------------------------------------
# SNVs only, QUAL >= MIN_QUAL, all FILTER values. One row per site per sample.
read_fb <- function(path) {
  smp <- system2(BCFTOOLS, c("query", "-l", shQuote(path)), stdout = TRUE)
  cmd <- paste(shQuote(BCFTOOLS), "view -v snps -i", shQuote(paste0("QUAL>=", MIN_QUAL)), shQuote(path), "|",
               shQuote(BCFTOOLS), "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%DP\\t%AD]\\n"))
  cn <- c("CHROM", "POS", "REF", "ALT", paste0(rep(c("GT_", "DP_", "AD_"), length(smp)), rep(seq_along(smp), each = 3)))
  x <- fread(cmd = cmd, sep = "\t", header = FALSE, col.names = cn, colClasses = "character")
  x <- x[sub("^chr", "", CHROM) %chin% as.character(1:19)]          # autosomes
  x[, site := paste(CHROM, POS, REF, ALT, sep = ":")]
  rbindlist(lapply(seq_along(smp), function(j) {
    ad <- tstrsplit(x[[paste0("AD_", j)]], ",", fixed = TRUE)
    data.table(site = x$site, sample = short_name(smp[j]), GT = x[[paste0("GT_", j)]],
               DP = suppressWarnings(as.numeric(x[[paste0("DP_", j)]])),
               REF_N = suppressWarnings(as.numeric(ad[[1]])),
               ALT_N = suppressWarnings(as.numeric(if (length(ad) > 1) ad[[2]] else NA)))
  }))
}
gt <- rbindlist(lapply(fb_files$path, read_fb))
gt[, .(sites = .N, called = sum(!grepl("\\.", GT)), dp_ok = sum(DP >= MIN_DP, na.rm = TRUE)), by = sample]

# Dosage and VAF; a genotype counts only with DP >= MIN_DP
gt[, GT := sub("|", "/", GT, fixed = TRUE)]
gt[, dose := fcase(GT == "0/0", 0L, GT %chin% c("0/1", "1/0"), 1L, GT == "1/1", 2L, default = NA_integer_)]
gt[is.na(DP) | DP < MIN_DP, dose := NA_integer_]
gt[, VAF := fifelse(REF_N + ALT_N > 0, ALT_N / (REF_N + ALT_N), NA_real_)]
gt <- unique(gt, by = c("site", "sample"))
gt[, .N, by = .(sample, dose)][order(sample, dose)]

# ---- 4. Site x sample table; informative sites ----------------------------------
SAMPLES <- c(names(TUMOURS), unname(KI_NORMALS))
SAMPLES <- intersect(SAMPLES, unique(gt$sample))
dose_w <- dcast(gt, site ~ sample, value.var = "dose")
vaf_w  <- dcast(gt[!is.na(dose)], site ~ sample, value.var = "VAF")
M <- as.matrix(dose_w[, ..SAMPLES])
n_usable <- rowSums(!is.na(M))
rmin <- suppressWarnings(do.call(pmin, c(as.data.frame(M), na.rm = TRUE)))
rmax <- suppressWarnings(do.call(pmax, c(as.data.frame(M), na.rm = TRUE)))
dose_w[, informative := n_usable >= 2 & rmax > rmin]
c(sites = nrow(dose_w), informative = sum(dose_w$informative))

# ---- 5. Pairwise concordance ---------------------------------------------------
pairs <- CJ(s1 = SAMPLES, s2 = SAMPLES)[s1 != s2]
info_sites <- dose_w[informative == TRUE]
pairwise <- pairs[, {
  a <- info_sites[[s1]]; b <- info_sites[[s2]]
  ok <- !is.na(a) & !is.na(b)
  va <- vaf_w[info_sites[ok], on = "site"][[s1]]; vb <- vaf_w[info_sites[ok], on = "site"][[s2]]
  .(n_sites = sum(ok), concordance = round(mean(a[ok] == b[ok]), 3),
    opposite_hom = sum(abs(a[ok] - b[ok]) == 2),
    vaf_cor = if (sum(ok) > 2) round(cor(va, vb, use = "complete.obs"), 3) else NA_real_)
}, by = .(s1, s2)]
pairwise[, labelled_pair := (s1 %chin% names(KI_NORMALS) & KI_NORMALS[s1] == s2) |
                            (s2 %chin% names(KI_NORMALS) & KI_NORMALS[s2] == s1)]
pairwise[, labelled_pair := fcoalesce(labelled_pair, FALSE)]
pairwise[s1 < s2][order(-concordance)]

# ---- 6. Verdict per tumour -----------------------------------------------------
tn <- pairwise[s1 %chin% names(TUMOURS) & s2 %chin% KI_NORMALS]
verdict <- tn[order(s1, -concordance)][, .(labelled_normal = KI_NORMALS[s1[1]],
                                           best_normal = s2[1], best_concordance = concordance[1], best_n = n_sites[1],
                                           labelled_concordance = concordance[s2 == KI_NORMALS[s1[1]]][1],
                                           second_concordance = concordance[2]), by = .(tumour = s1)]
verdict[, consistent := best_normal == labelled_normal]
verdict          # consistent = FALSE -> the labelled normal is not the best genotype match
fwrite(verdict, file.path(OUT, "verdict.tsv"), sep = "\t")

# ---- 7. Plots ------------------------------------------------------------------
# Concordance heatmap (informative sites); labelled tumour/normal pairs outlined
hm <- copy(pairwise)
hm[, `:=`(s1 = factor(s1, levels = SAMPLES), s2 = factor(s2, levels = rev(SAMPLES)))]
p_hm <- ggplot(hm, aes(s1, s2, fill = concordance)) +
  geom_tile(colour = "white") +
  geom_tile(data = hm[labelled_pair == TRUE], colour = "black", linewidth = 1, fill = NA) +
  geom_text(aes(label = sprintf("%.2f\nn=%s", concordance, format(n_sites, big.mark = ","))), size = 3) +
  scale_fill_gradient(low = "white", high = "#2E7D32", limits = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "concordance",
       title = "Genotype concordance at informative germline SNVs (FreeBayes, all records)",
       subtitle = sprintf("QUAL >= %d, DP >= %d in both samples, autosomes; black box = labelled tumour/normal pair",
                          MIN_QUAL, MIN_DP)) +
  coord_equal() + theme_minimal()
save_plot(p_hm, "concordance_heatmap", 8, 7)

# VAF scatter: each tumour vs each KI normal at informative sites
sc <- rbindlist(lapply(names(TUMOURS), function(t) rbindlist(lapply(unname(KI_NORMALS), function(n) {
  d <- vaf_w[info_sites[!is.na(get(t)) & !is.na(get(n)), .(site)], on = "site"]
  if (!nrow(d)) return(NULL)
  data.table(tumour = t, normal = n, vaf_t = d[[t]], vaf_n = d[[n]], labelled = KI_NORMALS[[t]] == n)
}))))
p_sc <- ggplot(sc, aes(vaf_n, vaf_t)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = labelled), alpha = 0.3, size = 0.8) +
  facet_grid(tumour ~ normal) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c(`TRUE` = "#C62828", `FALSE` = "grey30"), labels = c(`TRUE` = "labelled pair", `FALSE` = "other")) +
  labs(x = "KI normal VAF", y = "tumour VAF", colour = NULL,
       title = "VAF at informative germline SNVs: same mouse = points on the diagonal / at 0, 0.5, 1") + theme_bw()
save_plot(p_sc, "vaf_scatter_tumour_vs_normal", 9, 9)

# ---- 8. Outputs --------------------------------------------------------------
fwrite(pairwise, file.path(OUT, "pairwise_concordance.tsv"), sep = "\t")
fwrite(dose_w[, .N, by = informative], file.path(OUT, "sites_summary.tsv"), sep = "\t")
saveRDS(list(files = fb_files, gt = gt, dose_w = dose_w, pairwise = pairwise, verdict = verdict,
             params = list(MIN_QUAL = MIN_QUAL, MIN_DP = MIN_DP)), file.path(OUT, "sample_identity.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_sample_identity"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Parameters: MIN_QUAL=", MIN_QUAL, ", MIN_DP=", MIN_DP, ", SNVs, autosomes, all FILTER values"),
  paste0("- Inputs:     ", paste(fb_files$file, collapse = ", ")),
  paste0("- Counts/QC:  sites ", nrow(dose_w), ", informative ", sum(dose_w$informative)),
  paste0("- Verdict:    ", paste(verdict[, sprintf("%s best=%s (%.2f, labelled %s %.2f)", tumour, best_normal,
                                                   best_concordance, labelled_normal, labelled_concordance)], collapse = "; ")),
  paste0("- Outputs:    ", OUT),
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
