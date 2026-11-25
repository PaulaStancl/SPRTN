#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_identity_informativeness.R - how informative are the identity sites of mouse_wgs_somalier.R,
# how much did the tumours influence which sites were kept, and a short list of sites for a BAM check.
# Reads the outputs of the finished mouse_wgs_somalier.R run; does not rerun it (fast, no bcftools calling).
#
# Why: the sites come from the per-pair FreeBayes files (variant-only output). A site is in a pair's file
# only if the tumour OR the normal of that pair is non-ref, and a site had to be in all 3 files. So a
# tumour can decide whether a site is kept: where its labelled normal is 0/0, the site is in that pair's
# file only because the tumour is non-ref ("tumour-dependent retention"). The site filters themselves use
# only the normals; the site list still depends on the tumours through the files.
#
# Questions:
#   1. how many sites actually distinguish the three KI mice (genotype patterns across the KI normals)
#   2. per pair of normals: sites with different genotypes and opposite homozygotes (IBS0)
#   3. per normal: sites where it differs from both other normals (private)
#   4. per tumour: sites kept only because of the tumour (labelled normal 0/0)
#   5. the verdict on a TUMOUR-INDEPENDENT site set: all three KI normals non-ref, so every file holds the
#      site because of its normal - no site depends on a tumour. Metric: discordance at sites where the
#      normal is 1/1 (tumour not 1/1), all 9 tumour x normal combinations, exact 95% CIs.
#   6. chromosomal distribution of the informative sites
#   7. bam_check_sites_for_patrik.csv: opposite-homozygous sites between each pair of KI normals, for a
#      read-level check in the BAMs (IGV or bcftools mpileup). Tumours are NOT used to choose these sites
#      (only: tumour genotype not masked), so the list does not presuppose the swap.
#
# Site sets: "A" = all sites of panel A (masked_6samples.bcf, unthinned); "A_20kb" = the primary sites
# (sites_used_primary.tsv, thinned to one per 20 kb).
#
# Env on ssi: r-mutation (data.table, ggplot2, bcftools).
# Output: results/mouse_wgs/mouse_wgs_identity_informativeness/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
SOM      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_somalier")      # outputs of mouse_wgs_somalier.R
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_identity_informativeness")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_identity_informativeness.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")

TUMOURS    <- c("HCC_1", "HCC_2", "HCC_3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing (= which file)
HYPOTHESIS <- c(HCC_1 = "KI2", HCC_2 = "KI1", HCC_3 = "KI3")   # genetic match from mouse_wgs_somalier.R

N_PER_TEST <- 5     # BAM-check sites per pair of normals (at most one per chromosome)
EXAMPLE_SITES <- data.frame(CHROM = c("10", "1"), POS = c(3184122L, 5609923L))   # added to the CSV if in panel A

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
NORMALS <- unname(KI_NORMALS)
SAMPLES <- c(TUMOURS, NORMALS)
GT <- c("0/0", "0/1", "1/1")

MASKED <- file.path(SOM, "masked", "masked_6samples.bcf")
if (!file.exists(MASKED)) MASKED <- file.path(SOM, "masked", "masked_6samples.vcf")
PRIMARY_TSV <- file.path(SOM, "sites_used_primary.tsv")
for (f in c(MASKED, PRIMARY_TSV)) if (!file.exists(f)) stop("missing input from mouse_wgs_somalier.R: ", f)
if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS)

# exact (Clopper-Pearson) 95% CI
ci_lo <- function(x, n) fifelse(n > 0, qbeta(0.025, x, n - x + 1), NA_real_)
ci_hi <- function(x, n) fifelse(n > 0, qbeta(0.975, x + 1, n - x), NA_real_)
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150, limitsize = FALSE)
}

# ---- 2. Read the masked 6-sample genotypes -------------------------------------
smp <- system2(BCFTOOLS, c("query", "-l", shQuote(MASKED)), stdout = TRUE)
stopifnot(setequal(smp, SAMPLES))
fmt <- "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%AD]\\n"
cols <- c("CHROM", "POS", "REF", "ALT", as.vector(rbind(paste0("GT_", smp), paste0("AD_", smp))))
A <- fread(cmd = paste(shQuote(BCFTOOLS), "query -f", shQuote(fmt), shQuote(MASKED)), sep = "\t", header = FALSE,
           quote = "", col.names = cols, colClasses = c("character", "integer", rep("character", length(cols) - 2)))
for (s in SAMPLES) {
  A[, (s) := match(get(paste0("GT_", s)), GT) - 1L]                             # 0/1/2, NA = masked
  ad <- tstrsplit(A[[paste0("AD_", s)]], ",", fixed = TRUE)
  A[, paste0("ref_", s) := suppressWarnings(as.integer(ad[[1]]))]
  A[, paste0("alt_", s) := if (length(ad) > 1) suppressWarnings(as.integer(ad[[2]])) else NA_integer_]
}
A[, key := paste(CHROM, POS, REF, ALT, sep = ":")]
stopifnot(!anyNA(A[, ..NORMALS]))                             # normals are never masked by construction
prim <- fread(PRIMARY_TSV, colClasses = list(character = c("CHROM", "REF", "ALT")))
prim_keys <- prim[, paste(CHROM, POS, REF, ALT, sep = ":")]
stopifnot(all(prim_keys %chin% A$key))
SETS <- list(A = A$key, A_20kb = prim_keys)
message("panel A: ", nrow(A), " sites; A_20kb: ", length(prim_keys), " sites")
masked_per_sample <- A[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = SAMPLES]
print(masked_per_sample)

# ---- 3. Informativeness of the KI normals ------------------------------------
A[, pattern := paste(GT[KI1 + 1], GT[KI2 + 1], GT[KI3 + 1], sep = " | ")]
A[, n_distinct := apply(.SD, 1, uniqueN), .SDcols = NORMALS]
A[, opp_hom := do.call(pmin, .SD) == 0L & do.call(pmax, .SD) == 2L, .SDcols = NORMALS]
A[, class := fcase(n_distinct == 1, "identical in all 3 normals (not informative)",
                   opp_hom,         "differ, with opposite homozygotes (IBS0)",
                   default =        "differ, het vs hom only")]
A[, all_normals_nonref := KI1 > 0 & KI2 > 0 & KI3 > 0]

by_set <- function(fun) rbindlist(lapply(names(SETS), function(st) fun(A[key %chin% SETS[[st]]])[, set := st]), fill = TRUE)

# 3a. genotype patterns across the 3 normals
patterns <- by_set(function(d) d[, .N, by = .(pattern, class)][, pct := round(100 * N / sum(N), 1)][order(-N)])
setcolorder(patterns, "set")
print(patterns)

# 3b. summary
summary_inf <- by_set(function(d) d[, .(sites = .N,
  identical_normals = sum(n_distinct == 1), pct_identical = round(100 * mean(n_distinct == 1), 1),
  informative = sum(n_distinct > 1), pct_informative = round(100 * mean(n_distinct > 1), 1),
  with_opposite_homs = sum(opp_hom), all_normals_nonref = sum(all_normals_nonref))])
setcolorder(summary_inf, "set")
print(summary_inf)

# 3c. per pair of normals: different genotypes, opposite homozygotes
np <- as.data.table(t(combn(NORMALS, 2))); setnames(np, c("n1", "n2"))
normal_pairs <- by_set(function(d) np[, {
  a <- d[[n1]]; b <- d[[n2]]
  .(sites = length(a), different = sum(a != b), pct_different = round(100 * mean(a != b), 1),
    opposite_homs = sum(abs(a - b) == 2L), pct_opposite_homs = round(100 * mean(abs(a - b) == 2L), 2))
}, by = .(n1, n2)])
setcolorder(normal_pairs, "set")
print(normal_pairs)

# 3d. per normal: differs from both other normals (private genotype)
private <- by_set(function(d) rbindlist(lapply(NORMALS, function(n) {
  o <- setdiff(NORMALS, n)
  data.table(normal = n, private = sum(d[[n]] != d[[o[1]]] & d[[n]] != d[[o[2]]]),
             private_opposite_hom = sum(abs(d[[n]] - d[[o[1]]]) == 2L & abs(d[[n]] - d[[o[2]]]) == 2L))
})))
setcolorder(private, "set")
print(private)

# ---- 4. Tumour-dependent retention ------------------------------------------
# Labelled normal 0/0: the site is in that pair's file only because the tumour is non-ref.
tumour_dep <- by_set(function(d) rbindlist(lapply(TUMOURS, function(t) {
  n <- KI_NORMALS[[t]]; x <- d[get(n) == 0L]
  data.table(tumour = t, file_normal = n, sites = nrow(d), normal_0_0 = nrow(x),
             pct = round(100 * nrow(x) / nrow(d), 1),
             tumour_0_1 = sum(x[[t]] == 1L, na.rm = TRUE), tumour_1_1 = sum(x[[t]] == 2L, na.rm = TRUE),
             tumour_0_0 = sum(x[[t]] == 0L, na.rm = TRUE), tumour_masked = sum(is.na(x[[t]])))
})))
setcolorder(tumour_dep, "set")
print(tumour_dep)
A[, tumour_dependent := KI1 == 0L | KI2 == 0L | KI3 == 0L]   # some file holds the site only because of its tumour

# ---- 5. Verdict on the tumour-independent site set ---------------------------
# All three KI normals non-ref: each file holds the site because of its own normal. Discordance at sites
# where the normal is 1/1: tumour genotype (not masked) != 1/1. LOH cannot turn 1/1 into 0/1 or 0/0.
indep <- by_set(function(d) {
  d <- d[all_normals_nonref == TRUE]
  CJ(tumour = TUMOURS, normal = NORMALS)[, {
    tg <- d[[tumour]]; ng <- d[[normal]]; inf <- ng == 2L & !is.na(tg)
    .(sites = nrow(d), informative = sum(inf), discordant = sum(inf & tg != 2L))
  }, by = .(tumour, normal)]
})
indep[, `:=`(rate = round(discordant / informative, 4), ci_lo = round(ci_lo(discordant, informative), 4),
             ci_hi = round(ci_hi(discordant, informative), 4),
             pairing = fcase(KI_NORMALS[tumour] == normal & HYPOTHESIS[tumour] == normal, "labelled = genetic match",
                             KI_NORMALS[tumour] == normal, "labelled",
                             HYPOTHESIS[tumour] == normal, "genetic match (swap)", default = "other"))]
setcolorder(indep, "set")
print(indep[order(set, tumour, rate)])

p_indep <- ggplot(indep, aes(normal, rate, colour = pairing)) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.5) +
  geom_text(aes(label = paste0(discordant, "/", informative)), vjust = -1.2, size = 3, show.legend = FALSE) +
  facet_grid(set ~ tumour) +
  scale_colour_manual(values = c("labelled" = "#C62828", "genetic match (swap)" = "#2E7D32",
                                 "labelled = genetic match" = "#1565C0", "other" = "grey40"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = "KI normal", y = "discordance at normal 1/1 sites (95% CI)",
       title = "Identity on tumour-independent sites (all 3 KI normals non-ref)",
       subtitle = "Every file holds these sites because of its normal; tumours play no part in site selection") +
  theme_bw() + theme(legend.position = "bottom")
save_plot(p_indep, "tumour_independent_discordance", 10, 6)

# ---- 6. Chromosomal distribution of informative sites ------------------------
chr_inf <- by_set(function(d) d[, .(sites = .N, informative = sum(n_distinct > 1), opposite_homs = sum(opp_hom),
                                     tumour_dependent = sum(tumour_dependent)), by = CHROM])
chr_inf[, chr_n := as.integer(sub("^chr", "", CHROM))]
setorder(chr_inf, set, chr_n)
setcolorder(chr_inf, "set")
print(chr_inf)
cl <- melt(chr_inf, id.vars = c("set", "CHROM", "chr_n"), measure.vars = c("sites", "informative", "opposite_homs"))
p_chr <- ggplot(cl, aes(factor(chr_n), value, fill = variable)) + geom_col(position = "dodge") +
  facet_wrap(~ set, ncol = 1, scales = "free_y") +
  labs(x = "chromosome", y = "sites", fill = NULL, title = "Identity sites per chromosome: all, informative (normals differ), opposite homozygotes") +
  theme_bw() + theme(legend.position = "bottom")
save_plot(p_chr, "informative_sites_per_chromosome", 11, 6)

# ---- 7. Sites for a BAM check (Patrik) --------------------------------------
# Opposite homozygotes between each pair of KI normals; all tumour genotypes not masked (tumours not
# otherwise used, so the list does not presuppose the swap). Deepest sites first (minimum REF+ALT over
# the six samples), at most one per chromosome per test, N_PER_TEST per test, no site twice.
dp_cols <- paste0("dp_", SAMPLES)
for (s in SAMPLES) A[, paste0("dp_", s) := get(paste0("ref_", s)) + get(paste0("alt_", s))]
A[, min_dp := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = dp_cols]
cand <- A[!is.na(HCC_1) & !is.na(HCC_2) & !is.na(HCC_3)]
used <- character(); bam <- list()
for (i in seq_len(nrow(np))) {
  a <- np$n1[i]; b <- np$n2[i]
  x <- cand[abs(get(a) - get(b)) == 2L & !key %chin% used][order(-min_dp)]
  x <- x[!duplicated(CHROM)][seq_len(min(.N, N_PER_TEST))]
  used <- c(used, x$key)
  bam[[i]] <- x[, test := paste0(a, " vs ", b, " (opposite homozygotes)")]
}
ex <- A[as.data.table(EXAMPLE_SITES), on = .(CHROM, POS), nomatch = NULL][!key %chin% used]
if (nrow(ex)) bam[[length(bam) + 1]] <- ex[, test := "example site (earlier analysis)"]
bam <- rbindlist(bam, fill = TRUE)
bam_out <- bam[, .(test, locus = paste0(CHROM, ":", POS), CHROM, POS, REF, ALT)]
for (s in SAMPLES) set(bam_out, j = paste0("GT_", s), value = fifelse(is.na(bam[[s]]), "masked", GT[bam[[s]] + 1L]))
for (s in SAMPLES) set(bam_out, j = paste0("reads_ref_alt_", s), value = paste0(bam[[paste0("ref_", s)]], ",", bam[[paste0("alt_", s)]]))
bam_out[, `:=`(BAM_reads_ref_alt_HCC_1 = "", BAM_reads_ref_alt_HCC_2 = "", BAM_reads_ref_alt_HCC_3 = "",
               BAM_reads_ref_alt_KI1 = "", BAM_reads_ref_alt_KI2 = "", BAM_reads_ref_alt_KI3 = "", comment = "")]
print(bam_out[, 1:12])
fwrite(bam_out, file.path(OUT, "bam_check_sites_for_patrik.csv"))

# ---- 8. Outputs, RUN_LOG.md --------------------------------------------------
fwrite(patterns, file.path(OUT, "normal_genotype_patterns.tsv"), sep = "\t")
fwrite(summary_inf, file.path(OUT, "informativeness_summary.tsv"), sep = "\t")
fwrite(normal_pairs, file.path(OUT, "normal_pairs.tsv"), sep = "\t")
fwrite(private, file.path(OUT, "normal_private_sites.tsv"), sep = "\t")
fwrite(tumour_dep, file.path(OUT, "tumour_dependent_retention.tsv"), sep = "\t")
fwrite(indep, file.path(OUT, "tumour_independent_discordance.tsv"), sep = "\t")
fwrite(chr_inf, file.path(OUT, "informative_sites_per_chromosome.tsv"), sep = "\t")
fwrite(masked_per_sample, file.path(OUT, "masked_per_sample.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
in_md5 <- tools::md5sum(c(MASKED, PRIMARY_TSV))
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_identity_informativeness"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Command:    Rscript ", SCRIPT),
  paste0("- Inputs:     ", MASKED, " (md5 ", in_md5[[1]], "); ", PRIMARY_TSV, " (md5 ", in_md5[[2]], ")"),
  paste0("- Parameters: labelled ", paste(names(KI_NORMALS), KI_NORMALS, sep = "-", collapse = ", "),
         "; genetic match ", paste(names(HYPOTHESIS), HYPOTHESIS, sep = "-", collapse = ", "), "; N_PER_TEST ", N_PER_TEST),
  paste0("- Tools:      ", system2(BCFTOOLS, "--version", stdout = TRUE)[1], "; data.table ", packageVersion("data.table")),
  paste0("- Outputs:    ", OUT),
  paste0("- Counts:     panel A ", nrow(A), " sites, A_20kb ", length(prim_keys), "; tumour-independent (all normals non-ref) A ",
         A[all_normals_nonref == TRUE, .N], "; BAM-check sites ", nrow(bam_out)),
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
