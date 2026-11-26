#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-26
#
# mouse_wgs_identity_discordant_sites.R - the germline SNVs where a tumour's genotype differs from a KI
# normal's (the red triangles of the slide "Allele fractions at all germline SNPs", ~520 per labelled pair),
# as tables with genotypes, read counts and VAFs of all six samples. For choosing sites to check in the BAMs.
# Reads the outputs of mouse_wgs_somalier.R; no re-analysis.
#
# Sites: sites_used_primary.tsv (panel A, primary thinning). A site counts for a tumour x normal pair when
# both genotypes are known (not masked); it differs when the genotypes differ. All normal genotypes are
# used (0/0, 0/1, 1/1), as on the slide.
#
# Env on ssi: r-mutation (data.table, bcftools). Run:
#   Rscript scripts/mouse_wgs/mouse_wgs_identity_discordant_sites.R
# Output: results/mouse_wgs/mouse_wgs_identity_discordant_sites/
#   counts_3x3.tsv                     differing / compared sites per tumour x normal (compare with the slide)
#   discordant_sites_all_pairs.tsv     every differing site, all 9 tumour x normal pairs
#   discordant_<tumour>_vs_<normal>.tsv one file per labelled pair (HCC_1 vs KI1, HCC_2 vs KI2, HCC_3 vs KI3)

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
IN       <- file.path(REPO, "results/mouse_wgs/mouse_wgs_somalier")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_identity_discordant_sites")
PATRIK   <- file.path(REPO, "results/mouse_wgs/mouse_wgs_identity_informativeness/bam_check_sites_for_patrik.csv")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_identity_discordant_sites.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
# Console as wide as the terminal (falls back to 80 when there is no terminal, e.g. batch Rscript)
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
started <- Sys.time()
TUMOURS <- names(KI_NORMALS); NORMALS <- unname(KI_NORMALS); SAMPLES <- c(TUMOURS, NORMALS)
GT <- c("0/0", "0/1", "1/1")

# ---- 2. Genotypes (primary sites) and read counts (masked BCF) ----------------
su <- fread(file.path(IN, "sites_used_primary.tsv"), colClasses = list(character = c("CHROM", "REF", "ALT")))
su[, POS := as.integer(POS)]
mf <- file.path(IN, "masked", "masked_6samples.bcf")
if (!file.exists(mf)) mf <- file.path(IN, "masked", "masked_6samples.vcf")
smp <- system2(BCFTOOLS, c("query", "-l", shQuote(mf)), stdout = TRUE)
ad <- fread(cmd = paste(shQuote(BCFTOOLS), "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n"), shQuote(mf)),
            header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", paste0("AD_", smp)), colClasses = "character")
ad[, POS := as.integer(POS)]
for (s in SAMPLES) {
  p <- tstrsplit(ad[[paste0("AD_", s)]], ",", fixed = TRUE)
  r <- suppressWarnings(as.integer(p[[1]])); a <- suppressWarnings(as.integer(p[[2]]))
  set(ad, j = paste0("reads_ref_alt_", s), value = paste0(r, ",", a))
  set(ad, j = paste0("VAF_", s), value = round(a / (r + a), 3))
}
ad <- ad[, c("CHROM", "POS", "REF", "ALT", paste0("reads_ref_alt_", SAMPLES), paste0("VAF_", SAMPLES)), with = FALSE]
su <- ad[su, on = .(CHROM, POS, REF, ALT)]
message("primary sites: ", nrow(su), "; without read counts in the masked file: ", su[is.na(get(paste0("VAF_", SAMPLES[1]))), .N])

# ---- 3. Differing sites per tumour x normal -----------------------------------
long <- rbindlist(lapply(TUMOURS, function(t) rbindlist(lapply(NORMALS, function(n) {
  x <- su[!is.na(get(t)) & !is.na(get(n))]
  x[, `:=`(tumour = t, normal = n, compared = .N, g_t = get(t), g_n = get(n))]
  x[g_t != g_n]
}))))
long[, pairing := fcase(normal == KI_NORMALS[tumour], "labelled", default = "other")]
long[, type := fcase(abs(g_t - g_n) == 2L, "opposite homozygotes",
                     g_n != 1L, "normal hom, tumour het",
                     default = "normal het, tumour hom")]

counts <- CJ(tumour = TUMOURS, normal = NORMALS, sorted = FALSE)
counts[, pairing := fifelse(normal == KI_NORMALS[tumour], "labelled", "other")]
counts[, compared := mapply(function(t, n) sum(!is.na(su[[t]]) & !is.na(su[[n]])), tumour, normal)]
k <- long[, .(differ = .N, opposite_homs = sum(type == "opposite homozygotes"),
              normal_hom_tumour_het = sum(type == "normal hom, tumour het"),
              normal_het_tumour_hom = sum(type == "normal het, tumour hom")), by = .(tumour, normal)]
counts <- k[counts, on = .(tumour, normal)]
for (cl in c("differ", "opposite_homs", "normal_hom_tumour_het", "normal_het_tumour_hom")) set(counts, which(is.na(counts[[cl]])), cl, 0L)
setcolorder(counts, c("tumour", "normal", "pairing", "compared"))
message("Differing / compared sites (slide: HCC_1-KI1 522 / 2,972; HCC_2-KI2 519 / 2,979; HCC_3-KI3 1 / 2,955):")
print(counts)

# ---- 4. Output tables ----------------------------------------------------------
pk <- if (file.exists(PATRIK)) fread(PATRIK, colClasses = list(character = "CHROM"))[, paste(CHROM, POS)] else character()
long[, in_patrik_list := paste(CHROM, POS) %chin% pk]
long[, locus := paste0(CHROM, ":", POS)]
for (s in SAMPLES) set(long, j = paste0("GT_", s), value = fifelse(is.na(long[[s]]), "masked", GT[long[[s]] + 1L]))
cols <- c("tumour", "normal", "pairing", "type", "locus", "CHROM", "POS", "REF", "ALT",
          paste0("GT_", SAMPLES), paste0("reads_ref_alt_", SAMPLES), paste0("VAF_", SAMPLES), "in_patrik_list")
res <- long[order(match(tumour, TUMOURS), match(normal, NORMALS), suppressWarnings(as.integer(sub("^chr", "", CHROM))), POS), cols, with = FALSE]
fwrite(counts, file.path(OUT, "counts_3x3.tsv"), sep = "\t")
fwrite(res, file.path(OUT, "discordant_sites_all_pairs.tsv"), sep = "\t")
for (t in TUMOURS) {
  x <- res[tumour == t & normal == KI_NORMALS[[t]]]
  fwrite(x, file.path(OUT, sprintf("discordant_%s_vs_%s.tsv", t, KI_NORMALS[[t]])), sep = "\t")
  message(sprintf("%s vs %s: %d differing sites (%s)", t, KI_NORMALS[[t]], nrow(x),
                  paste(x[, .N, by = type][, paste(N, type)], collapse = ", ")))
}
print(res[pairing == "labelled", .N, by = .(tumour, normal, type)])

# ---- 5. RUN_LOG.md ---------------------------------------------------------------
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_identity_discordant_sites"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Command:    Rscript ", SCRIPT),
  paste0("- Inputs:     ", file.path(IN, "sites_used_primary.tsv"), " (md5 ", tools::md5sum(file.path(IN, "sites_used_primary.tsv")),
         "); ", mf, " (md5 ", tools::md5sum(mf), "); Patrik list: ", if (length(pk)) PATRIK else "not found"),
  paste0("- Parameters: labelled ", paste(names(KI_NORMALS), KI_NORMALS, sep = "-", collapse = ", "), "; all normal genotypes compared"),
  paste0("- Tools:      ", system2(BCFTOOLS, "--version", stdout = TRUE)[1], "; data.table ", packageVersion("data.table")),
  paste0("- Outputs:    ", OUT),
  paste0("- Counts:     primary sites ", nrow(su), "; labelled pairs differing: ",
         paste(counts[pairing == "labelled", sprintf("%s-%s %d/%d", tumour, normal, differ, compared)], collapse = ", ")),
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
