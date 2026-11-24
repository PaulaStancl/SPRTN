#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24
#
# mouse_wgs_somalier.R - sample identity with somalier (Pedersen et al., Genome Medicine 2020), from VCFs.
#
# somalier genotypes every sample at a fixed list of germline SNP sites and reports, per pair,
# relatedness, IBS0 (opposite homozygotes), IBS2 (identical genotypes) and homozygous concordance.
# Same animal: relatedness ~1 and IBS0 ~0. It works from VCF/BCF (no BAMs needed).
#
# somalier only ships HUMAN site lists, so the mouse site list is built here from the FreeBayes calls
# (setting A), keeping sites that look GERMLINE, not somatic:
#   - SNV, biallelic, autosome, FreeBayes QUAL >= MIN_QUAL
#   - non-reference (VAF >= 0.35) in at least one KI NORMAL liver   (present in normal tissue)
#   - no KI normal with an ambiguous VAF between 0.05 and 0.35        (no tumour-in-normal / artefact pattern)
#   - REF+ALT reads >= MIN_DP in at least MIN_SAMPLES_DP samples
#   - at most one site per SITE_SPACING bp (keeps the best-covered one), so one inherited block does not dominate
# The site list is selected from the same samples that are then compared. That is fine for detecting
# swaps (we ask which samples agree), but somalier's relatedness values assume population allele
# frequencies and are only rough here; IBS0 and IBS2 are the robust numbers.
#
# Blocks:
#   0  EDIT THIS BLOCK
#   1  setup, tool checks
#   2  build the germline site list from the FreeBayes files -> sites.vcf.gz
#   3  somalier extract (one run per FreeBayes file; all samples in it)
#   4  somalier relate
#   5  read the pairs table; per-tumour summary; plots
#   6  outputs, RUN_LOG.md
#
# Install somalier once (separate env, avoids conflicts with r-mutation):
#   micromamba create -y -p /common/WORK/pstancl/envs/somalier -c conda-forge -c bioconda somalier
# or the static binary: https://github.com/brentp/somalier/releases (chmod +x)
#
# Env on ssi: r-mutation (data.table, ggplot2, bcftools) + the somalier binary (SOMALIER below).
# Output: results/mouse_wgs/mouse_wgs_somalier/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
DIR_A    <- file.path(REPO, "data/mouse/wgs/tumor_vs_normal/snv_indel")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_somalier")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_somalier.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
SOMALIER <- "/common/WORK/pstancl/envs/somalier/bin/somalier"
# GRCm38 FASTA (+ .fai) with the same contig names as the VCFs ("1", "2", ...)
FASTA    <- "/common/WORK/pstancl/references/cross_species/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa"

TUMOURS    <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1", HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
                HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing

MIN_QUAL       <- 30
MIN_DP         <- 30
MIN_SAMPLES_DP <- 4       # of the 6 samples
SITE_SPACING   <- 20000   # bp

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# Console as wide as the terminal (falls back to 80 when there is no terminal, e.g. batch Rscript)
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS)
if (!file.exists(SOMALIER)) stop("somalier not found at ", SOMALIER, " - see the install lines in the header")
if (!file.exists(paste0(FASTA, ".fai"))) stop("FASTA index missing: ", FASTA, ".fai")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
B <- shQuote(BCFTOOLS)
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  invisible(p)
}
short_name <- function(s) fifelse(grepl("KI_HCC_[0-9]", s), sub(".*KI_(HCC_[0-9]).*", "\\1", s), sub(".*_(KI[0-9])$", "\\1", s))
run <- function(cmd) { message("$ ", cmd); st <- system(cmd); if (st != 0) stop("failed (exit ", st, "): ", cmd) }
system2(SOMALIER, "--version", stdout = TRUE)[1]

# ---- 2. Germline site list from the FreeBayes files ---------------------------
fb_files <- list.files(DIR_A, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf|vcf\\.gz)$")
fb_files <- fb_files[grepl("freebayes", fb_files, ignore.case = TRUE) & !grepl("union|merged", fb_files, ignore.case = TRUE)]
basename(fb_files)
stopifnot(length(fb_files) == length(TUMOURS))

read_fb <- function(path) {
  smp <- system2(BCFTOOLS, c("query", "-l", shQuote(path)), stdout = TRUE)
  cmd <- paste(B, "view -v snps -m2 -M2 -i", shQuote(paste0("QUAL>=", MIN_QUAL)), shQuote(path), "|",
               B, "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n"))
  x <- fread(cmd = cmd, sep = "\t", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", paste0("AD_", seq_along(smp))),
             colClasses = "character")
  x <- x[sub("^chr", "", CHROM) %chin% as.character(1:19)]
  rbindlist(lapply(seq_along(smp), function(j) {
    ad <- tstrsplit(x[[paste0("AD_", j)]], ",", fixed = TRUE)
    data.table(CHROM = x$CHROM, POS = as.integer(x$POS), REF = x$REF, ALT = x$ALT, sample = short_name(smp[j]),
               REF_N = suppressWarnings(as.numeric(ad[[1]])), ALT_N = suppressWarnings(as.numeric(if (length(ad) > 1) ad[[2]] else NA)))
  }))
}
ad <- rbindlist(lapply(fb_files, read_fb))
ad[, `:=`(DPA = REF_N + ALT_N, VAF = fifelse(REF_N + ALT_N > 0, ALT_N / (REF_N + ALT_N), NA_real_))]
ad <- unique(ad[order(-DPA)], by = c("CHROM", "POS", "REF", "ALT", "sample"))
ad[, .N, by = sample]

# Per site: depth in enough samples; germline evidence in a KI normal; no ambiguous normal
site_stats <- ad[, .(n_dp = sum(DPA >= MIN_DP, na.rm = TRUE),
                     normal_nonref = any(sample %chin% KI_NORMALS & DPA >= MIN_DP & VAF >= 0.35, na.rm = TRUE),
                     normal_ambiguous = any(sample %chin% KI_NORMALS & DPA >= MIN_DP & VAF > 0.05 & VAF < 0.35, na.rm = TRUE),
                     min_dp = min(DPA, na.rm = TRUE)),
                 by = .(CHROM, POS, REF, ALT)]
sites <- site_stats[n_dp >= MIN_SAMPLES_DP & normal_nonref & !normal_ambiguous]
c(candidate_sites = nrow(site_stats), germline_like = nrow(sites))
# Thin: best-covered site per SITE_SPACING window
sites[, win := POS %/% SITE_SPACING]
sites <- sites[order(-min_dp)][, .SD[1], by = .(CHROM, win)][order(as.integer(sub("^chr", "", CHROM)), POS)]
c(sites_after_thinning = nrow(sites))

# Write sites.vcf.gz (contigs from the FASTA index)
fai <- fread(paste0(FASTA, ".fai"), header = FALSE, select = 1:2, col.names = c("ID", "length"))
sites_vcf <- file.path(OUT, "sites.vcf")
writeLines(c("##fileformat=VCFv4.2", fai[, sprintf("##contig=<ID=%s,length=%s>", ID, length)],
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), sites_vcf)
fwrite(sites[, .(CHROM, POS, ID = ".", REF, ALT, QUAL = ".", FILTER = "PASS", INFO = ".")], sites_vcf,
       sep = "\t", col.names = FALSE, append = TRUE)
run(paste(B, "view -Oz -o", shQuote(paste0(sites_vcf, ".gz")), shQuote(sites_vcf)))
run(paste(B, "index -t -f", shQuote(paste0(sites_vcf, ".gz"))))
unlink(sites_vcf)

# ---- 3. somalier extract ---------------------------------------------------------
ext_dir <- file.path(OUT, "extracted"); dir.create(ext_dir, showWarnings = FALSE)
for (f in fb_files)
  run(paste(shQuote(SOMALIER), "extract -d", shQuote(ext_dir), "--sites", shQuote(paste0(sites_vcf, ".gz")),
            "-f", shQuote(FASTA), shQuote(f)))
list.files(ext_dir)

# ---- 4. somalier relate --------------------------------------------------------
prefix <- file.path(OUT, "somalier")
run(paste(shQuote(SOMALIER), "relate -o", shQuote(prefix), paste(shQuote(list.files(ext_dir, full.names = TRUE)), collapse = " ")))

# ---- 5. Results ----------------------------------------------------------------
pairs <- fread(paste0(prefix, ".pairs.tsv"))
setnames(pairs, 1, "sample_a")
pairs[, `:=`(a = short_name(sample_a), b = short_name(sample_b))]
pairs[, ibs0_per_1000 := round(1000 * ibs0 / n, 2)]
pairs[order(-relatedness), .(a, b, relatedness, ibs0, ibs2, hom_concordance, n, ibs0_per_1000)]

# Tumour vs KI normal
tn <- rbind(pairs[grepl("^HCC", a) & grepl("^KI", b), .(tumour = a, normal = b, relatedness, ibs0, ibs2, n, ibs0_per_1000, hom_concordance)],
            pairs[grepl("^HCC", b) & grepl("^KI", a), .(tumour = b, normal = a, relatedness, ibs0, ibs2, n, ibs0_per_1000, hom_concordance)])
tn[, labelled := KI_NORMALS[tumour] == normal]
verdict <- tn[order(tumour, -relatedness)][, .(labelled_normal = KI_NORMALS[tumour[1]], best_normal = normal[1],
                                               best_relatedness = relatedness[1], labelled_relatedness = relatedness[normal == KI_NORMALS[tumour[1]]][1],
                                               best_ibs0 = ibs0[1], labelled_ibs0 = ibs0[normal == KI_NORMALS[tumour[1]]][1]), by = tumour]
verdict[, consistent := best_normal == labelled_normal]
verdict

p_rel <- ggplot(tn, aes(normal, relatedness, fill = labelled)) +
  geom_col(width = 0.7) + geom_text(aes(label = sprintf("%.2f\nIBS0 %d", relatedness, ibs0)), vjust = -0.2, size = 3.5) +
  facet_wrap(~ tumour, labeller = labeller(tumour = function(x) paste("Tumour", x))) +
  scale_fill_manual(values = c(`TRUE` = "#C62828", `FALSE` = "grey60"), labels = c(`TRUE` = "labelled as its normal", `FALSE` = "other normal")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "normal liver", y = "somalier relatedness", fill = NULL,
       title = "somalier: which normal liver matches each tumour?",
       subtitle = sprintf("%s germline-like sites. Highest bar (and lowest IBS0) = same animal", format(nrow(sites), big.mark = ","))) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
save_plot(p_rel, "somalier_relatedness_per_tumour", 10, 5)

# IBS0 vs IBS2 for all pairs (somalier's standard view); same animal = right, near IBS0 = 0
pairs[, pair := paste(a, "-", b)]
pairs[, type := fcase(paste(a, b) %chin% c("HCC_1 KI1", "KI1 HCC_1", "HCC_2 KI2", "KI2 HCC_2", "HCC_3 KI3", "KI3 HCC_3"), "labelled pair",
                      grepl("^HCC", a) != grepl("^HCC", b), "other tumour-normal", default = "same tissue type")]
lab_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) ggrepel::geom_text_repel(aes(label = pair), size = 3) else
  geom_text(aes(label = pair), size = 3, vjust = -0.8)
p_ibs <- ggplot(pairs, aes(ibs2 / n, ibs0 / n, colour = type)) +
  geom_point(size = 3) + lab_layer +
  scale_colour_manual(values = c(`labelled pair` = "#C62828", `other tumour-normal` = "grey30", `same tissue type` = "#1565C0")) +
  labs(x = "IBS2 / n (identical genotypes)", y = "IBS0 / n (opposite homozygotes)", colour = NULL,
       title = "somalier IBS0 vs IBS2, all sample pairs", subtitle = "Same animal: bottom right") + theme_bw()
save_plot(p_ibs, "somalier_ibs0_vs_ibs2", 9, 7)

# ---- 6. Outputs --------------------------------------------------------------
fwrite(sites[, .(CHROM, POS, REF, ALT, min_dp)], file.path(OUT, "sites_used.tsv"), sep = "\t")
fwrite(tn, file.path(OUT, "tumour_vs_normal.tsv"), sep = "\t")
fwrite(verdict, file.path(OUT, "verdict.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_somalier"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ", ",
         system2(SOMALIER, "--version", stdout = TRUE)[1], ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  sprintf("- Parameters: MIN_QUAL=%d MIN_DP=%d MIN_SAMPLES_DP=%d SITE_SPACING=%d; sites: nonref (VAF>=0.35) in >=1 KI normal, no KI normal 0.05-0.35",
          MIN_QUAL, MIN_DP, MIN_SAMPLES_DP, SITE_SPACING),
  paste0("- Inputs:     ", paste(basename(fb_files), collapse = ", "), "; FASTA ", FASTA),
  paste0("- Sites:      ", nrow(sites), " (sites.vcf.gz)"),
  paste0("- Verdict:    ", paste(verdict[, sprintf("%s best=%s (rel %.2f; labelled %s %.2f)", tumour, best_normal, best_relatedness,
                                                   labelled_normal, labelled_relatedness)], collapse = "; ")),
  paste0("- Outputs:    ", OUT, " (somalier.html, somalier.pairs.tsv, plots)"),
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Open ", prefix, ".html for somalier's interactive report")
