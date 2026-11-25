#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_swap_checks.R - two extra checks of the setting-A sample swap, from the existing VCFs (no BAMs).
# mouse_wgs_somalier.R showed HCC_1 matches KI2 and HCC_2 matches KI1 (0 IBS0), HCC_3 matches KI3.
#
# Check 1 - Sprtn locus. Every variant in Sprtn +- FLANK, with per-sample GT / REF,ALT reads / VAF, from
#   the setting A FreeBayes and Mutect2 files (tumours + KI normals) and the setting B Mutect2 files (the WT
#   normals). The knock-in allele should be a variant present in all KI samples and absent from the WT
#   normals; a tumour must carry it (and any nearby donor-strain variants) like its true normal.
#   The KI mutation is not recorded in the project, so it is found from the data, not assumed. A large
#   knock-in (cassette, loxP) is not visible as a SNV/indel; then only nearby passenger variants inform.
#
# Check 2 - are the "somatic" calls of the labelled pairs germline in the other mouse? If HCC_1 (animal of
#   KI2) was called against KI1, many of its PASS "somatic" SNVs are KI2-animal germline variants: they
#   should be non-reference (het/hom VAF) in KI2. Per tumour x caller x KI normal: fraction of the tumour's
#   PASS SNVs that are germline in that normal, from the FreeBayes files (germline genotypes for every pair).
#   Expected if the swap is real: HCC_1 calls germline in KI2, HCC_2 calls germline in KI1, both high;
#   HCC_3 (correct pair) low in KI1 and KI2. Germline genotypes of a normal are only known where FreeBayes
#   wrote a record in that normal's own file; "absent" = not variant in that pair = most likely hom-ref.
#
# Blocks:
#   0  EDIT THIS BLOCK
#   1  setup, helpers
#   2  files
#   3  check 1: Sprtn locus
#   4  check 2: somatic calls germline in another normal
#   5  plots, outputs, RUN_LOG.md
#
# Env on ssi: r-mutation (data.table, ggplot2, bcftools). Run: Rscript scripts/mouse_wgs/mouse_wgs_swap_checks.R
# Output: results/mouse_wgs/mouse_wgs_swap_checks/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
DIR_A    <- file.path(REPO, "data/mouse/wgs/tumor_vs_normal/snv_indel")
DIR_B    <- file.path(REPO, "data/mouse/wgs/tumor_vs_joint_wt")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_swap_checks")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_swap_checks.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
# GRCm38 Ensembl 102 GTF (downloaded by scripts/cross_species/h2m 00_resources.R)
GTF      <- "/common/WORK/pstancl/references/cross_species/ensembl102/Mus_musculus.GRCm38.102.gtf.gz"
GENE     <- "Sprtn"
FLANK    <- 5000          # bp each side of the gene

TUMOURS    <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1", HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
                HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing
GENETIC    <- c(HCC_1 = "KI2", HCC_2 = "KI1", HCC_3 = "KI3")   # pairing from mouse_wgs_somalier.R
CALLERS    <- c("strelka", "mutect2", "freebayes")

MIN_DP      <- 20         # REF+ALT reads for a genotype in the other normal
HOM_REF_MAX <- 0.05
HET_MIN     <- 0.35       # VAF >= HET_MIN = germline non-reference (het or hom-alt)
MIN_QUAL_FB <- 20         # FreeBayes QUAL for the germline lookup

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS)
if (!file.exists(GTF)) stop("GTF not found: ", GTF, " (run scripts/cross_species/h2m 00_resources.R or set GTF)")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
B <- shQuote(BCFTOOLS)
bcf_version <- suppressWarnings(system2(BCFTOOLS, "--version", stdout = TRUE))[1]
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150, limitsize = FALSE)
  invisible(p)
}
# HCC_1..3 and KI1..3; any other sample (WT normals) keeps its name
short_name <- function(s) fifelse(grepl("KI_HCC_[0-9]", s), sub(".*KI_(HCC_[0-9]).*", "\\1", s),
                                  fifelse(grepl("_KI[0-9]$", s), sub(".*_(KI[0-9])$", "\\1", s), s))
header <- function(path) suppressWarnings(system2(BCFTOOLS, c("view", "-h", shQuote(path)), stdout = TRUE))
fmt_ids <- function(h) sub("^##FORMAT=<ID=([^,]+),.*", "\\1", grep("^##FORMAT=<ID=", h, value = TRUE))
samples_of <- function(h) strsplit(tail(grep("^#CHROM", h, value = TRUE), 1), "\t", fixed = TRUE)[[1]][-(1:9)]
TUMOUR_IDS <- names(TUMOURS); NORMAL_IDS <- unname(KI_NORMALS)

# ---- 2. Files ----------------------------------------------------------------
list_vcfs <- function(dir) {
  f <- list.files(dir, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf|vcf\\.gz)$")
  f[!grepl("union|merged", basename(f), ignore.case = TRUE)]
}
files <- rbind(data.table(setting = "A", path = list_vcfs(DIR_A)), data.table(setting = "B", path = list_vcfs(DIR_B)))
files[, tumour := NA_character_]
for (k in TUMOUR_IDS) files[grepl(TUMOURS[[k]], basename(path), fixed = TRUE), tumour := k]
files[, caller := NA_character_]
for (cl in CALLERS) files[is.na(caller) & grepl(cl, path, ignore.case = TRUE), caller := cl]
files[setting == "B" & is.na(caller), caller := "mutect2"]
hdrs <- setNames(lapply(files$path, header), files$path)
files[, has_GT_AD := vapply(path, function(p) all(c("GT", "AD") %in% fmt_ids(hdrs[[p]])), TRUE)]
files[, samples := vapply(path, function(p) paste(short_name(samples_of(hdrs[[p]])), collapse = ","), "")]
files[, chr_prefix := vapply(path, function(p) any(grepl("^##contig=<ID=chr", hdrs[[p]])), TRUE)]
print(files[, .(setting, tumour, caller, has_GT_AD, samples, file = basename(path))])
stopifnot(!anyNA(files[setting == "A", tumour]), !anyNA(files[setting == "A", caller]))

# ---- 3. Check 1: Sprtn locus -----------------------------------------------------
gene <- fread(cmd = paste("zcat", shQuote(GTF), "| awk -F'\\t' '$3==\"gene\" && $9 ~ /gene_name \"", GENE, "\";/'", sep = ""),
              header = FALSE, sep = "\t")
if (nrow(gene) != 1) stop(GENE, ": expected 1 gene record in the GTF, found ", nrow(gene))
G_CHR <- sub("^chr", "", gene$V1); G_START <- gene$V4; G_END <- gene$V5; G_STRAND <- gene$V7
sprintf("%s: %s:%d-%d (%s), window +-%d bp", GENE, G_CHR, G_START, G_END, G_STRAND, FLANK)

read_region <- function(i) {
  f <- files[i]
  if (!f$has_GT_AD) return(NULL)                       # e.g. Strelka (no GT/AD): skipped
  reg <- sprintf("%s%s:%d-%d", if (f$chr_prefix) "chr" else "", G_CHR, G_START - FLANK, G_END + FLANK)
  cmd <- paste(B, "view -t", reg, shQuote(f$path), "|", B, "query -f",
               shQuote("%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%SAMPLE=%GT=%AD]\\n"))
  x <- fread(cmd = cmd, sep = "\t", header = FALSE, colClasses = "character", fill = TRUE)
  if (!nrow(x)) return(NULL)
  setnames(x, 1:5, c("CHROM", "POS", "REF", "ALT", "FILTER"))
  long <- melt(x, id.vars = c("CHROM", "POS", "REF", "ALT", "FILTER"), value.name = "s")[nzchar(s) & !is.na(s)]
  long[, c("sample", "GT", "AD") := tstrsplit(s, "=", fixed = TRUE)]
  long[, ad := strsplit(AD, ",", fixed = TRUE)]
  long[, `:=`(REF_N = suppressWarnings(as.numeric(vapply(ad, `[`, "", 1))),
              ALT_N = suppressWarnings(as.numeric(vapply(ad, function(v) if (length(v) > 1) v[2] else NA_character_, ""))))]
  long[, .(setting = f$setting, caller = f$caller, file = basename(f$path), CHROM, POS = as.integer(POS), REF, ALT, FILTER,
           sample = short_name(sample), GT, REF_N, ALT_N)]
}
loc <- rbindlist(lapply(seq_len(nrow(files)), read_region))
if (!nrow(loc)) message("No records in the ", GENE, " window in any file with GT/AD")
if (nrow(loc)) {
  loc[, DP := REF_N + ALT_N]
  loc[, VAF := fifelse(DP > 0, ALT_N / DP, NA_real_)]
  loc[, in_gene := POS >= G_START & POS <= G_END]
  loc[, key := paste(CHROM, POS, REF, ALT, sep = ":")]
  # one value per variant x sample: the deepest record (FreeBayes and Mutect2 may both report it)
  best <- loc[order(-DP)][, .SD[1], by = .(key, sample)]
  best[, call := fcase(is.na(DP) | DP < MIN_DP, "low_depth", VAF <= HOM_REF_MAX, "ref",
                       VAF >= HET_MIN & VAF < 0.95, "het", VAF >= 0.95, "hom_alt", default = "ambiguous")]
  # Wide: VAF (depth) per sample; WT normals are the samples that are neither KI nor tumour
  ki_samples <- c(TUMOUR_IDS, NORMAL_IDS)
  wt_samples <- setdiff(unique(best$sample), ki_samples)
  locus_wide <- dcast(best, key + in_gene ~ sample, value.var = "call")
  locus_vaf  <- dcast(best[, .(key, sample, v = sprintf("%.2f (%d)", VAF, as.integer(DP)))], key ~ sample, value.var = "v")
  # Candidate knock-in / donor-strain variants: non-ref in every KI sample seen, ref in every WT normal seen
  cand <- best[, .(ki_seen = sum(sample %chin% ki_samples & call != "low_depth"),
                   ki_nonref = sum(sample %chin% ki_samples & call %chin% c("het", "hom_alt")),
                   wt_seen = sum(sample %chin% wt_samples & call != "low_depth"),
                   wt_nonref = sum(sample %chin% wt_samples & call %chin% c("het", "hom_alt"))), by = key]
  cand[, ki_specific := ki_seen >= 4 & ki_nonref == ki_seen & wt_nonref == 0]
  print(cand[order(-ki_specific, -ki_nonref)])
  # Per tumour: Sprtn-window genotypes concordant with the labelled vs the genetic normal
  conc <- rbindlist(lapply(TUMOUR_IDS, function(t) rbindlist(lapply(unique(c(KI_NORMALS[[t]], GENETIC[[t]])), function(n) {
    a <- best[sample == t & call %chin% c("ref", "het", "hom_alt"), .(key, ct = call)]
    b <- best[sample == n & call %chin% c("ref", "het", "hom_alt"), .(key, cn = call)]
    m <- a[b, on = "key", nomatch = NULL]
    data.table(tumour = t, normal = n, role = fifelse(n == KI_NORMALS[[t]], "labelled", "genetic"),
               sites_both = nrow(m), same_genotype = sum(m$ct == m$cn),
               opposite_homs = sum((m$ct == "ref" & m$cn == "hom_alt") | (m$ct == "hom_alt" & m$cn == "ref")))
  }))))
  print(conc)
} else { locus_wide <- data.table(); locus_vaf <- data.table(); cand <- data.table(); conc <- data.table(); best <- data.table() }

# ---- 4. Check 2: somatic calls germline in another normal --------------------------
# 4a. PASS SNVs per tumour x caller (setting A, labelled pairs)
pass_snvs <- function(path) {
  cmd <- paste(B, "view -f PASS -v snps", shQuote(path), "|", B, "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT\\n"))
  x <- fread(cmd = cmd, sep = "\t", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT"), colClasses = "character")
  x[, key := paste(sub("^chr", "", CHROM), POS, REF, ALT, sep = ":")][]
}
som <- rbindlist(lapply(which(files$setting == "A"), function(i)
  pass_snvs(files$path[i])[, `:=`(tumour = files$tumour[i], caller = files$caller[i])]))
som <- unique(som, by = c("tumour", "caller", "key"))
som[, .N, by = .(tumour, caller)][order(tumour, caller)]

# 4b. Germline genotypes of every KI normal from its own FreeBayes file (all FILTER values)
fb <- files[setting == "A" & caller == "freebayes"]
stopifnot(nrow(fb) >= length(TUMOUR_IDS), all(fb$has_GT_AD))
germ <- rbindlist(lapply(fb$path, function(p) {
  smp <- short_name(samples_of(hdrs[[p]]))
  q <- paste0("%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n")
  x <- fread(cmd = paste(B, "view -v snps -m2 -M2 -i", shQuote(paste0("QUAL>=", MIN_QUAL_FB)), shQuote(p), "|", B, "query -f", shQuote(q)),
             sep = "\t", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", smp), colClasses = "character")
  n <- intersect(smp, NORMAL_IDS)
  if (length(n) != 1) stop(basename(p), ": expected one KI normal, found ", paste(smp, collapse = ", "))
  ad <- tstrsplit(x[[n]], ",", fixed = TRUE)
  data.table(normal = n, key = paste(sub("^chr", "", x$CHROM), x$POS, x$REF, x$ALT, sep = ":"),
             REF_N = suppressWarnings(as.numeric(ad[[1]])), ALT_N = suppressWarnings(as.numeric(if (length(ad) > 1) ad[[2]] else NA)))
}))
germ[, DP := REF_N + ALT_N]
germ[, VAF := fifelse(DP > 0, ALT_N / DP, NA_real_)]
germ <- unique(germ[order(-DP)], by = c("normal", "key"))
germ[, .N, by = normal]

# 4c. Each tumour's PASS SNVs looked up in each KI normal
look <- rbindlist(lapply(NORMAL_IDS, function(n) {
  g <- germ[normal == n]
  x <- copy(som)[, normal := n]
  x[g, on = "key", `:=`(DP = i.DP, VAF = i.VAF)]
  x[, status := fcase(is.na(DP), "absent from normal's file", DP < MIN_DP, "low depth",
                      VAF >= HET_MIN, "germline non-ref", VAF <= HOM_REF_MAX, "ref", default = "ambiguous")]
  x
}))
germline_overlap <- look[, .(pass_snvs = .N, germline_nonref = sum(status == "germline non-ref"),
                             ref = sum(status == "ref"), absent = sum(status == "absent from normal's file"),
                             low_depth = sum(status == "low depth"), ambiguous = sum(status == "ambiguous")),
                         by = .(tumour, caller, normal)]
germline_overlap[, pct_germline_in_normal := round(100 * germline_nonref / pass_snvs, 1)]
germline_overlap[, role := fcase(normal == KI_NORMALS[tumour], "labelled normal (called against)",
                                 normal == GENETIC[tumour], "genetic normal", default = "other normal")]
setorder(germline_overlap, tumour, caller, normal)
germline_overlap
# VAF of the tumour's calls in its genetic normal: germline het ~0.5, hom ~1
vaf_genetic <- look[normal == GENETIC[tumour] & normal != KI_NORMALS[tumour] & !is.na(VAF) & DP >= MIN_DP]

# ---- 5. Plots, outputs, RUN_LOG.md ------------------------------------------------
p_overlap <- ggplot(germline_overlap, aes(normal, pct_germline_in_normal, fill = role)) +
  geom_col(width = 0.7) + geom_text(aes(label = paste0(germline_nonref, "/", pass_snvs)), vjust = -0.3, size = 3) +
  facet_grid(caller ~ tumour) +
  scale_fill_manual(values = c(`labelled normal (called against)` = "grey60", `genetic normal` = "#C62828", `other normal` = "#1565C0"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "KI normal", y = "% of the tumour's PASS SNVs that are germline non-ref in this normal",
       title = "Setting A 'somatic' SNVs that are germline variants of another mouse",
       subtitle = sprintf("Germline = VAF >= %.2f with >= %d reads in that normal's FreeBayes file (QUAL >= %d)", HET_MIN, MIN_DP, MIN_QUAL_FB)) +
  theme_bw() + theme(legend.position = "bottom")
save_plot(p_overlap, "somatic_calls_germline_in_other_normal", 11, 8)
if (nrow(vaf_genetic)) {
  p_vaf <- ggplot(vaf_genetic, aes(VAF)) + geom_histogram(binwidth = 0.02) + facet_grid(caller ~ tumour + normal, scales = "free_y") +
    labs(x = "VAF in the genetic normal", y = "PASS SNVs", title = "VAF of the swapped tumours' 'somatic' calls in their genetic normal",
         subtitle = "Germline variants peak at 0.5 (het) and 1 (hom); true somatic mutations sit at 0") + theme_bw()
  save_plot(p_vaf, "vaf_in_genetic_normal", 10, 7)
}
if (nrow(best)) {
  tile <- copy(best)
  tile[, sample := factor(sample, c(TUMOUR_IDS, NORMAL_IDS, sort(setdiff(unique(sample), c(TUMOUR_IDS, NORMAL_IDS)))))]
  tile[, key := factor(key, unique(key[order(as.integer(sub("^[^:]*:([0-9]+):.*", "\\1", key)))]))]
  p_loc <- ggplot(tile, aes(sample, key, fill = VAF)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f\n%d", VAF, as.integer(DP))), size = 2.4) +
    scale_fill_gradient(low = "white", high = "#C62828", limits = c(0, 1), na.value = "grey85") +
    labs(x = NULL, y = paste(GENE, "window variant"), title = sprintf("%s +- %d bp: VAF (and depth) per sample", GENE, FLANK),
         subtitle = "Knock-in / donor-strain alleles: non-ref in all KI samples, ref in WT. A tumour must match its true normal.") +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_loc, "sprtn_locus_vaf", 4 + 0.8 * uniqueN(tile$sample), 2 + 0.25 * uniqueN(tile$key))
}

fwrite(files[, .(setting, tumour, caller, has_GT_AD, samples, file = basename(path))], file.path(OUT, "files.tsv"), sep = "\t")
if (nrow(loc)) {
  fwrite(loc, file.path(OUT, "sprtn_locus_records.tsv"), sep = "\t")
  fwrite(locus_vaf, file.path(OUT, "sprtn_locus_vaf_wide.tsv"), sep = "\t")
  fwrite(locus_wide, file.path(OUT, "sprtn_locus_calls_wide.tsv"), sep = "\t")
  fwrite(cand, file.path(OUT, "sprtn_locus_ki_specific.tsv"), sep = "\t")
  fwrite(conc, file.path(OUT, "sprtn_locus_tumour_vs_normals.tsv"), sep = "\t")
}
fwrite(germline_overlap, file.path(OUT, "somatic_calls_germline_in_other_normal.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
git_hash <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)),
                     error = function(e) "unknown")
if (!length(git_hash)) git_hash <- "unknown"
md_table <- function(x) { x <- as.data.table(x); if (!nrow(x)) return("(none)")
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))) }
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_swap_checks"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ", ", bcf_version, ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  sprintf("- Parameters: GENE=%s (%s:%d-%d) FLANK=%d MIN_DP=%d HOM_REF_MAX=%s HET_MIN=%s MIN_QUAL_FB=%d; genetic pairing %s",
          GENE, G_CHR, G_START, G_END, FLANK, MIN_DP, HOM_REF_MAX, HET_MIN, MIN_QUAL_FB, paste(names(GENETIC), GENETIC, sep = "-", collapse = ", ")),
  paste0("- Inputs:     ", nrow(files), " VCF/BCFs from ", DIR_A, " and ", DIR_B, "; GTF ", GTF),
  paste0("- Outputs:    ", OUT),
  "", "### Check 2: tumour PASS SNVs that are germline non-ref in each KI normal", "", md_table(germline_overlap),
  "", "### Check 1: Sprtn window, tumour vs labelled / genetic normal", "", md_table(conc),
  "", "### Notes",
  "- Germline genotypes of a normal exist only where FreeBayes wrote a record in that normal's own pair file; 'absent' most likely means hom-ref in that normal.",
  "- SNV keys are CHROM:POS:REF:ALT; an SNV that FreeBayes wrote inside a complex/MNP record is missed (counted absent).",
  "- Strelka files have no GT/AD and are not used for check 1; their PASS SNVs are used in check 2.",
  "- A large knock-in (cassette, loxP) is not visible as SNV/indel calls; then only passenger variants near the gene inform.",
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
