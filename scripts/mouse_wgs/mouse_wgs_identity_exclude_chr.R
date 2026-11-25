#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_identity_exclude_chr.R - sample matching with the chromosome of the engineered gene removed.
# Control asked for in review: all KI mice carry the same Sprtn knock-in, but DNA flanking the construct
# (ES-cell / donor-strain passenger region on chr8) can differ between animals. If that region drove the
# HCC_1-KI2 / HCC_2-KI1 match, the match would weaken or disappear without chr8.
#
# Scenarios (EXCLUDE): "all"     = full primary site set (sanity check: must reproduce mouse_wgs_somalier.R)
#                      "no_chr8" = the same sites minus chromosome 8 (Sprtn, GRCm38 chr8)
# Methods on each scenario, all on the SAME masked genotypes as mouse_wgs_somalier.R (primary combo A_20kb):
#   custom    normal-hom discordance (tumour genotype != homozygous normal genotype) + paired exact sign
#             test best vs second normal - the primary metric of mouse_wgs_somalier.R
#   somalier  extract + relate on the reduced sites VCF; best = lowest IBS0 rate, then highest IBS2 rate
#   gtcheck   bcftools gtcheck -u GT,GT -E 0 on the reduced targets; best = lowest discordance rate
#   KING      plink2 --make-king-table on the reduced sites; best = highest kinship, then lowest IBS0
# Not independent evidence (same FreeBayes calls); it tests whether chr8 carries the signal.
# The BAM-based check (direct genotyping of a fixed SNP panel in all six BAMs) needs the BAMs - not here.
#
# Inputs (from mouse_wgs_somalier.R, run it first):
#   sites_used_primary.tsv, masked/masked_6samples.bcf (+ .csi), sites/sites_A_20kb.vcf.gz,
#   jackknife.tsv (its row for the dropped chromosome is shown next to this result)
# Env on ssi: r-mutation (data.table, ggplot2, bcftools 1.24); somalier and plink2 as in the sibling scripts.
# Run: Rscript scripts/mouse_wgs/mouse_wgs_identity_exclude_chr.R
# Output: results/mouse_wgs/mouse_wgs_identity_exclude_chr/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
SOM_OUT  <- file.path(REPO, "results/mouse_wgs/mouse_wgs_somalier")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_identity_exclude_chr")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_identity_exclude_chr.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
SOMALIER <- Sys.which("somalier"); if (!nzchar(SOMALIER)) SOMALIER <- "/common/WORK/pstancl/tools/somalier/somalier"
PLINK2   <- Sys.which("plink2");   if (!nzchar(PLINK2))   PLINK2   <- "/common/WORK/pstancl/envs/plink2/bin/plink2"
FASTA    <- "/common/WORK/pstancl/references/cross_species/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa"
COMBO    <- "A_20kb"                                    # primary panel x thinning of mouse_wgs_somalier.R
EXCLUDE  <- list(all = character(), no_chr8 = "8")      # chromosome names without "chr"

TUMOUR_IDS <- c("HCC_1", "HCC_2", "HCC_3")
NORMAL_IDS <- c("KI1", "KI2", "KI3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing
HYPOTHESIS <- c(HCC_1 = "KI2", HCC_2 = "KI1", HCC_3 = "KI3")   # swap hypothesis
MIN_DP <- 30; HET_MIN <- 0.35; AB_HOM_CUTOFF <- 0.05 + 1e-6    # as in mouse_wgs_somalier.R
PAIRED_ALPHA <- 0.05

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
started <- Sys.time()
SAMPLES <- c(TUMOUR_IDS, NORMAL_IDS)
GT_TSV  <- file.path(SOM_OUT, "sites_used_primary.tsv")
MASKED  <- file.path(SOM_OUT, "masked", "masked_6samples.bcf")
SITES   <- file.path(SOM_OUT, "sites", paste0("sites_", COMBO, ".vcf.gz"))
if (!file.exists(GT_TSV)) stop("not found: ", GT_TSV, " - run mouse_wgs_somalier.R first")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NOTES <- character()
note <- function(...) { m <- paste0(...); NOTES <<- c(NOTES, m); message("NOTE: ", m) }
run <- function(cmd) { message("$ ", cmd); st <- system(cmd); if (st != 0) stop("failed (exit ", st, "): ", cmd) }
capture <- function(bin, args = character()) suppressWarnings(system2(bin, args, stdout = TRUE, stderr = TRUE))
HAVE_BCF  <- file.exists(BCFTOOLS) && file.exists(MASKED) && file.exists(SITES)
HAVE_SOM  <- HAVE_BCF && file.exists(SOMALIER) && file.exists(FASTA)
HAVE_KING <- file.exists(PLINK2) && file.exists(MASKED)
if (!HAVE_BCF)  note("bcftools, ", MASKED, " or ", SITES, " missing: somalier and gtcheck skipped (custom only)")
if (HAVE_BCF && !HAVE_SOM) note("somalier or FASTA missing: somalier skipped")
if (!HAVE_KING) note("plink2 not found at ", PLINK2, ": KING skipped")
versions <- c(bcftools = if (HAVE_BCF) capture(BCFTOOLS, "--version")[1] else NA,
              somalier = if (HAVE_SOM) grep("somalier version", capture(SOMALIER), value = TRUE)[1] else NA,
              plink2   = if (HAVE_KING) capture(PLINK2, "--version")[1] else NA)
print(versions)

G <- fread(GT_TSV, colClasses = list(character = "CHROM"))
stopifnot(all(SAMPLES %in% names(G)))
G[, chr_n := sub("^chr", "", CHROM)]
message(nrow(G), " primary sites on ", uniqueN(G$chr_n), " chromosomes; per chromosome:")
print(G[, .N, by = chr_n][order(as.integer(chr_n))])

# ---- 2. Custom normal-hom discordance + paired sign test ------------------------------
# Same definitions as mouse_wgs_somalier.R (pair_counts / paired_counts): genotypes 0/1/2, NA = masked
tn_counts <- function(d, t, n) {
  x <- d[[t]]; a <- d[[n]]; inf <- !is.na(x) & !is.na(a) & a != 1L
  list(n_inf = sum(inf), disc = sum(inf & x != a), ibs0 = sum(!is.na(x) & !is.na(a) & abs(x - a) == 2L))
}
paired_p <- function(d, t, nj, nk) {
  x <- d[[t]]; a <- d[[nj]]; b <- d[[nk]]
  ok <- !is.na(x) & !is.na(a) & !is.na(b) & (a != 1L | b != 1L)
  dj <- ok & a != 1L & x != a; dk <- ok & b != 1L & x != b
  nb <- sum(dj & !dk); nc <- sum(!dj & dk)
  c(disc_best_only = nb, disc_second_only = nc, p = if (nb + nc > 0) binom.test(nb, nb + nc)$p.value else 1)
}
custom_run <- function(d, sc) {
  m <- rbindlist(lapply(TUMOUR_IDS, function(t) rbindlist(lapply(NORMAL_IDS, function(n) {
    pc <- tn_counts(d, t, n)
    ci <- if (pc$n_inf > 0) binom.test(pc$disc, pc$n_inf)$conf.int else c(NA, NA)
    data.table(scenario = sc, method = "custom", tumour = t, normal = n, n_sites = pc$n_inf, value = pc$disc / pc$n_inf,
               value_lo = ci[1], value_hi = ci[2], count = pc$disc, sort1 = pc$disc / pc$n_inf, sort2 = 0)
  }))))
  m
}

# ---- 3. somalier, gtcheck, KING on the reduced site set ---------------------------------
gt_usage <- if (HAVE_BCF) capture(BCFTOOLS, c("gtcheck", "-h")) else character()   # -h: never reads stdin
has_opt <- function(o) any(grepl(o, gt_usage, fixed = TRUE))
GT_ARGS <- if (!HAVE_BCF) NULL else if (has_opt("-u, --use"))
  c("-u GT,GT", if (has_opt("-E, --error-probability")) "-E 0", if (has_opt("--no-HWE-prob")) "--no-HWE-prob") else
  if (has_opt("-G, --GTs-only")) "-G 1" else { note("gtcheck lists neither -u nor -G: skipped"); NULL }

contig_of <- function(chrs) unique(G[chr_n %chin% chrs, CHROM])   # actual contig spelling (8 or chr8)

somalier_run <- function(sc, d_dir, ex) {
  sites <- SITES
  if (length(ex)) {
    sites <- file.path(d_dir, "sites.vcf.gz")
    run(paste(shQuote(BCFTOOLS), "view -t", shQuote(paste0("^", paste(contig_of(ex), collapse = ","))),
              "-Oz -o", shQuote(sites), shQuote(SITES)))
    run(paste(shQuote(BCFTOOLS), "index -t -f", shQuote(sites)))
  }
  ex_dir <- file.path(d_dir, "extracted"); unlink(ex_dir, recursive = TRUE); dir.create(ex_dir, recursive = TRUE)
  run(paste(shQuote(SOMALIER), "extract -d", shQuote(ex_dir), "-s", shQuote(sites), "-f", shQuote(FASTA), shQuote(MASKED)))
  got <- list.files(ex_dir, pattern = "\\.somalier$", full.names = TRUE)
  if (length(got) != 6) stop("somalier extract ", sc, ": expected 6 .somalier files, got ", length(got))
  prefix <- file.path(d_dir, "somalier")
  run(paste0("SOMALIER_AB_HOM_CUTOFF=", format(AB_HOM_CUTOFF, scientific = FALSE), " ", shQuote(SOMALIER), " relate -s ",
             shQuote(sites), " -d ", MIN_DP, " --min-ab ", HET_MIN, " -o ", shQuote(prefix), " ",
             paste(shQuote(got), collapse = " "), " > ", shQuote(paste0(prefix, ".log")), " 2>&1"))
  p <- fread(paste0(prefix, ".pairs.tsv")); setnames(p, 1, sub("^#", "", names(p)[1]))
  need <- c("sample_a", "sample_b", "ibs0", "ibs2", "n", "concordance")
  if (!all(need %in% names(p))) stop("somalier pairs.tsv lacks: ", paste(setdiff(need, names(p)), collapse = ", "))
  p[, .(s1 = sample_a, s2 = sample_b, n_sites = n, value = ibs0 / n, count = ibs0, ibs2_rate = ibs2 / n, concordance)]
}
gtcheck_run <- function(sc, d_dir, tg) {
  raw <- file.path(d_dir, "gtcheck.txt")
  run(paste(shQuote(BCFTOOLS), "gtcheck", paste(GT_ARGS, collapse = " "), "-T", shQuote(tg), shQuote(MASKED),
            ">", shQuote(raw), "2>", shQuote(paste0(raw, ".stderr"))))
  l <- readLines(raw)
  hdr <- grep("^#\\s*(DCv2|DC)\\b.*\\[2\\]", l, value = TRUE, perl = TRUE)[1]
  if (!is.na(hdr)) {
    tag <- sub("^#\\s*(DCv2|DC)\\b.*", "\\1", hdr, perl = TRUE)
    m <- regmatches(hdr, gregexpr("\\[([0-9]+)\\]([^\\t\\[]+)", hdr, perl = TRUE))[[1]]   # perl: \\t is a tab inside [ ]
    idx <- setNames(as.integer(sub("^\\[([0-9]+)\\].*", "\\1", m)), tolower(trimws(sub("[,;:]\\s*$", "", trimws(sub("^\\[[0-9]+\\]", "", m))))))
    find <- function(p) { i <- grep(p, names(idx)); if (length(i)) idx[[i[1]]] else NA_integer_ }
    ii <- c(find("query"), find("genotyped"), find("^discordance"), find("sites compared|number of sites"))
    if (anyNA(ii)) { note("gtcheck ", sc, ": header not recognised - raw kept in ", raw); return(NULL) }
    x <- as.data.table(do.call(rbind, strsplit(grep(paste0("^", tag, "\t"), l, value = TRUE), "\t", fixed = TRUE)))
    x <- data.table(s1 = x[[ii[1]]], s2 = x[[ii[2]]], count = as.numeric(x[[ii[3]]]), n_sites = as.numeric(x[[ii[4]]]))
  } else if (any(grepl("^# ERR", l))) {                                                     # bcftools 1.4.1
    e <- tstrsplit(grep("^ERR\t", l, value = TRUE), "\t", fixed = TRUE)
    x <- data.table(s1 = e[[4]], s2 = e[[5]], n_sites = as.numeric(e[[3]]))[, count := round(as.numeric(e[[2]]) * n_sites)]
  } else { note("gtcheck ", sc, ": unknown output format - raw kept in ", raw); return(NULL) }
  x[s1 != s2][, value := count / n_sites]
}
king_run <- function(sc, d_dir, tg) {
  pre <- file.path(d_dir, "king"); rng <- paste0(pre, ".range.txt")
  t <- fread(tg, header = FALSE, col.names = c("CHROM", "POS"), colClasses = list(character = 1))
  fwrite(t[, .(CHROM, POS, POS, paste0("s", .I))], rng, sep = "\t", col.names = FALSE)
  run(paste(shQuote(PLINK2), "--bcf", shQuote(MASKED), "--allow-extra-chr --chr-set 19 no-xy",
            "--extract bed1", shQuote(rng), "--set-all-var-ids @:#:\\$r:\\$a --make-king-table --out", shQuote(pre),
            ">", shQuote(paste0(pre, ".stdout")), "2>&1"))
  k <- fread(paste0(pre, ".kin0")); setnames(k, sub("^#", "", names(k)))
  k[, .(s1 = IID1, s2 = IID2, n_sites = NSNP, value = KINSHIP, count = IBS0, ibs0_rate = IBS0 / NSNP)]
}
# tumour x normal rows from a pairwise table (either sample order)
to_tn <- function(p, sc, meth) {
  if (is.null(p) || !nrow(p)) return(NULL)
  p <- p[(s1 %chin% TUMOUR_IDS) != (s2 %chin% TUMOUR_IDS) & s1 %chin% SAMPLES & s2 %chin% SAMPLES]
  p[, `:=`(tumour = fifelse(s1 %chin% TUMOUR_IDS, s1, s2), normal = fifelse(s1 %chin% TUMOUR_IDS, s2, s1))]
  p <- unique(p, by = c("tumour", "normal"))
  p[, `:=`(scenario = sc, method = meth)]
  if (meth == "somalier") p[, `:=`(sort1 = value, sort2 = -ibs2_rate)]
  if (meth == "gtcheck")  p[, `:=`(sort1 = value, sort2 = 0)]
  if (meth == "KING")     p[, `:=`(sort1 = -value, sort2 = ibs0_rate)]
  p[, .(scenario, method, tumour, normal, n_sites, value, count, sort1, sort2)]
}

# ---- 4. Run every scenario -------------------------------------------------------
res <- list(); site_tab <- list(); paired <- list()
for (sc in names(EXCLUDE)) {
  ex <- EXCLUDE[[sc]]
  d <- G[!chr_n %chin% ex]
  if (length(ex) && nrow(d) == nrow(G)) note(sc, ": no primary site on chromosome(s) ", paste(ex, collapse = ","), " - nothing removed")
  site_tab[[sc]] <- data.table(scenario = sc, excluded = if (length(ex)) paste(ex, collapse = ",") else "none",
                              n_sites = nrow(d), removed = nrow(G) - nrow(d), n_chromosomes = uniqueN(d$chr_n))
  d_dir <- file.path(OUT, sc); dir.create(d_dir, showWarnings = FALSE)
  tg <- file.path(d_dir, "targets.tsv"); fwrite(d[, .(CHROM, POS)], tg, sep = "\t", col.names = FALSE)
  res[[paste(sc, "custom")]] <- custom_run(d, sc)
  if (HAVE_SOM)  res[[paste(sc, "somalier")]] <- to_tn(somalier_run(sc, d_dir, ex), sc, "somalier")
  if (length(GT_ARGS)) res[[paste(sc, "gtcheck")]] <- to_tn(gtcheck_run(sc, d_dir, tg), sc, "gtcheck")
  if (HAVE_KING) res[[paste(sc, "KING")]] <- to_tn(king_run(sc, d_dir, tg), sc, "KING")
  # paired sign test (custom), best vs second normal
  cu <- res[[paste(sc, "custom")]][order(tumour, sort1)]
  paired[[sc]] <- cu[, { pp <- paired_p(d, tumour[1], normal[1], normal[2])
    .(scenario = sc, best_normal = normal[1], second_normal = normal[2], disc_best_only = pp[["disc_best_only"]],
      disc_second_only = pp[["disc_second_only"]], paired_p = pp[["p"]]) }, by = tumour]
}
tn <- rbindlist(res, fill = TRUE)
site_tab <- rbindlist(site_tab); paired <- rbindlist(paired)
print(site_tab)

# ---- 5. Best normal per scenario x method x tumour ----------------------------------------
setorder(tn, scenario, method, tumour, sort1, sort2)
tn[, rank := seq_len(.N), by = .(scenario, method, tumour)]
best <- tn[, .(best_normal = normal[1], best_value = value[1], second_normal = normal[2], second_value = value[2],
               n_sites_best = n_sites[1]), by = .(scenario, method, tumour)]
best[, `:=`(labelled_normal = KI_NORMALS[tumour], hypothesis_normal = HYPOTHESIS[tumour])]
best[, `:=`(matches_hypothesis = best_normal == hypothesis_normal, matches_label = best_normal == labelled_normal)]
best <- paired[, .(scenario, tumour, paired_p_custom = paired_p)][best, on = .(scenario, tumour)]
best[method != "custom", paired_p_custom := NA]
best[, ambiguous := method == "custom" & paired_p_custom >= PAIRED_ALPHA]
best[, scenario := factor(scenario, names(EXCLUDE))]
best[, method := factor(method, c("custom", "somalier", "gtcheck", "KING"))]
setorder(best, scenario, method, tumour)
print(best[, .(scenario, method, tumour, best_normal, best_value = signif(best_value, 3), second_normal,
               second_value = signif(second_value, 3), n_sites_best, paired_p_custom = signif(paired_p_custom, 3), matches_hypothesis)])
if (any(best$matches_hypothesis == FALSE)) note("best normal differs from the swap hypothesis in: ",
  paste(best[matches_hypothesis == FALSE, paste(scenario, method, tumour)], collapse = "; "))

# Compare "all" with mouse_wgs_somalier.R verdict.tsv (sanity check of this script)
vf <- file.path(SOM_OUT, "verdict.tsv")
if (file.exists(vf)) {
  v <- fread(vf)[analysis == "primary", .(tumour, verdict_best = best_normal, verdict_rate = best_disc_rate)]
  chk <- v[best[scenario == "all" & method == "custom", .(tumour, best_normal, best_value)], on = "tumour"]
  chk[, same := verdict_best == best_normal & abs(verdict_rate - best_value) < 1e-9]
  print(chk)
  if (!all(chk$same)) note("scenario 'all' does not reproduce verdict.tsv - check that the inputs come from the same run")
}
# The jackknife row of mouse_wgs_somalier.R for each dropped chromosome
jf <- file.path(SOM_OUT, "jackknife.tsv")
jk <- if (file.exists(jf)) fread(jf, colClasses = list(character = "dropped_chr"))[dropped_chr %chin% unlist(EXCLUDE)] else data.table()
if (nrow(jk)) print(jk)

# ---- 6. Figure: 3 x 3 per method, with and without chr8 ------------------------------------
lab_fmt <- function(m, v) ifelse(m == "custom" | m == "gtcheck" | m == "somalier", sprintf("%.1f%%", 100 * v), sprintf("%.2f", v))
tn[, label := lab_fmt(method, value)]
tn[, role := fcase(normal == HYPOTHESIS[tumour] & normal == KI_NORMALS[tumour], "labelled = hypothesis",
                   normal == HYPOTHESIS[tumour], "hypothesis (swap)", normal == KI_NORMALS[tumour], "labelled", default = "")]
tn[, method_lab := factor(method, c("custom", "somalier", "gtcheck", "KING"),
                          c("Discordance at normal-hom SNVs\n(custom, primary)", "IBS0 rate\n(somalier)",
                            "Genotype discordance\n(bcftools gtcheck)", "Kinship\n(KING-robust)"))]
tn[, scenario_lab := factor(scenario, names(EXCLUDE), paste0(names(EXCLUDE), " (", site_tab$n_sites, " sites)"))]
p <- ggplot(tn, aes(normal, factor(tumour, rev(TUMOUR_IDS)))) +
  geom_tile(aes(fill = factor(rank)), colour = "white", linewidth = 1) +
  geom_tile(data = tn[role != ""], aes(linetype = role), fill = NA, colour = "black", linewidth = 0.7,
            width = 0.92, height = 0.92) +
  geom_text(aes(label = label, colour = rank == 1), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  facet_grid(scenario_lab ~ method_lab) +
  scale_fill_manual(values = c(`1` = "#1757ac", `2` = "#9ec2f0", `3` = "#eef5fe"), labels = c("best", "2nd", "3rd")) +
  scale_linetype_manual(values = c(`hypothesis (swap)` = "solid", labelled = "dashed", `labelled = hypothesis` = "dotted")) +
  labs(x = "KI normal", y = "tumour", fill = "match rank\n(per tumour)", linetype = NULL,
       title = "Sample matching with and without chromosome 8 (Sprtn knock-in)",
       subtitle = paste0("Same masked FreeBayes genotypes, primary sites (", COMBO, "). Rank 1 = best-matching normal per tumour and method.")) +
  theme_bw(base_size = 10) + theme(plot.title.position = "plot", panel.grid = element_blank())
ggsave(file.path(OUT, "exclude_chr_matching.png"), p, width = 13, height = 2 + 2.4 * length(EXCLUDE), dpi = 300)
ggsave(file.path(OUT, "exclude_chr_matching.pdf"), p, width = 13, height = 2 + 2.4 * length(EXCLUDE))

# ---- 7. Outputs, RUN_LOG.md ----------------------------------------------------
fwrite(tn[, .(scenario, method, tumour, normal, n_sites, value, count, rank)], file.path(OUT, "tumour_normal_matrix.tsv"), sep = "\t")
fwrite(best, file.path(OUT, "best_normal.tsv"), sep = "\t")
fwrite(site_tab, file.path(OUT, "scenario_sites.tsv"), sep = "\t")
fwrite(paired, file.path(OUT, "paired_tests.tsv"), sep = "\t")
git_hash <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)),
                     error = function(e) "unknown")
if (!length(git_hash)) git_hash <- "unknown"
md_table <- function(x) { x <- as.data.table(x); if (!nrow(x)) return("(none)")
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))) }
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_identity_exclude_chr"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Tools:      ", paste(names(versions), versions, sep = ": ", collapse = "; ")),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Inputs:     ", GT_TSV, "; ", MASKED, "; ", SITES),
  paste0("- Scenarios:  ", paste(names(EXCLUDE), sapply(EXCLUDE, function(e) if (length(e)) paste0("minus chr ", paste(e, collapse = ",")) else "all sites"), sep = " = ", collapse = "; ")),
  paste0("- Settings:   somalier relate -d ", MIN_DP, " --min-ab ", HET_MIN, " SOMALIER_AB_HOM_CUTOFF=", AB_HOM_CUTOFF,
         "; gtcheck ", paste(GT_ARGS, collapse = " "), "; KING plink2 --make-king-table"),
  paste0("- Outputs:    ", OUT),
  "", "### Sites per scenario", "", md_table(site_tab),
  "", "### Best normal per tumour", "",
  md_table(best[, .(scenario, method, tumour, best_normal, best_value = signif(best_value, 3), second_normal,
                    second_value = signif(second_value, 3), n_sites_best, paired_p_custom = signif(paired_p_custom, 3), matches_hypothesis)]),
  if (nrow(jk)) c("", "### mouse_wgs_somalier.R jackknife rows for the dropped chromosome", "", md_table(jk)) else character(),
  "", "### Notes",
  "- Same masked FreeBayes genotypes as mouse_wgs_somalier.R: tests whether chr8 carries the signal, not independent evidence.",
  "- Scenario 'all' must reproduce verdict.tsv of mouse_wgs_somalier.R (checked in the script).",
  "- KING kinship is unstable with few hets (inbred mice): use the ranking.",
  if (length(NOTES)) paste0("- NOTE: ", NOTES) else "- no notes",
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
