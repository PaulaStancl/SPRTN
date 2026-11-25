#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_identity_strain_check.R - are the identity sites (and the "somatic-looking" differences between
# a tumour and its labelled normal) known germline polymorphisms, and do they behave like inherited blocks?
#
# For every panel A site of mouse_wgs_somalier.R (masked 6-sample genotypes):
#   A. Mouse Genomes Project (Sanger MGP, GRCm38 release, e.g. REL-1505): in MGP with the same REF/ALT, ID
#      (rsID where the release has one), number of strains carrying the ALT and which, C57BL/6NJ and the 129
#      strains (possible donor-strain blocks of the knock-in line). C57BL/6J is the reference genome, so its
#      genotype is 0/0 by definition (not a sample in the MGP VCF).
#   B. WT normals (WT1-3) from the setting-B Mutect2 joint files: genotype from AD (same VAF bins), where the
#      file has a record. Mutect2 only writes sites with tumour ALT evidence, so many sites have no record.
#   C. Neighbouring SNPs: for each site where a tumour differs from its labelled normal ("discordant", i.e.
#      it would be called somatic in setting A), the other panel A sites within +-NEIGHBOUR_WINDOW. Inherited
#      differences between two mice come in blocks: the neighbours are discordant too and the tumour matches
#      another mouse's normal there. A real somatic mutation is isolated: its neighbours are concordant.
#
# Summaries per tumour: discordant sites vs the labelled normal and vs the genetic match (HYPOTHESIS), the
# share that are known MGP polymorphisms, carried by a 129 strain, carried by the WT normals, and the share
# with discordant neighbours. Isolated discordant sites are listed (candidates for real somatic mutations
# or genotype errors).
# Not in a database does NOT mean somatic: private germline variants exist. MGI's SNP query and GenomeMUSter
# are manual (web) checks; for GRCm38 they largely draw on the same MGP / dbSNP calls.
#
# Env on ssi: r-mutation (data.table, ggplot2, bcftools with libcurl for a remote MGP file).
# Output: results/mouse_wgs/mouse_wgs_identity_strain_check/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
SOM      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_somalier")      # outputs of mouse_wgs_somalier.R
INF      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_identity_informativeness")
DIR_B    <- file.path(REPO, "data/mouse/wgs/tumor_vs_joint_wt")        # setting B: Mutect2 joint, WT normals
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_identity_strain_check")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_identity_strain_check.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")

# MGP SNP VCF on GRCm38 (bgzipped + .tbi next to it). A local path, or a URL that bcftools reads directly
# (only the needed blocks are fetched). Must be GRCm38: REL-1505 (default; ID = dbSNP142 rsIDs) or REL-1807
# (mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz), NOT REL-2112 (GRCm39). EBI has no REL-2005 folder
# (checked 2026-11-25). NA = skip the MGP part.
MGP_VCF <- "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"

TUMOURS    <- c("HCC_1", "HCC_2", "HCC_3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing
HYPOTHESIS <- c(HCC_1 = "KI2", HCC_2 = "KI1", HCC_3 = "KI3")   # genetic match from mouse_wgs_somalier.R

NEIGHBOUR_WINDOW <- 5e6      # bp either side (panel A has ~1 site per Mb; inherited blocks are usually several Mb)
MIN_NEIGHBOURS   <- 3        # an "isolated" discordant site needs at least this many genotyped neighbours
WT_MIN_DP        <- 20       # REF+ALT reads for a WT genotype from the Mutect2 AD
HOM_REF_MAX <- 0.05; HET_MIN <- 0.35; HET_MAX <- 0.65; HOM_ALT_MIN <- 0.95   # as mouse_wgs_somalier.R

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
NORMALS <- unname(KI_NORMALS)
SAMPLES <- c(TUMOURS, NORMALS)
GT <- c("0/0", "0/1", "1/1")
B <- shQuote(BCFTOOLS)

bcf <- function(args) {
  err <- tempfile()
  out <- suppressWarnings(system2(BCFTOOLS, args, stdout = TRUE, stderr = err))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0) stop("bcftools failed: bcftools ", paste(args, collapse = " "), "\n", paste(readLines(err, warn = FALSE), collapse = "\n"))
  out
}
vaf_gt <- function(ref, alt, min_dp) {
  dp <- ref + alt; v <- alt / dp
  fifelse(is.na(dp) | dp < min_dp, NA_integer_,
          fcase(v <= HOM_REF_MAX, 0L, v >= HET_MIN & v <= HET_MAX, 1L, v >= HOM_ALT_MIN, 2L, default = NA_integer_))
}
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150, limitsize = FALSE)
}
pct <- function(x, n) fifelse(n > 0, round(100 * x / n, 1), NA_real_)

# ---- 2. Panel A sites (masked 6-sample genotypes) ----------------------------
MASKED <- file.path(SOM, "masked", "masked_6samples.bcf")
if (!file.exists(MASKED)) MASKED <- file.path(SOM, "masked", "masked_6samples.vcf")
if (!file.exists(MASKED)) stop("missing input from mouse_wgs_somalier.R: ", MASKED)
smp <- bcf(c("query", "-l", shQuote(MASKED)))
stopifnot(setequal(smp, SAMPLES))
cols <- c("CHROM", "POS", "REF", "ALT", paste0("GT_", smp))
A <- fread(cmd = paste(B, "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"), shQuote(MASKED)), sep = "\t", header = FALSE,
           quote = "", col.names = cols, colClasses = c("character", "integer", rep("character", length(cols) - 2)))
for (s in SAMPLES) A[, (s) := match(get(paste0("GT_", s)), GT) - 1L]
A[, paste0("GT_", smp) := NULL]
A[, key := paste(CHROM, POS, REF, ALT, sep = ":")]
A[, chr_n := suppressWarnings(as.integer(sub("^chr", "", CHROM)))]
setorder(A, chr_n, POS)
A[, informative := apply(.SD, 1, uniqueN) > 1, .SDcols = NORMALS]
message("panel A: ", nrow(A), " sites")

# ---- 3. Mouse Genomes Project ------------------------------------------------
mgp <- data.table(); MGP_STATUS <- "skipped (MGP_VCF = NA)"; MGP_STRAINS <- character()
if (!is.na(MGP_VCF)) {
  hdr <- tryCatch(bcf(c("view", "-h", shQuote(MGP_VCF))), error = function(e) { message("MGP not readable: ", conditionMessage(e)); NULL })
  if (!is.null(hdr)) {
    MGP_STRAINS <- strsplit(sub("^#CHROM.*FORMAT\t", "", grep("^#CHROM", hdr, value = TRUE)), "\t", fixed = TRUE)[[1]]
    asm <- grep("^##(reference|assembly|contig=<ID=1,)", hdr, value = TRUE)
    message("MGP: ", length(MGP_STRAINS), " strains; header lines on the assembly:\n  ", paste(head(asm, 3), collapse = "\n  "))
    if (any(grepl("GRCm39|mm39", asm))) stop("MGP_VCF looks like GRCm39 - use a GRCm38 release (REL-2005 / REL-1505)")
    mgp_chr <- any(grepl("^##contig=<ID=chr", hdr))
    reg <- file.path(OUT, "regions_panelA.tsv")
    fwrite(A[, .(CHROM = paste0(if (mgp_chr) "chr" else "", sub("^chr", "", CHROM)), POS, POS2 = POS)], reg, sep = "\t", col.names = FALSE)
    cmd <- paste(B, "view -R", shQuote(reg), "-Ou", shQuote(MGP_VCF), "|", B, "norm -m -any -Ou |", B,
                 "query -f", shQuote("%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n"))
    mgp <- fread(cmd = cmd, sep = "\t", header = FALSE, quote = "", col.names = c("CHROM", "POS", "ID", "REF", "ALT", MGP_STRAINS),
                 colClasses = c("character", "integer", rep("character", 3 + length(MGP_STRAINS))))
    mgp[, CHROM := sub("^chr", "", CHROM)]
    MGP_STATUS <- sprintf("%s (%d strains, %d records at panel A positions)", MGP_VCF, length(MGP_STRAINS), nrow(mgp))
  } else MGP_STATUS <- paste("NOT READ:", MGP_VCF)
}
message("MGP: ", MGP_STATUS)

if (nrow(mgp)) {
  carry <- as.matrix(mgp[, lapply(.SD, function(g) grepl("1", g, fixed = TRUE)), .SDcols = MGP_STRAINS])
  called <- as.matrix(mgp[, lapply(.SD, function(g) !grepl(".", g, fixed = TRUE)), .SDcols = MGP_STRAINS])
  s129 <- grep("^129", MGP_STRAINS, value = TRUE); sB6 <- grep("^C57BL", MGP_STRAINS, value = TRUE)
  mgp[, `:=`(mgp_n_strains_alt = rowSums(carry), mgp_n_strains_called = rowSums(called),
             mgp_strains_alt = apply(carry, 1, function(v) paste(MGP_STRAINS[v], collapse = ",")),
             mgp_129_alt = if (length(s129)) rowSums(carry[, s129, drop = FALSE]) > 0 else NA,
             mgp_C57BL_alt = if (length(sB6)) rowSums(carry[, sB6, drop = FALSE]) > 0 else NA)]
  keep_gt <- intersect(c(sB6, s129), MGP_STRAINS)
  m <- mgp[, c("CHROM", "POS", "REF", "ALT", "ID", "mgp_n_strains_alt", "mgp_n_strains_called", "mgp_strains_alt",
               "mgp_129_alt", "mgp_C57BL_alt", keep_gt), with = FALSE]
  setnames(m, c("ID", keep_gt), c("mgp_ID", paste0("mgp_", keep_gt)))
  m <- unique(m, by = c("CHROM", "POS", "REF", "ALT"))
  A <- m[A, on = .(CHROM, POS, REF, ALT)]
  A[, mgp_same_alt := !is.na(mgp_n_strains_alt)]
  A[, mgp_position := paste(sub("^chr", "", CHROM), POS) %chin% mgp[, paste(CHROM, POS)]]   # any MGP record at the position
  A[, mgp_known_polymorphism := mgp_same_alt & mgp_n_strains_alt > 0]
} else A[, c("mgp_same_alt", "mgp_position", "mgp_known_polymorphism", "mgp_129_alt") := NA]

# ---- 4. WT normals from the setting-B Mutect2 joint files ----------------------
wt <- data.table(); WT_IDS <- character()
fB <- if (dir.exists(DIR_B)) list.files(DIR_B, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf\\.gz)$") else character()
fB <- fB[grepl("mutect2", basename(fB), ignore.case = TRUE) & !grepl("union|merged", fB, ignore.case = TRUE)]
if (length(fB)) {
  tgt <- file.path(OUT, "targets_panelA.tsv")
  fwrite(A[, .(CHROM, POS)], tgt, sep = "\t", col.names = FALSE)
  wt <- rbindlist(lapply(fB, function(f) {
    x <- system(paste(B, "view -T", shQuote(tgt), "-Ou", shQuote(f), "|", B, "norm -m -any -Ou |", B, "query -f",
                      shQuote("%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%AD]\\n")), intern = TRUE)
    if (!length(x)) return(NULL)
    d <- fread(text = x, sep = "\t", header = FALSE, fill = TRUE, colClasses = "character")
    setnames(d, 1:4, c("CHROM", "POS", "REF", "ALT"))
    d <- melt(d, id.vars = 1:4, value.name = "sa")[!is.na(sa) & sa != ""]
    d[, `:=`(sample = sub("=.*$", "", sa), ad = sub("^[^=]*=", "", sa))]
    d[, `:=`(ref = suppressWarnings(as.integer(sub(",.*$", "", ad))),
             alt = fifelse(grepl(",", ad, fixed = TRUE), suppressWarnings(as.integer(sub("^[^,]*,([^,]*).*$", "\\1", ad))), NA_integer_))]
    d[, .(CHROM, POS = as.integer(POS), REF, ALT, sample, ref, alt, file = basename(f))]
  }))
  if (nrow(wt)) {
    wt <- wt[!grepl("KI_HCC_[0-9]|KI[0-9]$", sample)]                  # keep WT normals only
    wt[, short := sub("^(id_[0-9]+_)?(sfb_dikic_)?([0-9]{4}_[0-9]{2}_[0-9]{2}_)?", "", sample)]
    wt <- wt[order(-(ref + alt))][!duplicated(wt[, .(CHROM, POS, REF, ALT, short)])]   # deepest record per WT and site
    wt[, g := vaf_gt(ref, alt, WT_MIN_DP)]
    WT_IDS <- sort(unique(wt$short))
    ww <- dcast(wt, CHROM + POS + REF + ALT ~ short, value.var = "g")
    setnames(ww, WT_IDS, paste0("WT_", WT_IDS))
    A <- ww[A, on = .(CHROM, POS, REF, ALT)]
    A[, wt_n_genotyped := rowSums(!is.na(.SD)), .SDcols = paste0("WT_", WT_IDS)]
    A[, wt_any_alt := rowSums(.SD > 0, na.rm = TRUE) > 0, .SDcols = paste0("WT_", WT_IDS)]
    A[wt_n_genotyped == 0, wt_any_alt := NA]
  }
}
if (!length(WT_IDS)) { message("no WT genotypes (setting-B Mutect2 files not found or no records)"); A[, c("wt_n_genotyped", "wt_any_alt") := NA] }
message("WT normals: ", if (length(WT_IDS)) paste(WT_IDS, collapse = ", ") else "none")

# ---- 5. Discordant sites per tumour and their neighbours ---------------------
long <- rbindlist(lapply(TUMOURS, function(t) {
  data.table(tumour = t, site = A$key, CHROM = A$CHROM, POS = A$POS, tg = A[[t]],
             lab = A[[KI_NORMALS[[t]]]], gen = A[[HYPOTHESIS[[t]]]], KI1 = A$KI1, KI2 = A$KI2, KI3 = A$KI3)
}))
long <- long[!is.na(tg)]
long[, `:=`(disc_labelled = tg != lab, disc_genetic = tg != gen)]
# neighbours within the window (same tumour, same chromosome, other site)
nb <- long[long, on = .(tumour, CHROM), allow.cartesian = TRUE, nomatch = NULL][
  site != i.site & abs(POS - i.POS) <= NEIGHBOUR_WINDOW,
  .(n_neighbours = .N, nb_disc_labelled = sum(disc_labelled), nb_match_genetic = sum(!disc_genetic)),
  by = .(tumour, site = i.site)]
long <- nb[long, on = .(tumour, site)]
for (cl in c("n_neighbours", "nb_disc_labelled", "nb_match_genetic")) set(long, which(is.na(long[[cl]])), cl, 0L)
long[, isolated := disc_labelled & n_neighbours >= MIN_NEIGHBOURS & nb_disc_labelled == 0L]
# Candidate = isolated AND the tumour genotype is in none of the KI normals: not explained by another
# mouse, so somatic, LOH (het normal -> hom tumour) or a genotype error. Isolated sites where the tumour
# matches another mouse's normal are inherited differences that just have no differing neighbour.
long[, in_no_normal := tg != KI1 & tg != KI2 & tg != KI3]
long[, candidate := isolated & in_no_normal]
long[, candidate_type := fifelse(!candidate, NA_character_,
                                 fifelse(lab == 1L & tg != 1L, "LOH-compatible (het normal, hom tumour)", "tumour gained an allele"))]
ann_cols <- intersect(c("informative", "mgp_same_alt", "mgp_known_polymorphism", "mgp_129_alt", "mgp_C57BL_alt",
                        "mgp_n_strains_alt", "mgp_ID", "wt_any_alt", "wt_n_genotyped"), names(A))
long <- setnames(A[, c("key", ann_cols), with = FALSE], "key", "site")[long, on = "site"]

# 5a. summary per tumour: discordant vs labelled normal (= somatic-looking in setting A) and vs genetic match
summ <- rbindlist(lapply(c("labelled", "genetic"), function(w) {
  dc <- paste0("disc_", w)
  long[, {
    d <- .SD[get(dc) == TRUE]
    .(compared_with = if (w == "labelled") KI_NORMALS[[.BY$tumour]] else HYPOTHESIS[[.BY$tumour]],
      comparison = if (w == "labelled") "labelled normal" else "genetic match",
      sites_genotyped = .N, discordant = nrow(d), pct_discordant = pct(nrow(d), .N),
      disc_in_MGP_same_alt = sum(d$mgp_same_alt, na.rm = TRUE),
      pct_disc_known_polymorphism = pct(sum(d$mgp_known_polymorphism, na.rm = TRUE), nrow(d)),
      pct_disc_129_carries_alt = pct(sum(d$mgp_129_alt, na.rm = TRUE), nrow(d)),
      disc_with_WT_genotype = sum(d$wt_n_genotyped > 0, na.rm = TRUE),
      pct_of_those_WT_carries_alt = pct(sum(d$wt_any_alt, na.rm = TRUE), sum(d$wt_n_genotyped > 0, na.rm = TRUE)),
      pct_disc_with_discordant_neighbour = pct(sum(d$nb_disc_labelled > 0), nrow(d)),
      isolated_discordant = sum(d$isolated), isolated_matching_another_normal = sum(d$isolated & !d$in_no_normal),
      candidates_in_no_normal = sum(d$candidate))
  }, by = tumour]
}))
print(summ)

# 5b. all panel A sites: how many are known polymorphisms (background for 5a)
bg <- A[, .(sites = .N, informative = sum(informative),
            pct_MGP_same_alt = pct(sum(mgp_same_alt, na.rm = TRUE), .N),
            pct_known_polymorphism = pct(sum(mgp_known_polymorphism, na.rm = TRUE), .N),
            pct_129_carries_alt = pct(sum(mgp_129_alt, na.rm = TRUE), .N),
            pct_informative_known = pct(sum(mgp_known_polymorphism & informative, na.rm = TRUE), sum(informative)))]
print(bg)

# 5c. isolated discordant sites; candidate = tumour genotype in none of the KI normals
isolated <- long[isolated == TRUE][order(-candidate, tumour, CHROM, POS)]
print(isolated[, .(tumour, site, tg = GT[tg + 1], labelled = GT[lab + 1], genetic = GT[gen + 1], n_neighbours,
                   matches_another_normal = !in_no_normal, candidate_type,
                   mgp_known_polymorphism = if ("mgp_known_polymorphism" %in% names(isolated)) mgp_known_polymorphism else NA)])

# ---- 6. Plots ------------------------------------------------------------------
pl <- long[disc_labelled == TRUE]
if (nrow(pl)) {
  pl[, chr_n := as.integer(sub("^chr", "", CHROM))]
  pl[, status := fcase(candidate, "isolated, tumour genotype in no KI normal",
                       isolated, "isolated, matches another mouse's normal",
                       mgp_known_polymorphism %in% TRUE, "known MGP polymorphism",
                       default = if (nrow(mgp)) "not in MGP (same ALT)" else "MGP not checked")]
  p_gen <- ggplot(pl, aes(POS / 1e6, factor(chr_n, levels = 19:1), colour = status)) +
    geom_point(shape = 124, size = 4) + facet_wrap(~ tumour, nrow = 1) +
    scale_colour_manual(values = c("known MGP polymorphism" = "#2E7D32", "not in MGP (same ALT)" = "#F9A825",
                                   "isolated, tumour genotype in no KI normal" = "#C62828",
                                   "isolated, matches another mouse's normal" = "#7B1FA2", "MGP not checked" = "grey50"), name = NULL) +
    labs(x = "position (Mb)", y = "chromosome",
         title = "Sites where the tumour differs from its LABELLED normal (would be called somatic in setting A)",
         subtitle = sprintf("Blocks = inherited differences between two mice; isolated = no discordant site within %g Mb", NEIGHBOUR_WINDOW / 1e6)) +
    theme_bw() + theme(legend.position = "bottom")
  save_plot(p_gen, "discordant_sites_genome", 14, 6)
}
sb <- melt(summ, id.vars = c("tumour", "comparison"), measure.vars = c("pct_disc_known_polymorphism", "pct_disc_129_carries_alt",
                                                                    "pct_disc_with_discordant_neighbour"))
p_sum <- ggplot(sb[comparison == "labelled normal"], aes(variable, value, fill = tumour)) + geom_col(position = "dodge") +
  geom_text(aes(label = value), position = position_dodge(0.9), vjust = -0.3, size = 3) +
  scale_x_discrete(labels = c("known MGP\npolymorphism", "a 129 strain\ncarries the ALT", "has a discordant\nneighbour")) +
  labs(x = NULL, y = "% of discordant sites (tumour vs labelled normal)", fill = NULL,
       title = "Are the tumour vs labelled-normal differences germline?") + theme_bw()
save_plot(p_sum, "discordant_sites_germline_evidence", 8, 5)

# ---- 7. Outputs, RUN_LOG.md --------------------------------------------------
fwrite(A, file.path(OUT, "panelA_sites_annotated.tsv"), sep = "\t")
fwrite(long, file.path(OUT, "tumour_sites_discordance_neighbours.tsv"), sep = "\t")
fwrite(summ, file.path(OUT, "summary_discordant_sites.tsv"), sep = "\t")
fwrite(bg, file.path(OUT, "summary_all_sites.tsv"), sep = "\t")
fwrite(isolated, file.path(OUT, "isolated_discordant_sites.tsv"), sep = "\t")
bam_csv <- file.path(INF, "bam_check_sites_for_patrik.csv")
if (file.exists(bam_csv)) {
  bc <- fread(bam_csv, colClasses = list(character = c("CHROM", "REF", "ALT")))
  ac <- intersect(c("mgp_ID", "mgp_n_strains_alt", "mgp_strains_alt", "mgp_129_alt", "mgp_C57BL_alt",
                    grep("^mgp_(129|C57BL)", names(A), value = TRUE), paste0("WT_", WT_IDS)), names(A))
  bc <- A[, c("CHROM", "POS", "REF", "ALT", ac), with = FALSE][bc, on = .(CHROM, POS, REF, ALT)]
  for (w in intersect(paste0("WT_", WT_IDS), names(bc))) set(bc, j = w, value = fifelse(is.na(bc[[w]]), "no record/masked", GT[bc[[w]] + 1L]))
  setcolorder(bc, names(fread(bam_csv, nrows = 0)))
  fwrite(bc, file.path(OUT, "bam_check_sites_for_patrik_annotated.csv"))
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_identity_strain_check"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Command:    Rscript ", SCRIPT),
  paste0("- Inputs:     ", MASKED, " (md5 ", tools::md5sum(MASKED), "); setting-B Mutect2 files: ", length(fB),
         "; MGP: ", MGP_STATUS),
  paste0("- Parameters: NEIGHBOUR_WINDOW ", NEIGHBOUR_WINDOW, ", MIN_NEIGHBOURS ", MIN_NEIGHBOURS, ", WT_MIN_DP ", WT_MIN_DP,
         ", VAF bins ", HOM_REF_MAX, "/", HET_MIN, "-", HET_MAX, "/", HOM_ALT_MIN,
         "; labelled ", paste(names(KI_NORMALS), KI_NORMALS, sep = "-", collapse = ", "),
         "; genetic match ", paste(names(HYPOTHESIS), HYPOTHESIS, sep = "-", collapse = ", ")),
  paste0("- Tools:      ", bcf("--version")[1], "; data.table ", packageVersion("data.table")),
  paste0("- Outputs:    ", OUT),
  paste0("- Counts:     panel A ", nrow(A), " sites; WT normals ", length(WT_IDS), "; isolated discordant ", nrow(isolated), " (candidates in no normal: ", sum(isolated$candidate), ")"),
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
