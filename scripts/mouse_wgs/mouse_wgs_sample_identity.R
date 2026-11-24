#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24
#
# mouse_wgs_sample_identity.R - are the tumour / matched KI normal pairs from the same mouse?
#
# Trigger: chr10:3184122 G>T - HCC_1 hom-alt vs its "matched" KI1 hom-ref, HCC_2 hom-ref vs KI2 hom-alt.
# Principle (genotype fingerprinting, as in Conpair, NGSCheckMate, somalier, Picard CrosscheckFingerprints,
# bcftools gtcheck): samples of one animal share their inherited genotypes; a confidently hom-ref vs
# hom-alt site (IBS0, "opposite homozygotes") is not expected between them, except through rarer events
# (genotype misclassification, mapping artefacts, somatic mutation + LOH/amplification, allele-specific
# coverage). Many IBS0 sites spread over many independent genomic regions indicate different animals.
#
# Terminology: VAF gives an INFERRED GENOTYPE under a diploid assumption (not copy number):
#   ref  = VAF <= HOM_REF_MAX      het = HET_MIN <= VAF <= HET_MAX      alt = VAF >= HOM_ALT_MIN
#   confident = one of these and REF+ALT reads >= MIN_DP; anything else is not used.
#
# Data (no BAMs exist, so everything comes from the VCFs):
#   FreeBayes, setting A (tumor_vs_normal/snv_indel): tumours HCC_1-3 and KI normals KI1-3 at every site
#     FreeBayes found variable in that pair (germline sites included, all FILTER values, QUAL >= MIN_QUAL)
#   Mutect2 joint, setting B (tumor_vs_joint_wt): the WT normals, at Mutect2 candidate sites only
#     (ascertained by the tumour; WT comparisons rest on fewer, differently chosen sites)
#   A pair is compared only at sites present - and confidently genotyped - in both samples' files.
#
# Blocks:
#   0  EDIT THIS BLOCK
#   1  setup
#   2  files
#   3  read genotypes (FreeBayes A; Mutect2 B for the WT normals)
#   4  inferred genotypes; site x sample tables
#   5  pairwise statistics: shared confident sites, IBS0 count and rate, concordance, LD-thinned rate,
#      independent IBS0 blocks and chromosomes
#   6  verdict per tumour
#   7  plots: IBS0-rate heatmap (all samples), per-tumour bars, VAF scatter, genome distribution of IBS0,
#      depth, genotype/mapping quality, raw vs LD-thinned rate
#   8  optional: IBS0 sites vs Mouse Genomes Project strains (MGP_VCF)
#   9  bcftools gtcheck on the merged FreeBayes genotypes (established tool, same question)
#   10 outputs, summary.md, RUN_LOG.md
#
# Not possible with the available data (listed in summary.md): BAM-level pileup, restriction to copy-neutral
# regions (no copy-number calls), strategies D (GM normal vs WT) and E (WT vs WT) as new caller runs.
#
# Env on ssi: r-mutation (data.table, ggplot2, bcftools).
# Output: results/mouse_wgs/mouse_wgs_sample_identity/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
DIR_A    <- file.path(REPO, "data/mouse/wgs/tumor_vs_normal/snv_indel")
DIR_B    <- file.path(REPO, "data/mouse/wgs/tumor_vs_joint_wt")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_sample_identity")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_sample_identity.R")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")

TUMOURS    <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1", HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
                HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing

MIN_QUAL    <- 30      # FreeBayes site QUAL
MIN_DP      <- 30      # REF+ALT reads for a confident genotype
HOM_REF_MAX <- 0.05
HET_MIN     <- 0.35
HET_MAX     <- 0.65
HOM_ALT_MIN <- 0.95
THIN_BP     <- 1e6     # LD thinning: keep one shared site per 1 Mb window
BLOCK_GAP   <- 2e6     # IBS0 sites closer than 2 Mb on one chromosome = one block

# Optional: Mouse Genomes Project SNP VCF on GRCm38 (tabix-indexed), e.g. REL-1505
# mgp.v5.merged.snps_all.dbSNP142.vcf.gz from the Sanger MGP FTP. NA = skip block 8.
MGP_VCF     <- NA
MGP_STRAINS <- c("129S1_SvImJ", "129P2_OlaHsd", "C57BL_6NJ")

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
B <- shQuote(BCFTOOLS)
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150, limitsize = FALSE)
  invisible(p)
}
short_name <- function(s) fcase(grepl("KI_HCC_[0-9]", s), sub(".*KI_(HCC_[0-9]).*", "\\1", s),
                                default = sub("^(id_[0-9]+_)?(sfb_dikic_)?([0-9]{4}_[0-9]{2}_[0-9]{2}_)?", "", s))
sample_class <- function(s) fcase(grepl("^HCC_[0-9]$", s), "tumour", grepl("^KI[0-9]$", s), "KI normal", default = "WT normal")
AUTOSOMES <- as.character(1:19)

# ---- 2. Files ----------------------------------------------------------------
list_vcf <- function(d) list.files(d, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf|vcf\\.gz)$")
fb_files <- data.table(path = list_vcf(DIR_A))[grepl("freebayes", path, ignore.case = TRUE) & !grepl("union|merged", path, ignore.case = TRUE)]
m2_files <- data.table(path = list_vcf(DIR_B))
add_tumour <- function(f) {
  f[, tumour := NA_character_]
  for (k in names(TUMOURS)) f[grepl(TUMOURS[[k]], basename(path), fixed = TRUE), tumour := k]
  f
}
fb_files <- add_tumour(fb_files); m2_files <- add_tumour(m2_files)
fb_files[, .(tumour, file = basename(path))]
m2_files[, .(tumour, file = basename(path))]
stopifnot(nrow(fb_files) == length(TUMOURS), !anyNA(fb_files$tumour))
# Contig naming ("1" or "chr1") from the first FreeBayes header
CHR_PREFIX <- if (any(grepl("^##contig=<ID=chr1,", system2(BCFTOOLS, c("view", "-h", shQuote(fb_files$path[1])), stdout = TRUE)))) "chr" else ""

# ---- 3. Read genotypes -----------------------------------------------------------
# Multiallelic records are split (bcftools norm -m -any) so AD is "ref,alt" for every record.
read_fb <- function(path) {
  smp <- system2(BCFTOOLS, c("query", "-l", shQuote(path)), stdout = TRUE)
  cmd <- paste(B, "norm -m -any -Ou", shQuote(path), "|", B, "view -v snps -i", shQuote(paste0("QUAL>=", MIN_QUAL)), "-Ou |",
               B, "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/MQM\\t%INFO/MQMR[\\t%GT\\t%AD\\t%GQ]\\n"))
  cn <- c("CHROM", "POS", "REF", "ALT", "QUAL", "MQM", "MQMR", paste0(rep(c("GT_", "AD_", "GQ_"), length(smp)), rep(seq_along(smp), each = 3)))
  x <- fread(cmd = cmd, sep = "\t", header = FALSE, col.names = cn, colClasses = "character")
  x <- x[sub("^chr", "", CHROM) %chin% AUTOSOMES]
  rbindlist(lapply(seq_along(smp), function(j) {
    ad <- tstrsplit(x[[paste0("AD_", j)]], ",", fixed = TRUE)
    data.table(CHROM = sub("^chr", "", x$CHROM), POS = as.numeric(x$POS), REF = x$REF, ALT = x$ALT,
               sample = short_name(smp[j]), source = "freebayes A", GT = x[[paste0("GT_", j)]],
               REF_N = suppressWarnings(as.numeric(ad[[1]])), ALT_N = suppressWarnings(as.numeric(if (length(ad) > 1) ad[[2]] else NA)),
               GQ = suppressWarnings(as.numeric(x[[paste0("GQ_", j)]])), QUAL = suppressWarnings(as.numeric(x$QUAL)),
               MQM = suppressWarnings(as.numeric(x$MQM)), MQMR = suppressWarnings(as.numeric(x$MQMR)))
  }))
}
read_m2_wt <- function(path) {            # WT normals only (tumours come from FreeBayes)
  smp <- system2(BCFTOOLS, c("query", "-l", shQuote(path)), stdout = TRUE)
  keep <- which(sample_class(short_name(smp)) == "WT normal")
  if (!length(keep)) return(NULL)
  cmd <- paste(B, "norm -m -any -Ou", shQuote(path), "|", B, "view -v snps -Ou |",
               B, "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n"))
  cn <- c("CHROM", "POS", "REF", "ALT", paste0("AD_", seq_along(smp)))
  x <- fread(cmd = cmd, sep = "\t", header = FALSE, col.names = cn, colClasses = "character")
  x <- x[sub("^chr", "", CHROM) %chin% AUTOSOMES]
  rbindlist(lapply(keep, function(j) {
    ad <- tstrsplit(x[[paste0("AD_", j)]], ",", fixed = TRUE)
    data.table(CHROM = sub("^chr", "", x$CHROM), POS = as.numeric(x$POS), REF = x$REF, ALT = x$ALT,
               sample = short_name(smp[j]), source = "mutect2 joint B",
               REF_N = suppressWarnings(as.numeric(ad[[1]])), ALT_N = suppressWarnings(as.numeric(if (length(ad) > 1) ad[[2]] else NA)))
  }))
}
calls <- rbindlist(c(lapply(fb_files$path, read_fb), lapply(m2_files$path, read_m2_wt)), fill = TRUE)
calls[, site := paste(CHROM, POS, REF, ALT, sep = ":")]
calls[, DPA := REF_N + ALT_N]                                      # reads supporting REF or ALT
setorder(calls, site, sample, -DPA)
calls <- unique(calls, by = c("site", "sample"))                   # WT appear in all 3 B files: keep deepest
calls[, .(records = .N), by = .(sample, source)][order(sample)]

# ---- 4. Inferred genotypes ---------------------------------------------------
calls[, VAF := fifelse(DPA > 0, ALT_N / DPA, NA_real_)]
calls[, call := fcase(is.na(VAF) | DPA < MIN_DP, NA_character_,
                      VAF <= HOM_REF_MAX, "ref", VAF >= HET_MIN & VAF <= HET_MAX, "het", VAF >= HOM_ALT_MIN, "alt",
                      default = NA_character_)]
calls[, .N, by = .(sample, call)][order(sample, call)]
conf <- calls[!is.na(call)]
SAMPLES <- intersect(c(names(TUMOURS), unname(KI_NORMALS), sort(unique(conf[source == "mutect2 joint B", sample]))), unique(conf$sample))
SAMPLES
hw  <- dcast(conf, site ~ sample, value.var = "call")
vw  <- dcast(conf, site ~ sample, value.var = "VAF")
dw  <- dcast(conf, site ~ sample, value.var = "DPA")
gqw <- dcast(conf[!is.na(GQ)], site ~ sample, value.var = "GQ")
mqw <- dcast(conf[!is.na(MQM)], site ~ sample, value.var = "MQM")
hw[, c("CHROM", "POS") := .(sub(":.*", "", site), as.numeric(tstrsplit(site, ":", fixed = TRUE)[[2]]))]

# ---- 5. Pairwise statistics ----------------------------------------------------
# n_shared    sites confidently genotyped in both samples (all genotype combinations)
# IBS0        one ref, the other alt;   ibs0_per_1000 = 1000 * IBS0 / n_shared
# concordance same genotype, among shared sites where at least one sample is not ref
# thinned     the same after keeping one shared site per THIN_BP window (reduces dependence between
#             neighbouring SNPs in an inherited block - a simple stand-in for LD pruning)
# blocks      IBS0 sites grouped per chromosome when < BLOCK_GAP apart; n_chrom = chromosomes with IBS0
pair_stats <- function(x, y) {
  a <- hw[[x]]; b <- hw[[y]]
  ok <- !is.na(a) & !is.na(b)
  ibs0 <- ok & ((a == "ref" & b == "alt") | (a == "alt" & b == "ref"))
  nonref <- ok & !(a == "ref" & b == "ref")
  thin <- ok & !duplicated(paste(hw$CHROM, floor(hw$POS / THIN_BP), ok))
  d <- data.table(CHROM = hw$CHROM[ibs0], POS = hw$POS[ibs0])[order(CHROM, POS)]
  d[, block := cumsum(c(TRUE, diff(POS) > BLOCK_GAP)), by = CHROM]
  bs <- d[, .N, by = .(CHROM, block)]
  list(n_shared = sum(ok), IBS0 = sum(ibs0), ibs0_per_1000 = round(1000 * sum(ibs0) / sum(ok), 2),
       concordance = round(sum(nonref & a == b) / sum(nonref), 3),
       n_thinned = sum(thin), IBS0_thinned = sum(thin & ibs0),
       ibs0_per_1000_thinned = round(1000 * sum(thin & ibs0) / sum(thin), 2),
       n_blocks = nrow(bs), n_chrom = uniqueN(d$CHROM),
       largest_block_share = if (nrow(bs)) round(max(bs$N) / sum(bs$N), 2) else NA_real_)
}
pairs <- as.data.table(t(combn(SAMPLES, 2)))
setnames(pairs, c("s1", "s2"))
pairwise <- pairs[, pair_stats(s1, s2), by = .(s1, s2)]
pairwise[, `:=`(class1 = sample_class(s1), class2 = sample_class(s2))]
pairwise[, labelled_pair := mapply(function(a, b) isTRUE(KI_NORMALS[a] == b) || isTRUE(KI_NORMALS[b] == a), s1, s2)]
pairwise[order(ibs0_per_1000)]

# IBS0 sites of every pair (for the genome plot and block 8)
ibs0_sites <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  x <- pairs$s1[i]; y <- pairs$s2[i]; a <- hw[[x]]; b <- hw[[y]]
  k <- which(!is.na(a) & !is.na(b) & ((a == "ref" & b == "alt") | (a == "alt" & b == "ref")))
  if (!length(k)) return(NULL)
  data.table(pair = paste(x, "-", y), s1 = x, s2 = y, site = hw$site[k], CHROM = hw$CHROM[k], POS = hw$POS[k])
}))

# ---- 6. Verdict per tumour -----------------------------------------------------
tn <- rbind(pairwise[class1 == "tumour" & class2 != "tumour", .(tumour = s1, normal = s2, ibs0_per_1000, IBS0, n_shared, ibs0_per_1000_thinned, n_blocks, n_chrom)],
            pairwise[class2 == "tumour" & class1 != "tumour", .(tumour = s2, normal = s1, ibs0_per_1000, IBS0, n_shared, ibs0_per_1000_thinned, n_blocks, n_chrom)])
tn[, labelled := KI_NORMALS[tumour] == normal]
verdict <- tn[normal %chin% KI_NORMALS][order(tumour, ibs0_per_1000)][
  , .(labelled_normal = KI_NORMALS[tumour[1]], best_normal = normal[1], best_rate = ibs0_per_1000[1],
      labelled_rate = ibs0_per_1000[normal == KI_NORMALS[tumour[1]]][1], second_rate = ibs0_per_1000[2],
      second_to_best = round(ibs0_per_1000[2] / ibs0_per_1000[1], 1)), by = tumour]
verdict[, consistent := best_normal == labelled_normal]
verdict      # rates are IBS0 per 1000 shared confident sites; second_to_best = separation of the best match
fwrite(verdict, file.path(OUT, "verdict.tsv"), sep = "\t")

# ---- 7. Plots ------------------------------------------------------------------
# 7a. IBS0 rate heatmap, all samples (lower = more alike); labelled tumour/KI pairs outlined
hm <- rbind(pairwise[, .(a = s1, b = s2, ibs0_per_1000, n_shared, labelled_pair)],
            pairwise[, .(a = s2, b = s1, ibs0_per_1000, n_shared, labelled_pair)])
hm[, `:=`(a = factor(a, levels = SAMPLES), b = factor(b, levels = rev(SAMPLES)))]
p_hm <- ggplot(hm, aes(a, b, fill = ibs0_per_1000)) +
  geom_tile(colour = "white") +
  geom_tile(data = hm[labelled_pair == TRUE], colour = "black", linewidth = 1, fill = NA) +
  geom_text(aes(label = sprintf("%.1f\nn=%s", ibs0_per_1000, format(n_shared, big.mark = ","))), size = 2.6) +
  scale_fill_gradient(low = "#2E7D32", high = "white", trans = "sqrt") +
  labs(x = NULL, y = NULL, fill = "IBS0 per 1000\nshared sites",
       title = "Opposite homozygous genotypes (IBS0) per 1000 confidently genotyped shared sites",
       subtitle = sprintf("Lower = more alike. DP >= %d; ref <= %.2f, het %.2f-%.2f, alt >= %.2f; black box = labelled pair. WT: Mutect2 joint sites only",
                          MIN_DP, HOM_REF_MAX, HET_MIN, HET_MAX, HOM_ALT_MIN)) +
  coord_equal() + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_hm, "ibs0_rate_heatmap", 10, 9)

# 7b. Per tumour: IBS0 rate against each normal (the slide figure)
p_bar <- ggplot(tn, aes(normal, ibs0_per_1000, fill = labelled)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f\n(%s / %s)", ibs0_per_1000, IBS0, format(n_shared, big.mark = ","))), vjust = -0.15, size = 2.8) +
  facet_wrap(~ tumour, labeller = labeller(tumour = function(x) paste("Tumour", x)), scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = "#C62828", `FALSE` = "grey60"), labels = c(`TRUE` = "labelled as its normal", `FALSE` = "other normal")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "normal liver", y = "opposite homozygous genotypes\nper 1000 shared confident sites", fill = NULL,
       title = "Which normal liver has the same inherited genotypes as each tumour?",
       subtitle = "Lowest bar = most similar germline fingerprint. Label: rate (IBS0 sites / shared sites)") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_bar, "ibs0_rate_per_tumour", 12, 5.5)

# 7c. VAF vs VAF, tumour x KI normal, confident sites; IBS0 sites in red (max 20,000 points per panel)
set.seed(1)
sc <- rbindlist(lapply(names(TUMOURS), function(t) rbindlist(lapply(unname(KI_NORMALS), function(n) {
  if (!all(c(t, n) %chin% names(hw))) return(NULL)
  k <- which(!is.na(hw[[t]]) & !is.na(hw[[n]]))
  d <- data.table(tumour = t, normal = n, vaf_t = vw[[t]][k], vaf_n = vw[[n]][k],
                  ibs0 = (hw[[t]][k] == "ref" & hw[[n]][k] == "alt") | (hw[[t]][k] == "alt" & hw[[n]][k] == "ref"))
  rbind(d[ibs0 == TRUE], d[ibs0 == FALSE][sample(.N, min(.N, 20000))])
}))))
sc[, labelled := KI_NORMALS[tumour] == normal]
p_sc <- ggplot(sc, aes(vaf_n, vaf_t)) +
  geom_point(data = sc[ibs0 == FALSE], colour = "grey40", alpha = 0.15, size = 0.5) +
  geom_point(data = sc[ibs0 == TRUE], colour = "#C62828", size = 0.9) +
  geom_rect(data = unique(sc[labelled == TRUE, .(tumour, normal)]), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 1) +
  facet_grid(tumour ~ normal) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "KI normal VAF", y = "tumour VAF", title = "VAF at confidently genotyped shared sites",
       subtitle = "Red = opposite homozygous (corners). Normal ~0.5 with tumour ~0 or 1 can be LOH and is not by itself a mismatch. Box = labelled pair") +
  theme_bw()
save_plot(p_sc, "vaf_scatter_tumour_vs_normal", 10, 10)

# 7d. Where are the IBS0 sites? Tumour x KI pairs and KI x KI pairs
key_pairs <- pairwise[(class1 == "tumour" & class2 == "KI normal") | (class1 == "KI normal" & class2 == "KI normal"), paste(s1, "-", s2)]
gd <- ibs0_sites[pair %chin% key_pairs]
gd[, CHROM := factor(CHROM, levels = rev(AUTOSOMES))]
p_gen <- ggplot(gd, aes(POS / 1e6, CHROM)) +
  geom_point(shape = "|", size = 3, colour = "#C62828") +
  facet_wrap(~ pair, ncol = 3) +
  labs(x = "position (Mb)", y = "chromosome", title = "Genomic distribution of opposite homozygous sites",
       subtitle = "Spread over many chromosomes = many independent differences; few dense blocks = inherited segments / problem regions") +
  theme_bw() + theme(axis.text.y = element_text(size = 6))
save_plot(p_gen, "ibs0_genome_distribution", 14, 10)
blocks <- pairwise[, .(pair = paste(s1, "-", s2), IBS0, n_blocks, n_chrom, largest_block_share)][order(-IBS0)]
blocks

# 7e. Depth: tumour vs KI normal at shared sites (all sites as density, IBS0 sites as points)
dp <- rbindlist(lapply(names(TUMOURS), function(t) rbindlist(lapply(unname(KI_NORMALS), function(n) {
  if (!all(c(t, n) %chin% names(hw))) return(NULL)
  k <- which(!is.na(hw[[t]]) & !is.na(hw[[n]]))
  data.table(tumour = t, normal = n, dp_t = dw[[t]][k], dp_n = dw[[n]][k],
             ibs0 = (hw[[t]][k] == "ref" & hw[[n]][k] == "alt") | (hw[[t]][k] == "alt" & hw[[n]][k] == "ref"))
}))))
p_dp <- ggplot(dp, aes(dp_n, dp_t)) +
  geom_bin2d(bins = 60) + geom_point(data = dp[ibs0 == TRUE], colour = "#C62828", size = 0.6) +
  facet_grid(tumour ~ normal) + scale_fill_viridis_c(trans = "log10") +
  coord_cartesian(xlim = c(0, quantile(dp$dp_n, 0.99)), ylim = c(0, quantile(dp$dp_t, 0.99))) +
  labs(x = "KI normal depth (REF+ALT)", y = "tumour depth (REF+ALT)", fill = "sites",
       title = "Depth at shared confident sites (red = IBS0 sites)") + theme_bw()
save_plot(p_dp, "depth_tumour_vs_normal", 11, 10)
dp_summary <- dp[, .(median_dp_t = median(dp_t), median_dp_n = median(dp_n)), by = .(tumour, normal, ibs0)]
dp_summary

# 7f. Genotype quality (FreeBayes GQ, lower of the two samples) and mapping quality (MQM, lower of the two
#     files), IBS0 vs other shared sites
qq <- rbindlist(lapply(names(TUMOURS), function(t) rbindlist(lapply(unname(KI_NORMALS), function(n) {
  if (!all(c(t, n) %chin% names(hw)) || !all(c(t, n) %chin% names(gqw))) return(NULL)
  k <- which(!is.na(hw[[t]]) & !is.na(hw[[n]]) & !(hw[[t]] == "ref" & hw[[n]] == "ref"))
  g <- gqw[hw[k, .(site)], on = "site"]; m <- mqw[hw[k, .(site)], on = "site"]
  data.table(tumour = t, normal = n,
             group = fifelse((hw[[t]][k] == "ref" & hw[[n]][k] == "alt") | (hw[[t]][k] == "alt" & hw[[n]][k] == "ref"), "IBS0",
                             fifelse(hw[[t]][k] == hw[[n]][k], "same genotype", "het vs hom")),
             GQ_min = pmin(g[[t]], g[[n]]), MQM_min = pmin(m[[t]], m[[n]]))
}))))
qq_long <- melt(qq, id.vars = c("tumour", "normal", "group"), measure.vars = c("GQ_min", "MQM_min"), variable.name = "metric")
p_q <- ggplot(qq_long[!is.na(value)], aes(paste(tumour, "-", normal), value, fill = group)) +
  geom_boxplot(outlier.size = 0.2) + facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = NULL, fill = NULL, title = "Genotype quality (GQ) and mapping quality (MQM) at shared non-ref sites",
       subtitle = "IBS0 sites with lower quality than the rest point to technical artefacts") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
save_plot(p_q, "quality_ibs0_vs_other", 12, 8)
q_summary <- qq[, .(n = .N, median_GQ = median(GQ_min, na.rm = TRUE), median_MQM = median(MQM_min, na.rm = TRUE)), by = .(tumour, normal, group)]
q_summary

# 7g. Raw vs LD-thinned IBS0 rate, tumour x KI pairs
thin_dt <- melt(tn[normal %chin% KI_NORMALS], id.vars = c("tumour", "normal", "labelled"),
                measure.vars = c("ibs0_per_1000", "ibs0_per_1000_thinned"), variable.name = "set", value.name = "rate")
thin_dt[, set := fifelse(set == "ibs0_per_1000", "all shared sites", sprintf("1 site per %g Mb", THIN_BP / 1e6))]
p_thin <- ggplot(thin_dt, aes(normal, rate, fill = set)) +
  geom_col(position = "dodge", width = 0.7) + facet_wrap(~ tumour) +
  scale_fill_manual(values = c("grey40", "#90CAF9")) +
  labs(x = "KI normal", y = "IBS0 per 1000 shared sites", fill = NULL,
       title = "IBS0 rate with and without LD thinning", subtitle = "The order of normals should not change if the signal is genome-wide") +
  theme_bw() + theme(legend.position = "bottom")
save_plot(p_thin, "ibs0_rate_thinned", 10, 4.5)

# ---- 8. Optional: IBS0 sites vs Mouse Genomes Project strains ---------------------
strain_overlap <- data.table()
if (!is.na(MGP_VCF) && file.exists(MGP_VCF)) {
  reg <- unique(ibs0_sites[, .(CHROM = paste0(CHR_PREFIX, CHROM), POS)])
  rf <- tempfile(fileext = ".tsv"); fwrite(reg, rf, sep = "\t", col.names = FALSE)
  cmd <- paste(B, "view -R", shQuote(rf), "-s", paste(MGP_STRAINS, collapse = ","), shQuote(MGP_VCF), "-Ou |",
               B, "norm -m -any -Ou |", B, "query -f", shQuote("%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"))
  mgp <- fread(cmd = cmd, sep = "\t", header = FALSE, col.names = c("CHROM", "POS", "REF", "ALT", MGP_STRAINS), colClasses = "character")
  mgp[, site := paste(sub("^chr", "", CHROM), POS, REF, ALT, sep = ":")]
  for (s in MGP_STRAINS) mgp[, (s) := get(s) %chin% c("1/1", "1|1")]          # strain homozygous for this ALT
  strain_overlap <- ibs0_sites[pair %chin% key_pairs][mgp[, c("site", MGP_STRAINS), with = FALSE], on = "site", nomatch = NULL][
    , lapply(.SD, function(v) round(100 * mean(v), 1)), by = pair, .SDcols = MGP_STRAINS]
  strain_overlap <- ibs0_sites[pair %chin% key_pairs, .(IBS0 = .N), by = pair][strain_overlap, on = "pair"]
  print(strain_overlap)   # % of a pair's IBS0 sites where the strain is hom-alt for the same ALT
} else message("MGP_VCF not set: strain comparison (block 8) skipped")

# ---- 9. bcftools gtcheck (established tool) ---------------------------------------
# The FreeBayes files are filtered (SNVs, QUAL >= MIN_QUAL, autosomes) into OUT/tmp_gtcheck, indexed, merged
# and compared with bcftools gtcheck using the GT tag. Sites missing in a file are missing, not hom-ref.
gt_dir <- file.path(OUT, "tmp_gtcheck"); dir.create(gt_dir, showWarnings = FALSE)
autos <- paste(paste0(CHR_PREFIX, AUTOSOMES), collapse = ",")
tmp_bcf <- file.path(gt_dir, paste0(fb_files$tumour, ".bcf"))
for (i in seq_len(nrow(fb_files))) {
  system(paste(B, "view -v snps -t", autos, "-i", shQuote(paste0("QUAL>=", MIN_QUAL)), "-Ob -o", shQuote(tmp_bcf[i]), shQuote(fb_files$path[i])))
  system(paste(B, "index -f", shQuote(tmp_bcf[i])))
}
merged_bcf <- file.path(gt_dir, "freebayes_merged.bcf")
system(paste(B, "merge -Ob -o", shQuote(merged_bcf), paste(shQuote(tmp_bcf), collapse = " ")))
# bcftools >= 1.11: "-u GT", output DCv2 lines (query, genotyped, discordance, -logP(HWE), n_sites, n_matching)
# older bcftools:   "-G 1" (GTs only), output CN lines (discordance, n_sites, avg min depth, sample i, sample j)
bcf_ver <- as.numeric(strsplit(sub("^bcftools ([0-9]+\\.[0-9]+).*", "\\1", system2(BCFTOOLS, "--version", stdout = TRUE)[1]), ".", fixed = TRUE)[[1]])
new_gtcheck <- bcf_ver[1] > 1 || bcf_ver[2] >= 11
gtcheck_raw <- system(paste(B, "gtcheck", if (new_gtcheck) "-u GT" else "-G 1", shQuote(merged_bcf)), intern = TRUE)
writeLines(gtcheck_raw, file.path(OUT, "bcftools_gtcheck_raw.txt"))
gtcheck <- data.table()
dc <- grep("^DCv2\t", gtcheck_raw, value = TRUE)
cn <- grep("^CN\t", gtcheck_raw, value = TRUE)
if (length(dc)) {
  gtcheck <- fread(text = paste(dc, collapse = "\n"), header = FALSE, sep = "\t")
  setnames(gtcheck, 1:7, c("tag", "s1", "s2", "discordance", "avg_neg_log_p_hwe", "n_sites", "n_matching"))
} else if (length(cn)) {
  print(grep("^# \\[1\\]CN", gtcheck_raw, value = TRUE))          # column header of this bcftools version
  gtcheck <- fread(text = paste(cn, collapse = "\n"), header = FALSE, sep = "\t")
  setnames(gtcheck, 1:6, c("tag", "discordance", "n_sites", "avg_min_depth", "s1", "s2"))
}
if (nrow(gtcheck)) {
  gtcheck[, `:=`(s1 = short_name(s1), s2 = short_name(s2), discordance_per_1000 = round(1000 * as.numeric(discordance) / as.numeric(n_sites), 2))]
  print(gtcheck[order(discordance_per_1000), .(s1, s2, n_sites, discordance, discordance_per_1000)])
} else message("No DCv2 or CN lines in the gtcheck output - see bcftools_gtcheck_raw.txt")
unlink(gt_dir, recursive = TRUE)

# ---- 10. Outputs -------------------------------------------------------------
fwrite(pairwise, file.path(OUT, "pairwise_ibs0.tsv"), sep = "\t")
fwrite(tn, file.path(OUT, "tumour_vs_normal_ibs0.tsv"), sep = "\t")
fwrite(ibs0_sites, file.path(OUT, "ibs0_sites.tsv"), sep = "\t")
fwrite(blocks, file.path(OUT, "ibs0_blocks.tsv"), sep = "\t")
fwrite(dp_summary, file.path(OUT, "depth_summary.tsv"), sep = "\t")
fwrite(q_summary, file.path(OUT, "quality_summary.tsv"), sep = "\t")
if (nrow(strain_overlap)) fwrite(strain_overlap, file.path(OUT, "ibs0_strain_overlap.tsv"), sep = "\t")
if (nrow(gtcheck)) fwrite(gtcheck, file.path(OUT, "bcftools_gtcheck.tsv"), sep = "\t")
saveRDS(list(fb_files = fb_files, m2_files = m2_files, pairwise = pairwise, tn = tn, verdict = verdict, ibs0_sites = ibs0_sites,
             blocks = blocks, dp_summary = dp_summary, q_summary = q_summary, strain_overlap = strain_overlap, gtcheck = gtcheck,
             params = list(MIN_QUAL = MIN_QUAL, MIN_DP = MIN_DP, HOM_REF_MAX = HOM_REF_MAX, HET_MIN = HET_MIN, HET_MAX = HET_MAX,
                           HOM_ALT_MIN = HOM_ALT_MIN, THIN_BP = THIN_BP, BLOCK_GAP = BLOCK_GAP)),
        file.path(OUT, "sample_identity.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

md_table <- function(x) {
  x <- as.data.table(x)
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
writeLines(c(
  "# Sample identity (tumour / KI normal pairing)", "",
  "Genotype fingerprinting from VCF genotypes (no BAMs). Inferred genotypes under a diploid assumption.",
  "Rates are opposite homozygous genotypes (IBS0) per 1000 sites confidently genotyped in both samples.", "",
  "## Verdict per tumour", "", md_table(verdict), "",
  "## Tumour vs each normal", "", md_table(tn[order(tumour, ibs0_per_1000)]), "",
  "## IBS0 distribution (blocks, chromosomes)", "", md_table(blocks), "",
  if (nrow(gtcheck)) c("## bcftools gtcheck", "", md_table(gtcheck[order(discordance_per_1000), .(s1, s2, n_sites, discordance, discordance_per_1000)]), "") else "",
  "## Interpretation rules", "",
  "- A putative mismatch needs: a clearly higher IBS0 rate than the best match, IBS0 spread over many chromosomes and",
  "  independent blocks, no lower genotype/mapping quality at IBS0 sites, and the same order after LD thinning.",
  "- The baseline IBS0 in matched pairs is unexplained (genotype misclassification, mapping artefacts, somatic",
  "  mutation + LOH/amplification, allele-specific coverage), not necessarily sequencing error.",
  "- Many germline differences between mice of one line mean the colony is not isogenic (residual donor/129",
  "  background, B6J/B6N, incomplete backcrossing, segregating heterozygosity) or problem regions.",
  "- Until confirmed by sample records: a PUTATIVE sample mismatch.", "",
  "## Not done (data not available)", "",
  "- BAM-level pileup of the same sites (no BAMs).",
  "- Restriction to copy-neutral regions (no copy-number calls yet).",
  "- Strategies D (GM normal vs WT) and E (WT vs WT) as new caller runs (need BAMs).",
  if (!nrow(strain_overlap)) "- Mouse Genomes Project strain comparison (MGP_VCF not set)." else ""),
  file.path(OUT, "summary.md"))

git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_sample_identity"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  sprintf("- Parameters: MIN_QUAL=%d MIN_DP=%d ref<=%.2f het %.2f-%.2f alt>=%.2f THIN_BP=%g BLOCK_GAP=%g MGP_VCF=%s",
          MIN_QUAL, MIN_DP, HOM_REF_MAX, HET_MIN, HET_MAX, HOM_ALT_MIN, THIN_BP, BLOCK_GAP, MGP_VCF),
  paste0("- Inputs:     ", paste(basename(c(fb_files$path, m2_files$path)), collapse = ", ")),
  paste0("- Verdict:    ", paste(verdict[, sprintf("%s best=%s (%.2f/1000; labelled %s %.2f)", tumour, best_normal, best_rate,
                                                   labelled_normal, labelled_rate)], collapse = "; ")),
  paste0("- Outputs:    ", OUT),
  "", old_log), log_file)
message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
