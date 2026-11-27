#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-27
#
# for_patrick.R - three tables for Patrik:
#   1_identity_sites_all.csv         all germline identity sites (VAF / genome slides), wide: GT, REF and ALT reads per sample
#   2_igv_sites.csv                  20 discordant + 10 concordant sites for IGV snapshots
#   3_driver_genes_mouse_orthologs.csv  IntOGen driver genes (tumourevo) with their mouse orthologs
# Env on ssi: r-mutation. Output: results/mouse_wgs/for_patrick/

library(data.table)
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))

REPO      <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
MASKED    <- file.path(REPO, "results/mouse_wgs/mouse_wgs_somalier/masked/masked_6samples.bcf")
INTOGEN   <- "/common/WORK/pstancl/references/intogen/Compendium_Cancer_Genes.tsv"
ORTHOLOGS <- "/common/WORK/pstancl/references/cross_species/ensembl110/orthologs_hsapiens_mmusculus.tsv"
BCFTOOLS  <- "/common/WORK/pstancl/envs/r-mutation/bin/bcftools"
OUT       <- file.path(REPO, "results/mouse_wgs/for_patrick")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- c("HCC_1", "HCC_2", "HCC_3", "KI1", "KI2", "KI3")

# ---- 1. All identity sites, wide ----------------------------------------------
smp <- system(paste(BCFTOOLS, "query -l", MASKED), intern = TRUE)
smp
x <- fread(cmd = paste(BCFTOOLS, "query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%AD]\\n'", MASKED),
           header = FALSE, colClasses = "character",
           col.names = c("CHROM", "POS", "REF", "ALT", paste0(c("GT_", "AD_"), rep(smp, each = 2))))
x[, POS := as.integer(POS)]
for (s in SAMPLES) {
  ad <- tstrsplit(x[[paste0("AD_", s)]], ",", fixed = TRUE)
  x[, paste0("REF_reads_", s) := suppressWarnings(as.integer(ad[[1]]))]
  x[, paste0("ALT_reads_", s) := suppressWarnings(as.integer(ad[[2]]))]
}
# genotype "./." = masked (tumour VAF ambiguous); discordant = tumour genotype differs from its labelled normal
x[, disc_HCC_1_vs_KI1 := fifelse(GT_HCC_1 == "./.", NA, GT_HCC_1 != GT_KI1)]
x[, disc_HCC_2_vs_KI2 := fifelse(GT_HCC_2 == "./.", NA, GT_HCC_2 != GT_KI2)]
x[, disc_HCC_3_vs_KI3 := fifelse(GT_HCC_3 == "./.", NA, GT_HCC_3 != GT_KI3)]
sites <- x[order(suppressWarnings(as.integer(CHROM)), CHROM, POS),
           c("CHROM", "POS", "REF", "ALT", paste0(c("GT_", "REF_reads_", "ALT_reads_"), rep(SAMPLES, each = 3)),
             "disc_HCC_1_vs_KI1", "disc_HCC_2_vs_KI2", "disc_HCC_3_vs_KI3"), with = FALSE]
sites
nrow(sites)
# should be 522 / 2969, 519 / 2976, 1 / 2952 (genome slide)
sites[, .(HCC_1_vs_KI1 = paste(sum(disc_HCC_1_vs_KI1, na.rm = TRUE), "/", sum(!is.na(disc_HCC_1_vs_KI1))),
          HCC_2_vs_KI2 = paste(sum(disc_HCC_2_vs_KI2, na.rm = TRUE), "/", sum(!is.na(disc_HCC_2_vs_KI2))),
          HCC_3_vs_KI3 = paste(sum(disc_HCC_3_vs_KI3, na.rm = TRUE), "/", sum(!is.na(disc_HCC_3_vs_KI3))))]
fwrite(sites, file.path(OUT, "1_identity_sites_all.csv"))

# ---- 2. Sites for IGV --------------------------------------------------------------
# Criteria: all six samples genotyped; read fractions close to the genotype (max deviation <= 0.1);
# deepest first (lowest depth over the six samples); spread over chromosomes (one per chromosome first).
# Discordant: opposite homozygotes in BOTH labelled pairs (HCC_1 vs KI1 and HCC_2 vs KI2), so one snapshot
# shows both - the clearest evidence. Concordant: all six samples 0/1 (both alleles visible), no pair differs.
sel <- copy(sites)
dp  <- sapply(SAMPLES, function(s) sel[[paste0("REF_reads_", s)]] + sel[[paste0("ALT_reads_", s)]])
vaf <- sapply(SAMPLES, function(s) sel[[paste0("ALT_reads_", s)]]) / dp
gt  <- sapply(SAMPLES, function(s) sel[[paste0("GT_", s)]])
expected <- matrix(c(`0/0` = 0, `0/1` = 0.5, `1/1` = 1)[gt], nrow = nrow(gt))
sel[, min_depth := apply(dp, 1, min)]
sel[, max_vaf_dev := round(apply(abs(vaf - expected), 1, max), 3)]
TUM <- c("HCC_1", "HCC_2", "HCC_3"); NOR <- c("KI1", "KI2", "KI3")
sel[, min_depth_tumours := apply(dp[, TUM], 1, min)]
sel[, min_depth_normals := apply(dp[, NOR], 1, min)]
sel[, max_vaf_dev_tumours := round(apply(abs(vaf - expected)[, TUM], 1, max), 3)]
sel[, max_vaf_dev_normals := round(apply(abs(vaf - expected)[, NOR], 1, max), 3)]
sel[, opp_hom_both := GT_HCC_1 %chin% c("0/0", "1/1") & GT_KI1 %chin% c("0/0", "1/1") & disc_HCC_1_vs_KI1 %in% TRUE &
                      GT_HCC_2 %chin% c("0/0", "1/1") & GT_KI2 %chin% c("0/0", "1/1") & disc_HCC_2_vs_KI2 %in% TRUE]
sel[, all_het_concordant := apply(gt == "0/1", 1, all) %in% TRUE]
sel[, .(all_genotyped = sum(!is.na(min_depth)), opp_hom_both = sum(opp_hom_both), all_het_concordant = sum(all_het_concordant))]

pick <- function(d, n) {
  d <- d[!is.na(min_depth) & max_vaf_dev <= 0.1][order(-min_depth)]
  d[, rank_in_chr := seq_len(.N), by = CHROM]
  d[order(rank_in_chr, -min_depth)][seq_len(min(.N, n))]
}
disc <- pick(sel[opp_hom_both == TRUE], 20)
conc <- pick(sel[all_het_concordant == TRUE], 10)
igv <- rbind(disc[, group := "discordant (opposite homozygotes)"], conc[, group := "concordant (all 0/1)"])
igv[, igv_locus := paste0(CHROM, ":", POS - 50, "-", POS + 50)]
igv <- igv[, c("group", "igv_locus", names(sites), "min_depth_tumours", "min_depth_normals",
                "max_vaf_dev_tumours", "max_vaf_dev_normals"), with = FALSE]
igv[, comment_Patrik := ""]
igv[, .(group, igv_locus, GT_HCC_1, GT_KI1, GT_HCC_2, GT_KI2, GT_HCC_3, GT_KI3,
        min_depth_tumours, min_depth_normals, max_vaf_dev_tumours, max_vaf_dev_normals)]
fwrite(igv, file.path(OUT, "2_igv_sites.csv"))

# ---- 3. Driver genes with mouse orthologs ----------------------------------------
ig <- fread(INTOGEN)
genes <- sort(unique(ig$SYMBOL))
length(genes)
orth <- fread(ORTHOLOGS, colClasses = "character", na.strings = c("", "NA"))
drv <- orth[symbol_h %chin% genes & !is.na(gene_id_m),
            .(human_gene = symbol_h, human_ensembl = gene_id_h, mouse_gene = symbol_m, mouse_ensembl = gene_id_m, homology_type)]
drv <- merge(data.table(human_gene = genes), drv, by = "human_gene", all.x = TRUE)
drv
drv[, .N, by = homology_type]
drv[is.na(mouse_ensembl), human_gene]
fwrite(drv, file.path(OUT, "3_driver_genes_mouse_orthologs.csv"))

writeLines(c(paste("##", format(Sys.time(), "%Y-%m-%d %H:%M"), "- for_patrick.R, git",
                   system(paste("git -C", REPO, "rev-parse --short HEAD"), intern = TRUE)),
             paste("- Inputs:", MASKED, ";", INTOGEN, ";", ORTHOLOGS),
             paste("- Rows: sites", nrow(sites), "; IGV", nrow(igv), "; driver genes", length(genes), "(rows", nrow(drv), ")")),
           file.path(OUT, "RUN_LOG.md"))
