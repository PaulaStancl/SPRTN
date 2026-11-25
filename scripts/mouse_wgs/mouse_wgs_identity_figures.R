#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-25
#
# mouse_wgs_identity_figures.R - slide figures for the tumour / KI-normal identity result, from the outputs
# of mouse_wgs_somalier.R (no re-analysis). Primary combination only (panel A, 20 kb thinning).
#
#   fig0_site_selection   how the germline SNV panel was built (sites left after each filter step)
#   fig1_matrix           3 x 3 tumour x normal discordance at normal-homozygous sites (%, k / n, 95% CI);
#                         dashed = labelled pairing, solid green = best genetic match
#   fig2_vaf_scatter      VAF in the normal vs VAF in the tumour at the informative sites of every
#                         tumour x normal pair; same animal = points on the diagonal (0/0, 1/1)
#   fig3_genome           where the discordant sites lie along chromosomes 1-19, per tumour x normal pair
#
# Inputs (results/mouse_wgs/mouse_wgs_somalier/): pairwise_metrics.tsv, verdict.tsv, filter_counts.tsv,
# sites_used_primary.tsv, masked/masked_6samples.bcf (or .vcf; REF/ALT read counts).
# Env on ssi: r-mutation (data.table, ggplot2, bcftools). Run:
#   Rscript scripts/mouse_wgs/mouse_wgs_identity_figures.R [input folder]
# Output: <input folder>/slides/ (PNG 300 dpi + PDF, 16:9 friendly sizes)

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
IN       <- file.path(REPO, "results/mouse_wgs/mouse_wgs_somalier")
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
KI_NORMALS <- c(HCC_1 = "KI1", HCC_2 = "KI2", HCC_3 = "KI3")   # labelled pairing

# Colours: one sequential blue for the matrix; discordant = red AND a filled triangle (never colour alone)
BLUE_LO <- "#eef5fe"; BLUE_HI <- "#1757ac"
RED     <- "#dc3c3d"; GREY <- "#bdbdbd"; INK <- "#191919"; INK2 <- "#666666"
MATCH   <- "#008a17"                                             # outline of the genetic match

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) IN <- args[1]
OUT <- file.path(IN, "slides"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
}
theme_slide <- function(base = 13) theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e7e7e7", linewidth = 0.3),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = INK2),
        plot.caption = element_text(colour = INK2, hjust = 0), axis.text = element_text(colour = INK2),
        strip.text = element_text(colour = INK, face = "bold"), legend.position = "bottom",
        plot.title.position = "plot", plot.caption.position = "plot",
        plot.background = element_rect(fill = "white", colour = NA))
TUMOURS <- names(KI_NORMALS); NORMALS <- unname(KI_NORMALS)

pw <- fread(file.path(IN, "pairwise_metrics.tsv"))
vd <- fread(file.path(IN, "verdict.tsv"))[analysis == "primary"]
fc <- fread(file.path(IN, "filter_counts.tsv"))
su <- fread(file.path(IN, "sites_used_primary.tsv"), colClasses = list(character = "CHROM"))
PRIMARY <- pw[primary == TRUE, combo[1]]
GENETIC <- setNames(vd$best_normal, vd$tumour)[TUMOURS]
role_of <- function(t, n) fcase(n == KI_NORMALS[t] & n == GENETIC[t], "labelled = genetic match",
                                n == GENETIC[t], "genetic match", n == KI_NORMALS[t], "labelled", default = "other")

# ---- 2. Fig 0: site selection ------------------------------------------------
thin_lab <- sub("^A_", "", PRIMARY)
steps <- fc[scope == "all files" & !grepl("\\(info\\)", step) & !grepl("^panel [BC]", step) &
              (!grepl("^panel A, thinned", step) | step == paste0("panel A, thinned ", thin_lab))]
steps[, label := fcase(
  grepl("^union", step), "FreeBayes SNVs (3 pair files, biallelic, QUAL >= 30)",
  grepl("^present in", step), "site in all 3 pair files",
  grepl("^REF\\+ALT", step), sub("^REF\\+ALT >= ([0-9]+) in all ([0-9]+) samples.*", ">= \\1 reads in all \\2 samples", step),
  grepl("deeper than", step), "no sample > 2.5x its median depth (duplications)",
  grepl("MASK_BED", step), "outside masked regions",
  grepl("PASS record", step), "provider PASS filter",
  grepl("no KI normal ambiguous", step), "all 3 KI normals confidently genotyped",
  grepl(">= 1 KI normal", step), ">= 1 KI normal carries the variant (germline)",
  grepl("^panel A, thinned", step), paste0("thinned to >= ", thin_lab, " apart (final panel)"),
  default = step)]
steps <- steps[c(TRUE, diff(n_sites) != 0)]                     # drop steps that removed nothing (e.g. not applied)
steps[, label := factor(label, rev(unique(label)))]
p0 <- ggplot(steps, aes(n_sites, label)) +
  geom_col(fill = "#257de0", width = 0.7) +
  geom_text(aes(label = format(n_sites, big.mark = ",")), hjust = -0.1, size = 4, colour = INK) +
  scale_x_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.12))) +
  labs(x = "SNV sites", y = NULL, title = "Germline SNV panel used for the identity test",
       subtitle = "Genotype rules use the normal livers only; tumour genotypes are never used to keep or drop a site",
       caption = paste0("Genotypes from allele fractions: <= 0.05 hom-ref, 0.35-0.65 het, >= 0.95 hom-alt; ambiguous tumour genotypes masked.\n",
                        "Steps that removed no sites are not shown.")) +
  theme_slide() + theme(panel.grid.major.y = element_blank())
save_fig(p0, "fig0_site_selection", 11, 5.5)

# ---- 3. Fig 1: 3 x 3 discordance matrix ---------------------------------------
m <- pw[primary == TRUE & method == "custom" & type == "tumour-normal"]
m[, `:=`(tumour = fifelse(s1 %chin% TUMOURS, s1, s2), normal = fifelse(s1 %chin% TUMOURS, s2, s1))]
m[, role := role_of(tumour, normal)]
m[, lab := sprintf("%.1f%%\n%s / %s\n[%.1f-%.1f]", 100 * disc_rate, format(disc, big.mark = ","), format(n_inf, big.mark = ","), 100 * disc_lo, 100 * disc_hi)]
m[, `:=`(tumour = factor(tumour, rev(TUMOURS)), normal = factor(normal, NORMALS))]
ol_lab <- m[role %chin% c("labelled", "labelled = genetic match"), .(tumour, normal, what = "labelled pairing")]
ol_gen <- m[role %chin% c("genetic match", "labelled = genetic match"), .(tumour, normal, what = "best genetic match")]
p1 <- ggplot(m, aes(normal, tumour)) +
  geom_tile(aes(fill = 100 * disc_rate), colour = "white", linewidth = 1.5) +
  geom_text(aes(label = lab, colour = disc_rate > 0.12), size = 4.2, lineheight = 0.9, show.legend = FALSE) +
  geom_tile(data = ol_gen, aes(linetype = what), colour = MATCH, fill = NA, linewidth = 1.6, width = 0.94, height = 0.94) +
  geom_tile(data = ol_lab, aes(linetype = what), colour = INK, fill = NA, linewidth = 1.1, width = 0.86, height = 0.86) +
  scale_fill_gradient(low = BLUE_LO, high = BLUE_HI, name = "discordant sites (%)") +
  scale_colour_manual(values = c(`FALSE` = INK, `TRUE` = "white")) +
  scale_linetype_manual(values = c(`labelled pairing` = "22", `best genetic match` = "solid"), name = NULL) +
  guides(linetype = guide_legend(override.aes = list(colour = c(MATCH, INK)))) +
  labs(x = "KI normal liver", y = "Tumour",
       title = "Germline discordance: tumour vs each normal liver",
       subtitle = "SNVs where the normal is homozygous; % with a different tumour genotype\n(discordant / informative sites, [95% CI])",
       caption = sprintf("0%% = same animal. Dashed: sample labels. Green: best genetic match\n(paired exact sign test vs next best normal: %s).",
                         paste(sprintf("%s p = %.0e", vd$tumour, vd$paired_p_best_vs_second), collapse = ", "))) +
  coord_equal() + theme_slide() + theme(panel.grid = element_blank())
save_fig(p1, "fig1_matrix", 9, 7.5)

# ---- 4. Fig 2: VAF scatter per tumour x normal ---------------------------------
mf <- if (file.exists(file.path(IN, "masked", "masked_6samples.bcf"))) file.path(IN, "masked", "masked_6samples.bcf") else file.path(IN, "masked", "masked_6samples.vcf")
if (grepl("\\.bcf$", mf)) {
  smp <- system2(BCFTOOLS, c("query", "-l", shQuote(mf)), stdout = TRUE)
  ad <- fread(cmd = paste(shQuote(BCFTOOLS), "query -f", shQuote("%CHROM\\t%POS[\\t%AD]\\n"), shQuote(mf)),
              header = FALSE, col.names = c("CHROM", "POS", smp), colClasses = "character")
} else {
  v <- fread(mf, skip = "#CHROM", colClasses = "character"); setnames(v, 1, "CHROM")
  smp <- names(v)[-(1:9)]
  ad <- v[, c("CHROM", "POS", smp), with = FALSE]
  for (s in smp) set(ad, j = s, value = sub("^[^:]*:", "", ad[[s]]))        # GT:AD -> AD
}
ad[, POS := as.integer(POS)]
vaf_of <- function(x) { p <- tstrsplit(x, ",", fixed = TRUE); r <- suppressWarnings(as.numeric(p[[1]])); a <- suppressWarnings(as.numeric(p[[2]])); a / (r + a) }
su[, POS := as.integer(POS)]
sv <- ad[su[, c("CHROM", "POS", TUMOURS, NORMALS), with = FALSE], on = .(CHROM, POS)]
pts <- rbindlist(lapply(TUMOURS, function(t) rbindlist(lapply(NORMALS, function(n) {
  g_t <- sv[[paste0("i.", t)]]; g_n <- sv[[paste0("i.", n)]]
  inf <- !is.na(g_t) & !is.na(g_n) & g_n != 1L
  data.table(tumour = t, normal = n, chr = sv$CHROM[inf], pos = sv$POS[inf],
             vaf_n = vaf_of(sv[[n]][inf]), vaf_t = vaf_of(sv[[t]][inf]), discordant = (g_t != g_n)[inf])
}))))
pts[, role := role_of(tumour, normal)]
pts[, panel := paste0(tumour, " vs ", normal, "\n(", role, ")")]
ann <- pts[, .(n = .N, k = sum(discordant)), by = .(tumour, normal, panel)][, txt := sprintf("%s / %s discordant", format(k, big.mark = ","), format(n, big.mark = ","))]
pts[, `:=`(tumour = factor(tumour, TUMOURS), normal = factor(normal, NORMALS))]
ann[, `:=`(tumour = factor(tumour, TUMOURS), normal = factor(normal, NORMALS))]
set.seed(1)
p2 <- ggplot(pts, aes(vaf_n, vaf_t)) +
  geom_abline(slope = 1, intercept = 0, colour = "#c2c2c2", linewidth = 0.4) +
  geom_jitter(data = pts[discordant == FALSE], colour = GREY, size = 0.9, alpha = 0.6, width = 0.01, height = 0, shape = 16) +
  geom_jitter(data = pts[discordant == TRUE], colour = RED, size = 1.6, alpha = 0.8, width = 0.01, height = 0, shape = 17) +
  geom_text(data = ann, aes(x = 0.5, y = 1.12, label = txt), size = 3.4, colour = INK) +
  facet_grid(tumour ~ normal, labeller = labeller(tumour = function(x) paste("Tumour", x), normal = function(x) paste("Normal", x))) +
  geom_rect(data = unique(pts[role != "other", .(tumour, normal, role)]),
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, linetype = fifelse(grepl("genetic", role), "best genetic match", "labelled pairing")),
            inherit.aes = FALSE, fill = NA, colour = INK, linewidth = 0.8) +
  scale_linetype_manual(values = c(`labelled pairing` = "22", `best genetic match` = "solid"), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(-0.05, 1.05)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.05, 1.18)) +
  labs(x = "allele fraction in the normal liver (sites where the normal is homozygous)", y = "allele fraction in the tumour",
       title = "Allele fractions at informative germline SNVs",
       subtitle = "Same animal: tumour stays at 0 or 1 like the normal. Different animal: red triangles = tumour het or opposite homozygote",
       caption = "Grey circles: same genotype; red triangles: discordant. Solid frame: best genetic match; dashed frame: labelled pairing.") +
  theme_slide(12) + theme(panel.border = element_rect(colour = "#e7e7e7", fill = NA))
save_fig(p2, "fig2_vaf_scatter", 11, 9)

# ---- 5. Fig 3: discordant sites along the genome -------------------------------
chrs <- su[, .(len = max(POS)), by = .(chr_n = as.integer(sub("^chr", "", CHROM)))][order(chr_n)]
chrs[, offset := c(0, cumsum(as.numeric(len))[-.N])]
chrs[, mid := offset + len / 2]
pts[, chr_n := as.integer(sub("^chr", "", chr))]
pts[chrs, on = "chr_n", gpos := pos + i.offset]
rows <- CJ(tumour = TUMOURS, normal = NORMALS, sorted = FALSE)
rows[, row := paste0(tumour, " vs ", normal, "  ", fcase(role_of(tumour, normal) == "genetic match", "(genetic match)",
                                                         role_of(tumour, normal) == "labelled", "(labelled)",
                                                         role_of(tumour, normal) == "labelled = genetic match", "(labelled = match)", default = ""))]
pts[rows, on = .(tumour, normal), row := i.row]
pts[, row := factor(row, rev(rows$row))]
bands <- chrs[seq(1, .N, by = 2)]
p3 <- ggplot() +
  geom_rect(data = bands, aes(xmin = offset, xmax = offset + len, ymin = -Inf, ymax = Inf), fill = "#f3f3f3") +
  geom_point(data = pts[discordant == TRUE], aes(gpos, row), shape = "|", size = 4, colour = RED) +
  scale_x_continuous(breaks = chrs$mid, labels = chrs$chr_n, expand = c(0, 0)) +
  labs(x = "chromosome", y = NULL, title = "Where the discordant germline SNVs lie",
       subtitle = "Each tick = one SNV where the tumour genotype differs from the normal (normal homozygous). Same animal: an empty row",
       caption = sprintf("Panel A, %s thinning. Chromosomes with discordant SNVs: %s.", sub("^A_", "", PRIMARY),
                         paste(pts[discordant == TRUE, .(n = uniqueN(chr_n)), by = .(tumour, normal)][order(tumour, normal), sprintf("%s-%s %d", tumour, normal, n)], collapse = ", "))) +
  theme_slide(12) + theme(panel.grid = element_blank(), axis.text.y = element_text(colour = INK))
save_fig(p3, "fig3_genome", 12, 5)

message("Figures written to ", OUT)
