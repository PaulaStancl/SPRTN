# ---------------------------------------------------------------------------
# 07_roq_damage.R - are Mutect2's low-ROQ calls DNA damage? Deamination (C>T/G>A) vs oxidation
# (8-oxoG, G>T/C>A) vs a general library / orientation problem, for RJALS (one tumour, one normal).
# Run:  Rscript 07_roq_damage.R   (r-mutation env)
#   1st run: variant-level steps; writes read_level/sites.tsv for the read-level step
#   then:    cd ../wgs && qsub qsub_damage.sh   (Picard + GIV on tumour and normal, read positions)
#   2nd run: the same, now with the sample metrics and read-level evidence in the report
#
# Groups (tumour SNVs, Mutect2 FORMAT/F1R2 and F2R1 = reads per allele in each pair orientation):
#   PASS ROQ < cut / PASS ROQ >= cut   cut = SPRTN_ROQ_CUT (default 20, as in 01)
#   orientation-filtered               FILTER contains "orientation" (alone or with other filters)
# Steps:
#   0  header check: ROQ, F1R2/F2R1, AD, AF, the orientation filter, FilterMutectCalls run with
#      --ob-priors (else ROQ is not interpreted). SNVs only; multiallelics split (06's bcftools norm).
#   1  per-variant table: context, CpG, mutation12 (REF-based, C>T and G>A apart), mutation6
#   2  orientation: signed_orient = F1R2_alt / (F1R2_alt + F2R1_alt), 0.5 = balanced.
#      Main test: counts pooled per group and class (binomial vs 0.5; Fisher vs the REF reads'
#      orientation at the same sites, which removes any library-wide imbalance; sign test over
#      variants, which is robust to read overdispersion). Per-variant binomial only with >= 5
#      ALT orientation reads (BH, orient_biased at q < 0.05). Damage on one strand gives
#      OPPOSITE skew in complementary classes (C>T vs G>A, G>T vs C>A) - orientation_mirror.tsv.
#      C>T / G>A also split CpG vs non-CpG: a CpG excess alone is what real SBS1 looks like, so
#      it only counts together with an orientation skew.
#   3  plots (below)
#   4  read level (after qsub_damage.sh): ALT vs REF reads - distance to fragment end, cycle, BQ, MQ
#   6  tumour vs normal: Picard pre-adapter Q (Deamination, OxoG) and GIV (qsub_damage.sh)
#   7  96-context per group; cosine to SBS1, SBS5, SBS45 (COSMIC, from the sigprofiler env), the
#      Guo et al. 2022 FFPE signatures (ref/FFPE_signatures_Guo2022.tsv, repaired / unrepaired)
#      and to the other groups
# Output: <OUT>/qc_damage/
#   header_check.txt, variant_table.tsv, orientation_by_class.tsv, orientation_mirror.tsv,
#   orientation_cpg.tsv, class_excess.tsv, signature_cosine.tsv, sbs96_groups.txt,
#   sample_summary.tsv, read_level_summary.tsv, report.md,
#   chance_one_sided_by_class.tsv / _by_n.tsv / .pdf, tc_ag_by_context.tsv, alt_filter_vs_roq.tsv (section 2b)
#   orientation_pooled.pdf (main), orientation_hist.pdf, vaf_vs_roq.pdf, roq_by_class.pdf,
#   orient_biased_by_class.pdf, vaf_ct_by_group.pdf, signature_cosine.pdf, snv_96context_by_group.pdf,
#   read_level.pdf;  read_level/sites.tsv (input for 13_damage_metrics.sh)
# ---------------------------------------------------------------------------
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
od <- file.path(OUT, "qc_damage"); dir.create(file.path(od, "read_level"), recursive = TRUE, showWarnings = FALSE)
ROQ_CUT <- as.numeric(Sys.getenv("SPRTN_ROQ_CUT", "20"))
MIN_N   <- 5                                    # ALT orientation reads for the per-variant test
SITES_MAX <- 2000                               # per group, for the read-level step
MIRROR_MIN <- 0.2                               # F1R2 shares of the two classes at least this far apart
G_LO <- sprintf("PASS ROQ < %g", ROQ_CUT); G_HI <- sprintf("PASS ROQ >= %g", ROQ_CUT); G_OR <- "orientation-filtered"
GROUPS <- c(G_LO, G_HI, G_OR)
GCOL <- setNames(c("#d62728", "grey50", "#1f77b4"), GROUPS)
PAIRS <- list(c("C>T", "G>A"), c("G>T", "C>A"), c("T>C", "A>G"), c("T>A", "A>T"), c("C>G", "G>C"), c("T>G", "A>C"))
M12 <- c(vapply(PAIRS, `[`, "", 1), vapply(PAIRS, `[`, "", 2))   # facet order: row 1 one strand, row 2 its complement
comp <- c(A = "T", C = "G", G = "C", T = "A")
# two-sided exact binomial p vs 0.5 (symmetric, so 2 x the smaller tail), vectorised
binom_p <- function(k, n) fifelse(n > 0, pmin(1, 2 * pmin(pbinom(k, n, 0.5), pbinom(k - 1, n, 0.5, lower.tail = FALSE))), NA_real_)

# ---- 0. input + header check ---------------------------------------------------------------
nd  <- file.path(SAREK, "normalized_bcftools", PAIR)
vcf <- if (dir.exists(nd)) find_one(nd, "\\.mutect2\\.filtered\\.norm\\.vcf\\.gz$", required = FALSE) else NA
if (is.na(vcf)) vcf <- find_one(file.path(SAREK, "variant_calling", "mutect2", PAIR), "\\.mutect2\\.filtered\\.vcf\\.gz$")
message("Mutect2: ", vcf, if (!grepl("\\.norm\\.", vcf)) "  (not normalised: multiallelic records are dropped, not split)")
con <- gzfile(vcf, "r"); hdr <- character(0)
repeat { l <- readLines(con, 2000); hdr <- c(hdr, l); if (!length(l) || any(startsWith(l, "#CHROM"))) break }
close(con); hdr <- hdr[startsWith(hdr, "##")]
def <- function(kind, id) grep(sprintf("^##%s=<ID=%s,", kind, id), hdr, value = TRUE)
fmc <- grep("^##GATKCommandLine=<ID=FilterMutectCalls", hdr, value = TRUE)
chk <- data.table(field = c("INFO/ROQ", "FORMAT/F1R2", "FORMAT/F2R1", "FORMAT/AD", "FORMAT/AF", "FILTER/orientation",
                            "FilterMutectCalls with --ob-priors"),
                  found = c(length(def("INFO", "ROQ")) > 0, length(def("FORMAT", "F1R2")) > 0, length(def("FORMAT", "F2R1")) > 0,
                            length(def("FORMAT", "AD")) > 0, length(def("FORMAT", "AF")) > 0, length(def("FILTER", "orientation")) > 0,
                            any(grepl("--ob-priors|--orientation-bias-artifact-priors", fmc))))
OB_UNKNOWN <- !length(fmc)                     # no FilterMutectCalls line kept in the header
ROQ_OK <- chk[field == "INFO/ROQ", found] && (chk[field == "FilterMutectCalls with --ob-priors", found] ||
                                               (OB_UNKNOWN && chk[field == "FILTER/orientation", found]))
writeLines(c(paste("#", vcf), "# header definitions used (no records):",
             unlist(lapply(c("ROQ", "F1R2", "F2R1", "AD", "AF"), function(i) c(def("INFO", i), def("FORMAT", i)))),
             def("FILTER", "orientation"), "",
             sprintf("%-36s %s", chk$field, fifelse(chk$found, "yes", "NO")),
             if (OB_UNKNOWN) "FilterMutectCalls command line not in the header - --ob-priors assumed from the ROQ / orientation definitions",
             sprintf("ROQ interpretable: %s", ROQ_OK)), file.path(od, "header_check.txt"))
cat(readLines(file.path(od, "header_check.txt")), sep = "\n")
if (!all(chk[field %like% "F1R2|F2R1", found])) stop("no F1R2 / F2R1 in the Mutect2 VCF - orientation cannot be tested", call. = FALSE)
if (!ROQ_OK) warning("ROQ not interpretable (no --ob-priors): groups by ROQ are reported but not interpreted")

# PASS + orientation-filtered single-base records only (the full VCF has ~1.25 M records)
tmp <- tempfile(fileext = ".vcf")
system(sprintf("gzip -dc %s | awk -F'\\t' '/^#CHROM/ || (!/^##/ && ($7 == \"PASS\" || $7 ~ /orientation/) && length($4) == 1 && length($5) == 1)' > %s",
               shQuote(vcf), shQuote(tmp)))
d <- read_vcf_table(tmp)[TYPE == "SNV"]; unlink(tmp)

# ---- 1. per-variant table -----------------------------------------------------------------
two <- function(x, k) suppressWarnings(as.numeric(tstrsplit(as.character(x), ",", fixed = TRUE, keep = k)[[1]]))
d[, `:=`(ROQ = suppressWarnings(as.numeric(info_ROQ)),
         F1R2_ref = two(t_F1R2, 1), F1R2_alt = two(t_F1R2, 2), F2R1_ref = two(t_F2R1, 1), F2R1_alt = two(t_F2R1, 2),
         DP = suppressWarnings(as.numeric(t_DP)))]
d[, group := fcase(FILTER == "PASS" & !is.na(ROQ) & ROQ <  ROQ_CUT, G_LO,
                   FILTER == "PASS" & !is.na(ROQ) & ROQ >= ROQ_CUT, G_HI,
                   grepl("orientation", FILTER, fixed = TRUE),      G_OR)]
if (d[FILTER == "PASS" & is.na(ROQ), .N]) message(d[FILTER == "PASS" & is.na(ROQ), .N], " PASS SNVs without ROQ - left out")
d <- d[!is.na(group)][, group := factor(group, levels = GROUPS)]
d[, orientation_only := FILTER == "orientation"]
d[, NC_3 := trinuc_context(d)]
d[, `:=`(mutation12 = paste0(REF, ">", ALT1),
         mutation6  = fifelse(REF %chin% c("C", "T"), paste0(REF, ">", ALT1), paste0(comp[REF], ">", comp[ALT1])),
         CpG = (REF == "C" & ALT1 == "T" & substr(NC_3, 3, 3) == "G") | (REF == "G" & ALT1 == "A" & substr(NC_3, 1, 1) == "C"))]
d[, IN_DBS := in_doublet(.SD, "group")]
d[!IN_DBS & substr(NC_3, 2, 2) == REF, SBS96 := sbs96(REF, ALT1, NC_3)]

# ---- 2. orientation --------------------------------------------------------------------------
d[, `:=`(orient_n = F1R2_alt + F2R1_alt, signed_orient = F1R2_alt / (F1R2_alt + F2R1_alt),
         ref_orient = F1R2_ref / (F1R2_ref + F2R1_ref))]
d[, tested := !is.na(orient_n) & orient_n >= MIN_N]
d[tested == TRUE, binom_p := binom_p(F1R2_alt, orient_n)]
d[tested == TRUE, binom_q := p.adjust(binom_p, "BH")]
d[, orient_biased := tested & !is.na(binom_q) & binom_q < 0.05]
excl <- d[, .(n = .N, tested = sum(tested), excluded = sum(!tested)), by = group][order(group)]
message("per-variant test needs >= ", MIN_N, " ALT orientation reads:"); print(excl)
vt <- d[, .(sample = TUMOUR, chrom = CHROM, pos = POS, REF, ALT = ALT1, FILTER, group, orientation_only, AF = VAF, DP,
            alt_depth = T_ALT, ROQ, F1R2_ref, F1R2_alt, F2R1_ref, F2R1_alt, context = NC_3, CpG, mutation12, mutation6,
            SBS96, signed_orient, ref_orient, orient_n, binom_p = signif(binom_p, 3), binom_q = signif(binom_q, 3), orient_biased)]
fwrite(vt, file.path(od, "variant_table.tsv"), sep = "\t")

# pooled per group x class: is the ALT skewed, and more than the REF reads at the same sites?
pool <- function(x) x[, {
  a1 <- sum(F1R2_alt, na.rm = TRUE); a2 <- sum(F2R1_alt, na.rm = TRUE)
  r1 <- sum(F1R2_ref, na.rm = TRUE); r2 <- sum(F2R1_ref, na.rm = TRUE)
  maj1 <- sum(F1R2_alt > F2R1_alt, na.rm = TRUE); maj2 <- sum(F2R1_alt > F1R2_alt, na.rm = TRUE)
  .(n_var = .N, n_tested = sum(tested), F1R2_alt = a1, F2R1_alt = a2, frac_F1R2_alt = a1 / (a1 + a2),
    ci_lo = qbeta(0.025, a1, a2 + 1), ci_hi = qbeta(0.975, a1 + 1, a2),
    frac_F1R2_ref = r1 / (r1 + r2), p_binom = binom_p(a1, a1 + a2),
    p_vs_ref = if (a1 + a2 > 0 && r1 + r2 > 0) fisher.test(matrix(c(a1, a2, r1, r2), 2))$p.value else NA_real_,
    var_F1R2_major = maj1, var_F2R1_major = maj2, p_sign = binom_p(maj1, maj1 + maj2),
    pct_orient_biased = round(100 * mean(orient_biased[tested]), 1))
}, by = c("group", setdiff(names(x)[names(x) %chin% c("mutation12", "CpG_class")], "group"))]
adj <- function(t) { for (p in c("p_binom", "p_vs_ref", "p_sign")) set(t, j = sub("^p_", "q_", p), value = signif(p.adjust(t[[p]], "BH"), 3)); t }
oc <- adj(pool(d))[order(group, match(mutation12, M12))]
print(oc[, .(group, mutation12, n_var, n_tested, frac_F1R2_alt = round(frac_F1R2_alt, 3), frac_F1R2_ref = round(frac_F1R2_ref, 3), q_vs_ref, q_sign)])
fwrite(oc, file.path(od, "orientation_by_class.tsv"), sep = "\t")

# mirrored pairs: opposite skew in the two complementary classes
mir <- rbindlist(lapply(PAIRS, function(pr) oc[mutation12 %chin% pr, {
  a <- .SD[mutation12 == pr[1]]; b <- .SD[mutation12 == pr[2]]
  if (!nrow(a) || !nrow(b)) NULL else
  .(pair = paste(pr, collapse = " | "), frac_1 = a$frac_F1R2_alt, frac_2 = b$frac_F1R2_alt,
    n_1 = a$n_var, n_2 = b$n_var,
    p = fisher.test(matrix(c(a$F1R2_alt, a$F2R1_alt, b$F1R2_alt, b$F2R1_alt), 2))$p.value)
}, by = group]))
mir[, q := signif(p.adjust(p, "BH"), 3)][, opposite := (frac_1 - 0.5) * (frac_2 - 0.5) < 0]
mir[, mirrored := opposite & q < 0.05 & abs(frac_1 - frac_2) >= MIRROR_MIN]   # pooled reads: tiny skews are "significant"
print(mir); fwrite(mir, file.path(od, "orientation_mirror.tsv"), sep = "\t")

# C>T / G>A at CpG vs elsewhere
ct <- d[mutation12 %chin% c("C>T", "G>A")][, CpG_class := fifelse(CpG, "CpG", "non-CpG")]
ocg <- adj(pool(ct))[order(group, mutation12, CpG_class)]
fwrite(ocg, file.path(od, "orientation_cpg.tsv"), sep = "\t")

# class excess vs the clean group (pyrimidine classes, as in 01)
cl <- dcast(d[, .N, by = .(group, mutation6)], mutation6 ~ group, value.var = "N", fill = 0)
for (g in setdiff(GROUPS, names(cl))) set(cl, j = g, value = 0L)
ex <- rbindlist(lapply(c(G_LO, G_OR), function(g) {
  tg <- sum(cl[[g]]); th <- sum(cl[[G_HI]])
  cl[, .(group = g, mutation6, n = get(g), pct = round(100 * get(g) / tg, 1), pct_clean = round(100 * get(G_HI) / th, 1),
         p = mapply(function(a, b) if (tg && th) fisher.test(matrix(c(a, tg - a, b, th - b), 2))$p.value else NA_real_, get(g), get(G_HI)))]
}))
ex[, q := signif(p.adjust(p, "BH"), 3)]
print(ex); fwrite(ex, file.path(od, "class_excess.tsv"), sep = "\t")

# ---- 2b. follow-up checks ------------------------------------------------------------------------
# (a) Chance: with n ALT orientation reads split 50/50, P(all on one side) = 2 * 0.5^n. Observed vs
#     expected one-sided calls per class (n >= 2; n = 1 is always one-sided) and per n. Both groups
#     are SELECTED on one-sidedness (ROQ / the orientation filter), so some excess is built in; the
#     telling part is n >= 6, where chance one-sidedness is <= 3%.
ch <- d[group %chin% c(G_LO, G_OR) & !is.na(orient_n) & orient_n >= 1]
ch[, `:=`(one_sided = F1R2_alt == 0 | F2R1_alt == 0, p_exp = pmin(1, 2 * 0.5^orient_n),
          n_bin = factor(fcase(orient_n <= 5, as.character(orient_n), orient_n <= 9, "6-9", default = ">=10"),
                         levels = c(1:5, "6-9", ">=10")))]
chance <- function(x, by) x[, .(n_var = .N, observed = sum(one_sided), expected = sum(p_exp),
                                sd = sqrt(sum(p_exp * (1 - p_exp)))), by = by][
  , `:=`(obs_exp = round(observed / expected, 2), z = fifelse(sd > 0, round((observed - expected) / sd, 1), NA_real_),
         expected = round(expected, 1))][, sd := NULL][]
ch_cls <- rbind(chance(ch[orient_n >= 2], c("group", "mutation12")), chance(ch[orient_n >= 2], "group")[, mutation12 := "all"])
ch_cls <- ch_cls[order(group, match(mutation12, c(M12, "all")))]
ch_n <- chance(ch, c("group", "n_bin"))[order(group, n_bin)]
print(ch_n); fwrite(ch_cls, file.path(od, "chance_one_sided_by_class.tsv"), sep = "\t")
fwrite(ch_n, file.path(od, "chance_one_sided_by_n.tsv"), sep = "\t")
save_plot(ggplot(ch_n, aes(n_bin, observed / n_var, colour = group, group = group)) + geom_line() + geom_point(aes(size = n_var)) +
            geom_line(aes(y = expected / n_var), linetype = 2) + scale_colour_manual(values = GCOL) +
            scale_y_continuous(labels = function(x) paste0(100 * x, "%"), limits = c(0, 1)) +
            labs(x = "ALT orientation reads (F1R2_alt + F2R1_alt)", y = "calls with all ALT reads on one side", colour = NULL, size = "calls",
                 title = "One-sided ALT support: observed (solid) vs chance 2 x 0.5^n (dashed)",
                 subtitle = "both groups are selected on one-sidedness; chance is <= 3% from n = 6"),
          "chance_one_sided", od, w = 8, h = 4.5)

# (b) T>C / A>G calls by trinucleotide context (+ strand, REF-based), orientation-filtered vs clean
tc <- d[group %chin% c(G_OR, G_HI) & mutation12 %chin% c("T>C", "A>G") & !is.na(NC_3),
        .(n = .N, F1R2_alt = sum(F1R2_alt, na.rm = TRUE), F2R1_alt = sum(F2R1_alt, na.rm = TRUE)), by = .(group, mutation12, context = NC_3)]
tc[, `:=`(pct_of_class = round(100 * n / sum(n), 1), frac_F1R2_alt = round(F1R2_alt / (F1R2_alt + F2R1_alt), 3)), by = .(group, mutation12)]
tc <- dcast(tc, mutation12 + context ~ group, value.var = c("n", "pct_of_class", "frac_F1R2_alt"), fill = NA)
setnames(tc, gsub(paste0("_", G_OR, "$"), "_filtered", gsub(paste0("_", G_HI, "$"), "_clean", names(tc), fixed = FALSE), fixed = FALSE))
if ("n_filtered" %in% names(tc)) setorder(tc, mutation12, -n_filtered)
fwrite(tc, file.path(od, "tc_ag_by_context.tsv"), sep = "\t")

# (c) ROQ decision: a minimum-ALT-read filter vs ROQ < 30 on the PASS SNVs - how many each removes,
#     their overlap, the share C>T, and (with 03's mutation_callers.csv) how many removed calls
#     another caller confirms. Nothing is filtered here.
d[, KEY := paste(CHROM, POS, REF, ALT1, sep = ":")]
mcf <- file.path(OUT, "comparison", "overlap", "mutation_callers.csv")
if (file.exists(mcf)) d[fread(mcf, select = c("KEY", "N_CALLERS")), on = "KEY", confirmed := i.N_CALLERS >= 2]
if (!"confirmed" %in% names(d)) { d[, confirmed := NA]; message("no ", mcf, " - confirmation left out (run 03)") }
ps <- d[FILTER == "PASS"]
r30 <- ps$ROQ < 30
rules <- c(list("ROQ < 30" = r30, "ROQ < 30 unless confirmed" = r30 & !(ps$confirmed %in% TRUE)),
           setNames(lapply(3:6, function(k) ps$T_ALT < k), sprintf("ALT reads < %d", 3:6)))
af <- rbindlist(lapply(names(rules), function(r) { x <- rules[[r]] %in% TRUE
  data.table(rule = r, removed = sum(x), pct_of_pass = round(100 * mean(x), 1),
             removed_confirmed = sum(x & ps$confirmed %in% TRUE), pct_CT = round(100 * mean(ps$mutation6[x] == "C>T"), 1),
             also_ROQ30 = sum(x & r30), only_this = sum(x & !r30), only_ROQ30 = sum(!x & r30)) }))
print(af); fwrite(af, file.path(od, "alt_filter_vs_roq.tsv"), sep = "\t")

# ---- 3. plots ---------------------------------------------------------------------------------
d[, m12 := factor(mutation12, levels = M12)]; oc[, m12 := factor(mutation12, levels = M12)]
sub_excl <- paste(excl[, sprintf("%s: %d of %d", group, tested, n)], collapse = "; ")
save_plot(ggplot(oc, aes(group, frac_F1R2_alt, colour = group)) +
            geom_hline(yintercept = 0.5, linetype = 2) +
            geom_point(aes(y = frac_F1R2_ref), shape = 4, colour = "black") +
            geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi)) +
            geom_text(aes(y = 1.04, label = n_var), size = 2.5, colour = "black") +
            facet_wrap(~ m12, nrow = 2) + scale_colour_manual(values = GCOL) + coord_cartesian(ylim = c(0, 1.06)) +
            labs(x = NULL, y = "F1R2 share of ALT reads (pooled, 95% CI)", colour = NULL,
                 title = "ALT read orientation per class (REF-based), pooled per group",
                 subtitle = "row 2 = complement of row 1: one-strand damage skews the two in opposite directions\nx = REF reads at the same sites; number = variants") +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom"),
          "orientation_pooled", od, w = 11, h = 6.5)
save_plot(ggplot(d[tested == TRUE], aes(signed_orient, fill = group)) +
            geom_histogram(binwidth = 0.05, boundary = 0, position = "identity", alpha = 0.55) +
            geom_vline(xintercept = 0.5, linetype = 2) + facet_wrap(~ m12, nrow = 2, scales = "free_y") +
            scale_fill_manual(values = GCOL) +
            labs(x = "signed_orient = F1R2_alt / (F1R2_alt + F2R1_alt)", y = "variants", fill = NULL,
                 title = sprintf("Per-variant ALT orientation (variants with >= %d ALT orientation reads)", MIN_N),
                 subtitle = paste("tested:", sub_excl)) + theme(legend.position = "bottom"),
          "orientation_hist", od, w = 12, h = 6)
save_plot(ggplot(d[!is.na(signed_orient) & VAF > 0], aes(ROQ, VAF, colour = signed_orient)) +
            geom_point(size = 0.6, alpha = 0.7) + scale_y_log10() + facet_wrap(~ m12, nrow = 2) +
            scale_colour_gradient2(midpoint = 0.5, low = "#2166ac", mid = "grey85", high = "#b2182b", limits = c(0, 1)) +
            labs(colour = "signed\norient", title = "VAF vs ROQ, coloured by ALT orientation (PASS + orientation-filtered)"),
          "vaf_vs_roq", od, w = 12, h = 6)
d[, roq_class := factor(fifelse(mutation12 %chin% c("C>T", "G>A"), paste(mutation12, fifelse(CpG, "CpG", "non-CpG")), mutation12),
                        levels = c(paste(rep(c("C>T", "G>A"), each = 2), c("CpG", "non-CpG")), setdiff(M12, c("C>T", "G>A"))))]
d[, filt := factor(fifelse(group == G_OR, G_OR, "PASS"), levels = c("PASS", G_OR))]
save_plot(ggplot(d[!is.na(ROQ)], aes(roq_class, ROQ, fill = filt)) + geom_boxplot(outlier.size = 0.3) +
            geom_hline(yintercept = ROQ_CUT, linetype = 2) + scale_fill_manual(values = c(PASS = "grey70", setNames(GCOL[G_OR], G_OR))) +
            labs(x = NULL, fill = NULL, title = "ROQ per class (C>T / G>A split CpG vs non-CpG)") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom"),
          "roq_by_class", od, w = 10, h = 5)
fb <- rbind(d[tested == TRUE, .(pct = 100 * mean(orient_biased), n = .N, what = "% orient_biased (of tested)"), by = .(group, m12)],
            d[, .(pct = 100 * mean(group == G_OR), n = .N, group = "all", what = "% orientation-filtered (of PASS + orientation-filtered)"), by = m12])
save_plot(ggplot(fb, aes(m12, pct, fill = group)) + geom_col(position = position_dodge(preserve = "single")) +
            facet_wrap(~ what, ncol = 1, scales = "free_y") + scale_fill_manual(values = c(GCOL, all = "grey30")) +
            labs(x = NULL, y = "%", fill = NULL, title = "Orientation-biased and orientation-filtered calls per class") +
            theme(legend.position = "bottom"),
          "orient_biased_by_class", od, w = 10, h = 6)
save_plot(ggplot(d[mutation6 == "C>T" & VAF > 0], aes(group, VAF, fill = group)) + geom_violin(scale = "width") +
            geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) + scale_y_log10() +
            stat_summary(fun = median, geom = "text", aes(label = sprintf("%.3f", 10^after_stat(y))), vjust = -0.8, size = 3) +
            scale_fill_manual(values = GCOL, guide = "none") +
            labs(x = NULL, title = "VAF of C>T / G>A calls: orientation-filtered vs PASS by ROQ", subtitle = "number = median VAF"),
          "vaf_ct_by_group", od, w = 7, h = 4.5)

# ---- 7. signatures ----------------------------------------------------------------------------
cand <- d[!is.na(SBS96), .(CHROM, POS, REF, ALT1, SBS96, group = as.character(group), CALLER = "Mutect2 SNVs by group")]
sc <- data.table(group = character(0), reference = character(0), cosine = numeric(0))
if (!nrow(cand)) message("no SBS96 classes (no trinucleotide context) - signature step skipped") else {
plot_96context(cand, od, rowsplit = "group", name = "snv_96context_by_group", show_all = FALSE, roworder = GROUPS)
m96 <- dcast(cand, SBS96 ~ group, fun.aggregate = length, value.var = "POS")[data.table(SBS96 = SBS96_TYPES), on = "SBS96"]
for (g in setdiff(names(m96), "SBS96")) set(m96, which(is.na(m96[[g]])), g, 0L)
fwrite(setnames(copy(m96), "SBS96", "MutationType"), file.path(od, "sbs96_groups.txt"), sep = "\t")
refs <- data.table(SBS96 = SBS96_TYPES)
ffpe_f <- Filter(file.exists, c("ref/FFPE_signatures_Guo2022.tsv", "analysis/ref/FFPE_signatures_Guo2022.tsv",
  "scripts/analysis/ref/FFPE_signatures_Guo2022.tsv",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/ref/FFPE_signatures_Guo2022.tsv"))[1]
if (!is.na(ffpe_f)) refs <- fread(ffpe_f)[, .(SBS96 = MutationType, FFPE_repaired, FFPE_unrepaired)][refs, on = "SBS96"]
cosmic_f <- Sys.getenv("SPRTN_COSMIC_SBS", "")
if (!nzchar(cosmic_f)) cosmic_f <- sort(Sys.glob(file.path("/common/WORK/pstancl/envs/sigprofiler/lib/python3*/site-packages",
  "SigProfilerAssignment/data/Reference_Signatures/GRCh38/COSMIC_v3*_SBS_GRCh38.txt")), decreasing = TRUE)[1]
if (!is.na(cosmic_f) && file.exists(cosmic_f)) {
  cs <- fread(cosmic_f); setnames(cs, 1, "SBS96")
  refs <- cs[, intersect(c("SBS96", "SBS1", "SBS5", "SBS45"), names(cs)), with = FALSE][refs, on = "SBS96"]
  message("COSMIC: ", cosmic_f)
} else message("no COSMIC SBS file (set SPRTN_COSMIC_SBS) - SBS1 / SBS5 / SBS45 left out")
cosine <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))
grp_cols <- intersect(GROUPS, names(m96))
ref_cols <- c(setdiff(names(refs), "SBS96"), grp_cols)
allm <- refs[m96, on = "SBS96"]
sc <- CJ(group = grp_cols, reference = ref_cols, sorted = FALSE)[group != reference]
sc[, `:=`(n = vapply(group, function(g) sum(allm[[g]]), 0),
          cosine = round(mapply(function(g, r) cosine(allm[[g]], allm[[r]]), group, reference), 3))]
print(dcast(sc, group ~ reference, value.var = "cosine")); fwrite(sc, file.path(od, "signature_cosine.tsv"), sep = "\t")
save_plot(ggplot(sc, aes(reference, group, fill = cosine)) + geom_tile(colour = "white") +
            geom_text(aes(label = sprintf("%.2f", cosine)), size = 3.2) +
            scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
            labs(x = NULL, y = NULL, title = "Cosine similarity of each group's SBS96 profile",
                 subtitle = "SBS45 = 8-oxoG artefact; FFPE = Guo et al. 2022 (deamination-type)") +
            theme(axis.text.x = element_text(angle = 30, hjust = 1)),
          "signature_cosine", od, w = 8, h = 3.5)
}

# sites for the read-level step (13_damage_metrics.sh): C>T/G>A and G>T/C>A, up to SITES_MAX per group
set.seed(1)
st <- d[mutation6 %chin% c("C>T", "C>A")][, .SD[sample(.N, min(.N, SITES_MAX))], by = group]
fwrite(st[, .(chrom = CHROM, pos = POS, ref = REF, alt = ALT1, group)], file.path(od, "read_level", "sites.tsv"), sep = "\t")

# ---- 4. read level (after qsub_damage.sh) -----------------------------------------------------
rf <- file.path(od, "read_level", paste0("reads_", TUMOUR, ".tsv.gz"))
rl_sum <- NULL
if (file.exists(rf)) {
  r <- fread(rf)[allele %chin% c("REF", "ALT")]
  r[, pair := fifelse(paste0(ref, ">", alt) %chin% c("C>T", "G>A"), "C>T / G>A", "G>T / C>A")]
  r[, group := factor(group, levels = GROUPS)]
  rl_sum <- r[, .(reads = .N, median_dist_frag_end = as.numeric(median(dist_frag_end, na.rm = TRUE)),
                  pct_within_10bp = round(100 * mean(dist_frag_end <= 10, na.rm = TRUE), 1),
                  median_cycle = as.numeric(median(cycle)), median_bq = as.numeric(median(bq)),
                  median_mq = as.numeric(median(mq)), pct_R1 = round(100 * mean(mate == "R1"), 1)),
              by = .(group, pair, allele)][order(group, pair, allele)]
  wt <- r[, .(p_dist_ALT_vs_REF = if (uniqueN(allele) == 2) wilcox.test(dist_frag_end ~ allele)$p.value else NA_real_), by = .(group, pair)]
  rl_sum <- wt[rl_sum, on = c("group", "pair")]
  print(rl_sum); fwrite(rl_sum, file.path(od, "read_level_summary.tsv"), sep = "\t")
  save_plot(ggplot(r[!is.na(dist_frag_end) & dist_frag_end <= 300], aes(dist_frag_end, colour = allele)) +
              stat_ecdf() + facet_grid(pair ~ group) + scale_colour_manual(values = c(ALT = "#d62728", REF = "grey40")) +
              labs(x = "bases to the nearer fragment end", y = "cumulative share of reads", colour = NULL,
                   title = "Where the ALT base sits in the fragment: ALT vs REF reads at the same sites",
                   subtitle = "deamination during library prep -> ALT closer to fragment ends (weak expectation in fresh-frozen DNA)"),
            "read_level", od, w = 10, h = 6)
} else message("no ", rf, " - read-level step not run yet (cd ../wgs && qsub qsub_damage.sh)")

# ---- 6. tumour vs normal (after qsub_damage.sh) ----------------------------------------------------
sm <- file.path(od, "sample_metrics")
ss <- rbindlist(lapply(c(TUMOUR, NORMAL), function(s) {
  pf <- file.path(sm, paste0(s, ".pre_adapter_summary_metrics")); gf <- file.path(sm, paste0(s, ".giv.tsv"))
  if (!file.exists(pf) && !file.exists(gf)) return(NULL)
  out <- data.table(sample = s)
  if (file.exists(pf)) {
    p <- fread(pf, skip = "SAMPLE_ALIAS", fill = TRUE)[!is.na(TOTAL_QSCORE)]
    p <- p[, .(Q = min(TOTAL_QSCORE)), by = .(sub = paste0(REF_BASE, "_", ALT_BASE))]   # worst library
    for (i in seq_len(nrow(p))) set(out, j = paste0("Q_", p$sub[i]), value = p$Q[i])
  }
  if (file.exists(gf)) { g <- fread(gf); for (i in seq_len(nrow(g))) set(out, j = paste0("GIV_", g$type[i]), value = g$GIV[i]) }
  out
}), fill = TRUE)
TN_DIFF <- NA; tn_txt <- "Picard / GIV sample metrics not available yet (qsub_damage.sh)."
if (nrow(ss)) {
  setnames(ss, c("Q_C_T", "Q_G_T"), c("Q_Deamination_C_T", "Q_OxoG_G_T"), skip_absent = TRUE)
  fwrite(ss, file.path(od, "sample_summary.tsv"), sep = "\t"); print(ss)
  if (nrow(ss) == 2) {
    # our rule for "clearly different": pre-adapter Q differs by >= 5 (~3x error rate), or GIV on
    # opposite sides of 1.5, for deamination (C>T) or oxidation (G>T)
    dq <- function(col) if (col %in% names(ss)) abs(diff(ss[[col]])) >= 5 else FALSE
    dg <- function(col) if (col %in% names(ss)) length(unique(ss[[col]] > 1.5)) == 2 else FALSE
    hits <- c(Deamination_Q = dq("Q_Deamination_C_T"), OxoG_Q = dq("Q_OxoG_G_T"), GIV_C_T = dg("GIV_C_T"), GIV_G_T = dg("GIV_G_T"))
    TN_DIFF <- any(hits)
    tn_txt <- if (TN_DIFF) paste("Tumour and normal differ clearly in:", paste(names(hits)[hits], collapse = ", "),
                                 "(rule: Q differs by >= 5, or GIV on opposite sides of 1.5).")
              else "Tumour and normal do not differ clearly in Picard pre-adapter Q or GIV (rule: Q differs by >= 5, or GIV on opposite sides of 1.5)."
  }
}

# ---- report --------------------------------------------------------------------------------------
call_group <- function(g) {
  e <- ex[group == g]; m <- mir[group == g]
  up <- function(cls) nrow(e[mutation6 == cls & pct > pct_clean & q < 0.05]) > 0
  mr <- function(pr) isTRUE(m[pair == pr, mirrored])
  deam <- up("C>T") && mr("C>T | G>A"); oxo <- up("C>A") && mr("G>T | C>A")
  others <- sum(m[!pair %chin% c("C>T | G>A", "G>T | C>A"), mirrored], na.rm = TRUE)
  general <- others >= 3
  kind <- fcase(general, "general library / orientation problem (mirrored skew in most classes)",
                deam & oxo, "deamination-type and oxidation-type orientation artefact",
                deam, "deamination-type orientation artefact", oxo, "oxidation-type orientation artefact",
                up("C>T") || up("C>A"), "class excess without mirrored orientation skew - not explained by one-strand damage",
                default = "no meaningful damage signal")
  if (grepl("^(deamination|oxidation)", kind))
    kind <- paste0(kind, if (isTRUE(TN_DIFF)) "; tumour and normal differ in the read-level damage metrics (see below)"
                         else "; storage vs library prep not separable with one sample")
  list(kind = kind, deam = deam, oxo = oxo, others = others)
}
fmt_pair <- function(g, pr) { x <- mir[group == g & pair == pr]
  if (!nrow(x)) "n/a" else sprintf("F1R2 share %.2f vs %.2f (n %d / %d), q = %s%s", x$frac_1, x$frac_2, x$n_1, x$n_2,
                                   format(x$q, digits = 2), if (x$mirrored) ", mirrored" else "") }
med_vaf <- d[, .(v = median(VAF, na.rm = TRUE)), by = group]
rep_lines <- c(
  sprintf("# Damage QC of Mutect2 low-ROQ calls - %s (%s)", TUMOUR, format(Sys.Date())), "",
  "One tumour and one normal: storage time and batch effects cannot be modelled. Assumptions and open points are listed at the end.", "",
  "## ROQ", "",
  if (ROQ_OK) "ROQ present and FilterMutectCalls used orientation priors (header_check.txt)." else
    "**ROQ cannot be trusted** (no ROQ or no --ob-priors in the header) - groups by ROQ are descriptive only.", "",
  sprintf("Groups: %s", paste(excl[, sprintf("%s %d SNVs (%d with >= %d ALT orientation reads)", group, n, tested, MIN_N)], collapse = "; ")), "",
  "## Verdict per group", "")
for (g in intersect(c(G_LO, G_OR), as.character(excl$group))) {
  cg <- call_group(g)
  rep_lines <- c(rep_lines, sprintf("**%s: %s.**", g, cg$kind), "",
    sprintf("- C>T/G>A: %s%% of SNVs (clean group %s%%); orientation %s",
            ex[group == g & mutation6 == "C>T", pct], ex[group == g & mutation6 == "C>T", pct_clean], fmt_pair(g, "C>T | G>A")),
    sprintf("- G>T/C>A: %s%% (clean %s%%); orientation %s",
            ex[group == g & mutation6 == "C>A", pct], ex[group == g & mutation6 == "C>A", pct_clean], fmt_pair(g, "G>T | C>A")),
    sprintf("- other mirrored pairs: %d of 4", cg$others),
    sprintf("- median VAF %.3f (clean group %.3f)", med_vaf[group == g, v], med_vaf[group == G_HI, v]),
    sprintf("- C>T/G>A orientation at CpG vs non-CpG: %s",
            paste(ocg[group == g, sprintf("%s %s %.2f (n %d)", mutation12, CpG_class, frac_F1R2_alt, n_var)], collapse = "; ")),
    sprintf("- cosine to: %s", paste(sc[group == g, sprintf("%s %.2f", reference, cosine)], collapse = ", ")), "")
}
hi_n <- ch[orient_n >= 6, .(n = .N, obs = 100 * mean(one_sided), exp = 100 * mean(p_exp)), by = group]
rep_lines <- c(rep_lines, "## Chance one-sidedness (chance_one_sided_*.tsv, chance_one_sided.pdf)", "",
  "Expected P(all ALT reads on one side) = 2 x 0.5^n. Both groups are selected on one-sidedness, so some excess over chance is built in; n >= 6 is the telling part (chance <= 3%).", "",
  paste0("- ", ch_cls[mutation12 == "all", sprintf("%s, n >= 2: %d calls, %d one-sided vs %.1f expected (obs/exp %.2f)", group, n_var, observed, expected, obs_exp)]),
  paste0("- ", hi_n[, sprintf("%s, n >= 6: %d calls, %.1f%% one-sided vs %.1f%% expected", group, n, obs, exp)]), "",
  "## T>C / A>G by context (tc_ag_by_context.tsv)", "",
  if ("n_filtered" %in% names(tc)) paste0("- ", head(tc[!is.na(n_filtered)][order(-n_filtered)], 6)[
    , sprintf("%s %s: %d filtered calls (%s%% of the class; clean %s%%), F1R2 share %s", mutation12, context, n_filtered,
              pct_of_class_filtered, fifelse(is.na(pct_of_class_clean), "0", as.character(pct_of_class_clean)), frac_F1R2_alt_filtered)]),
  "", "Flagged as unexplained (possibly adenine deamination: hypoxanthine pairs with C, giving A>G / T>C) unless a context stands out.", "",
  "## Minimum ALT reads vs ROQ < 30, PASS SNVs (alt_filter_vs_roq.tsv; nothing applied)", "",
  paste0("- ", af[, sprintf("%s: removes %d (%.1f%% of PASS; %d confirmed by another caller; %.1f%% C>T); overlap with ROQ < 30: %d, only this rule: %d, only ROQ < 30: %d",
                            rule, removed, pct_of_pass, removed_confirmed, pct_CT, also_ROQ30, only_this, only_ROQ30)]), "")
rep_lines <- c(rep_lines, "## Tumour vs normal", "", tn_txt, "",
  if (!is.null(rl_sum)) c("## Read level (tumour)", "",
    paste0("- ", rl_sum[, sprintf("%s, %s, %s reads: n %d, median distance to fragment end %s, within 10 bp %s%%, median BQ %s, R1 %s%% (ALT vs REF distance p = %s)",
                                  group, pair, allele, reads, median_dist_frag_end, pct_within_10bp, median_bq, pct_R1,
                                  format(p_dist_ALT_vs_REF, digits = 2))]), "") else c("## Read level", "", "Not run yet.", ""),
  "## Assumptions and caveats", "",
  "- Orientation from Mutect2's F1R2/F2R1 (reads per allele by pair orientation); classes REF-based on the + strand, so a one-strand lesion gives opposite skews in a class and its complement.",
  sprintf("- Low-ROQ calls have few ALT reads, so the pooled test per class is the main test; the per-variant binomial covers only calls with >= %d ALT orientation reads. Pooled read counts are overdispersed (reads of one variant are not independent) - the sign test over variants (orientation_by_class.tsv, q_sign) does not assume that.", MIN_N),
  sprintf("- Mirrored skew = Fisher test of the two classes' F1R2/F2R1 ALT counts, q < 0.05 (BH over all pairs and groups), skews on opposite sides of 0.5 and F1R2 shares >= %.1f apart (pooled reads make tiny skews significant; the threshold is ours).", MIRROR_MIN),
  "- C>T at CpG alone is what real SBS1 (clock) mutations look like; it counts as damage only together with an orientation skew. Whether the library was PCR-free is not known; cluster amplification copies damaged templates either way.",
  "- orientation-filtered includes calls that also failed other filters (orientation_only flags those that failed only orientation).",
  "- Picard / GIV on a region slice (sample_metrics/run_info.tsv), duplicates excluded; GIV = R1 / R2 mismatch rate per type in read orientation, > 1.5 = damaged (Chen et al. 2017). The tumour-normal rule (Q differs by >= 5, GIV across 1.5) is ours.",
  "- Read position near fragment ends is a weak expectation for fresh-frozen DNA; a negative result there is uninformative.",
  "- No sample metadata (storage time, extraction / library batch, kit, polymerase) available.")
writeLines(rep_lines, file.path(od, "report.md"))
message("done: ", od, "  (report.md)")
