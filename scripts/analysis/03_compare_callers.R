# ---------------------------------------------------------------------------
# 03_compare_callers.R - PASS SNVs / indels of every somatic caller that ran, compared:
#   sarek: Mutect2, Strelka2 (main run) + MuSE and any other caller of the extra run
#          (sarek/RJALS_vc, SAREK_STEP=variant_calling)
#   oncoanalyser: SAGE, as in PURPLE's final somatic VCF
# A caller that has not run yet is simply absent; re-run this once it has.
# Run:  Rscript 03_compare_callers.R   (or line by line in R)
#
# Before it: scripts/wgs/06_normalize_vcfs.sh, so every caller's indels are left-aligned the
# same way (otherwise identical indels written differently count as different mutations).
# MNVs are split into SNVs for every comparison (Strelka and MuSE call none), mutations are
# matched on CHROM:POS:REF:ALT. The pipelines aligned separately (sarek: GATK GRCh38 +
# BQSR; oncoanalyser: Hartwig's GRCh38), so some disagreement at low VAF is expected.
# Output: <OUT>/comparison/
#   inputs_used.csv                which VCF per caller, normalised or not
#   pass_counts.csv/.pdf           PASS calls per caller and mutation type
#   snv_indel_concordance.csv/.pdf every combination of callers (MNVs split)
#   mutation_callers.csv           one row per PASS mutation: callers, n callers, each one's VAF
#   n_callers_per_caller.csv/.pdf  of each caller's calls, how many callers found them
#   pairwise_overlap.csv, pairwise_overlap_<SNV|INDEL>.pdf   % of row caller's calls in column caller
#   support_by_vaf.csv/.pdf        share of each caller's calls confirmed by >= 1 other, by VAF
#   vaf_agreement.csv/.pdf         tumour VAF of shared SNVs, caller vs caller (n, Pearson, Spearman)
#   consensus_summary.csv          union, >= 2 callers, all callers - per type
#   concordance_groups.csv/.pdf    exact combinations + summary bars two_plus and all_callers
#   snv_96context_groups.pdf/.csv  96-context profile of two_plus / all_callers SNVs (callers named)
#   snv_96context_per_caller.pdf/.csv  per caller: all / shared / unique SNVs
#   sbs96_groups.txt               the same counts, SigProfiler matrix format
#   snv_vaf_hist.pdf, snv_spectrum.pdf, snv_substitution_spectrum.csv, snv_indel_pass*.csv
# ---------------------------------------------------------------------------
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
od <- file.path(OUT, "comparison"); dir.create(od, recursive = TRUE, showWarnings = FALSE)

# ---- 1. input files ---------------------------------------------------------
vcfs <- find_caller_vcfs()
print(vcfs); fwrite(vcfs, file.path(od, "inputs_used.csv"))
if (uniqueN(vcfs$CALLER) < 2) stop("fewer than two callers found - nothing to compare", call. = FALSE)
if (!all(vcfs$normalised)) message("not normalised (indels may not match): ",
  paste(unique(vcfs[normalised == FALSE, CALLER]), collapse = ", "), " - run scripts/wgs/06_normalize_vcfs.sh")
calls <- rbindlist(lapply(split(vcfs, by = "CALLER"), function(x)
  rbindlist(lapply(x$file, read_vcf_table), fill = TRUE)[, CALLER := x$CALLER[1]]), fill = TRUE)

# ---- 2. PASS calls: counts, concordance of all callers, VAF, spectrum --------------
pass <- snv_indel_summary(calls, od)        # writes the concordance, VAF histogram, spectrum, PASS tables
calls_pass <- pass$pass; atom <- pass$atom
cnt <- calls_pass[, .N, by = .(CALLER, MUTTYPE)][order(CALLER, MUTTYPE)]
fwrite(dcast(cnt, CALLER ~ MUTTYPE, value.var = "N", fill = 0), file.path(od, "pass_counts.csv"))
save_plot(ggplot(cnt, aes(CALLER, N, fill = MUTTYPE)) + geom_col(position = position_dodge(0.9)) +
            geom_text(aes(label = N), position = position_dodge(0.9), vjust = -0.3, size = 3) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            labs(x = NULL, y = "PASS calls (as called)", title = "PASS calls per caller"),
          "pass_counts", od, w = 8, h = 5)

# ---- 3. per mutation: which callers found it -----------------------------------
atom[, KEY := paste(CHROM, POS, REF, ALT1, sep = ":")]
atom <- unique(atom, by = c("KEY", "CALLER"))
mut <- atom[, .(CALLERS = paste(sort(CALLER), collapse = "+"), N_CALLERS = .N),
            by = .(KEY, CHROM, POS, REF, ALT = ALT1, TYPE)]
vaf_w <- dcast(atom, KEY ~ CALLER, value.var = "VAF")
setnames(vaf_w, setdiff(names(vaf_w), "KEY"), paste0("VAF_", setdiff(names(vaf_w), "KEY")))
mut <- vaf_w[mut, on = "KEY"][order(match(CHROM, STD_CHR), POS)]
fwrite(mut, file.path(od, "mutation_callers.csv"))
atom[mut, on = "KEY", N_CALLERS := i.N_CALLERS]

# callers that report this mutation type at all (MuSE: SNVs only) - "all callers" means these
able <- atom[, .(callers = list(sort(unique(CALLER)))), by = TYPE]
cons <- mut[, .(union = .N, two_plus = sum(N_CALLERS >= 2)), by = TYPE]
cons[able, on = "TYPE", `:=`(n_callers = lengths(i.callers), callers = vapply(i.callers, paste, "", collapse = "+"))]
mut[cons, on = "TYPE", n_able := i.n_callers]
cons[mut[, .(all_callers = sum(N_CALLERS == n_able)), by = TYPE], on = "TYPE", all_callers := i.all_callers]
mut[, n_able := NULL]
# labels for plots: "sage" is SAGE's calls as in PURPLE's final VCF
disp <- function(x) gsub("\\bsage\\b", "sage/purple", x)
cons[, `:=`(two_plus_label    = paste0("two_plus: >= 2 of ", disp(gsub("+", ", ", callers, fixed = TRUE))),
            all_callers_label = paste0("all_callers: ", disp(callers)))]
print(cons); fwrite(cons, file.path(od, "consensus_summary.csv"))

# concordance bars: every exact combination of callers (each mutation once, in the set of
# callers that found it) plus two summary groups - found by >= 2 callers (two_plus) and by
# every caller of that type (all_callers); the summary bars overlap the combination bars
cg <- rbind(mut[, .(n = .N), by = .(TYPE, group = disp(CALLERS))][, kind := "exact combination"],
            cons[, .(TYPE, group = two_plus_label, n = two_plus, kind = "summary")],
            cons[, .(TYPE, group = all_callers_label, n = all_callers, kind = "summary")])
fwrite(cg[order(TYPE, kind, -n)], file.path(od, "concordance_groups.csv"))
# summary bars on top, combinations below them by size
# (ordered within each panel: the label carries the type, stripped again on the axis)
setorder(cg, TYPE, -kind, n)
cg[, group_f := factor(paste0(group, "___", TYPE), levels = paste0(group, "___", TYPE))]
save_plot(ggplot(cg, aes(group_f, n, fill = kind)) + geom_col() +
            geom_text(aes(label = n), hjust = -0.15, size = 3.2) + coord_flip() +
            facet_wrap(~TYPE, scales = "free") + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
            scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
            scale_fill_manual(values = c("exact combination" = "grey55", summary = "steelblue")) +
            labs(x = NULL, y = "PASS mutations (MNVs split)", fill = NULL,
                 title = "Concordance: exact caller combinations and summary groups"),
          "concordance_groups", od, w = 11, h = 6)

# of each caller's calls, how many callers found them (1 = this caller only)
nc <- atom[, .N, by = .(CALLER, TYPE, N_CALLERS)][, pct := round(100 * N / sum(N), 1), by = .(CALLER, TYPE)]
fwrite(nc[order(CALLER, TYPE, N_CALLERS)], file.path(od, "n_callers_per_caller.csv"))
save_plot(ggplot(nc, aes(CALLER, N, fill = factor(N_CALLERS))) + geom_col(position = "fill") +
            geom_text(aes(label = N), position = position_fill(vjust = 0.5), size = 3) +
            facet_wrap(~TYPE, scales = "free_x") + scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
            labs(x = NULL, y = "share of the caller's PASS calls", fill = "found by\nn callers",
                 title = "How many callers found each caller's PASS calls"),
          "n_callers_per_caller", od, w = 9, h = 5)

# ---- 4. pairwise overlap: % of the row caller's calls also PASS in the column caller --------
pw <- rbindlist(lapply(c("SNV", "INDEL"), function(ty) {
  a <- atom[TYPE == ty, .(KEY, CALLER)]
  cl <- sort(unique(a$CALLER))
  if (length(cl) < 2) return(NULL)
  rbindlist(lapply(cl, function(x) rbindlist(lapply(cl, function(y) {
    kx <- a[CALLER == x, KEY]
    data.table(TYPE = ty, caller = x, in_caller = y, n = length(kx), shared = sum(kx %chin% a[CALLER == y, KEY]))
  }))))
}))
pw[, pct := round(100 * shared / n, 1)]
fwrite(pw, file.path(od, "pairwise_overlap.csv"))
for (ty in unique(pw$TYPE))
  save_plot(ggplot(pw[TYPE == ty], aes(in_caller, caller, fill = pct)) + geom_tile(colour = "white") +
              geom_text(aes(label = sprintf("%.1f%%\n%d", pct, shared)), size = 3.2) +
              scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100)) +
              labs(x = "... also PASS in", y = "PASS calls of", fill = "%",
                   title = paste0(ty, ": % of each caller's calls found by the other")),
            paste0("pairwise_overlap_", ty), od, w = 6.5, h = 5)

# ---- 5. confirmation by another caller, by VAF ---------------------------------------
# among callers able to call the type; low-VAF calls are where callers disagree most
atom[able, on = "TYPE", n_able := lengths(i.callers)]
sv <- atom[n_able >= 2 & !is.na(VAF)]
sv[, VAF_BIN := cut(VAF, c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 1), include.lowest = TRUE)]
sv <- sv[, .(n = .N, confirmed_pct = round(100 * mean(N_CALLERS >= 2), 1)), by = .(TYPE, CALLER, VAF_BIN)][order(TYPE, CALLER, VAF_BIN)]
fwrite(sv, file.path(od, "support_by_vaf.csv"))
save_plot(ggplot(sv, aes(VAF_BIN, confirmed_pct, colour = CALLER, group = CALLER)) + geom_line() + geom_point(aes(size = n)) +
            facet_wrap(~TYPE) + scale_y_continuous(limits = c(0, 100)) + scale_size_area(max_size = 4) +
            labs(x = "tumour VAF (the caller's own)", y = "% also PASS in >= 1 other caller", size = "calls",
                 title = "Caller agreement by VAF"),
          "support_by_vaf", od, w = 10, h = 5)

# ---- 6. VAF of shared SNVs, caller vs caller -----------------------------------------
vc <- grep("^VAF_", names(mut), value = TRUE)
if (length(vc) >= 2) {
  prs <- combn(vc, 2, simplify = FALSE)
  va <- rbindlist(lapply(prs, function(p) mut[TYPE == "SNV" & !is.na(get(p[1])) & !is.na(get(p[2])),
    .(pair = paste(sub("VAF_", "", p), collapse = " vs "), x = get(p[1]), y = get(p[2]))]))
  if (nrow(va)) {
    agr <- va[, .(n = .N, pearson = round(cor(x, y), 3), spearman = round(cor(x, y, method = "spearman"), 3),
                  median_diff = round(median(y - x), 4)), by = pair]
    print(agr); fwrite(agr, file.path(od, "vaf_agreement.csv"))
    agr[, label := sprintf("n = %d\nPearson r = %.2f\nSpearman rho = %.2f", n, pearson, spearman)]
    save_plot(ggplot(va, aes(x, y)) + geom_point(size = 0.3, alpha = 0.3) + geom_abline(colour = "red", linetype = 2) +
                geom_text(data = agr, aes(x = 0.02, y = 0.98, label = label), hjust = 0, vjust = 1, size = 3,
                          inherit.aes = FALSE) +
                facet_wrap(~pair) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
                labs(x = "tumour VAF, first caller", y = "tumour VAF, second caller",
                     title = "Tumour VAF of SNVs found by both callers", subtitle = "red: y = x"),
              "vaf_agreement", od, w = 9, h = 7)
  }
}
# ---- 7. 96-context profiles: two_plus and all_callers SNVs --------------------------------
# Each SNV once per group (not once per caller); all_callers is a subset of two_plus. The plot
# names the callers of each group; sbs96_groups.txt (SigProfiler matrix format) uses the short
# names two_plus / all_callers. Context from the GATK GRCh38 FASTA (add_sbs96(), as in 01/02).
unlink(file.path(od, c("snv_96context.pdf", "snv_96context.csv")))     # the earlier all/two_plus/all_callers version
snv <- mut[TYPE == "SNV", .(CHROM, POS, REF, ALT1 = ALT, MUTTYPE = "SNV", N_CALLERS)]
add_sbs96(snv)
sc  <- cons[TYPE == "SNV"]
grp <- rbindlist(list(
  two_plus    = snv[N_CALLERS >= 2],
  all_callers = snv[N_CALLERS == sc$n_callers]), idcol = "GROUP")
grp[, GROUP := factor(GROUP, levels = c("two_plus", "all_callers"))]
print(grp[, .(snvs = .N, with_context = sum(!is.na(SBS96))), by = GROUP])
m96 <- dcast(grp[!is.na(SBS96)], SBS96 ~ GROUP, fun.aggregate = length, value.var = "POS", drop = FALSE)
m96 <- m96[data.table(SBS96 = SBS96_TYPES), on = "SBS96"]
for (cl in setdiff(names(m96), "SBS96")) set(m96, which(is.na(m96[[cl]])), cl, 0L)
setnames(m96, "SBS96", "MutationType")
fwrite(m96, file.path(od, "sbs96_groups.txt"), sep = "\t")
# plot_96context() draws one column per CALLER: here the group, labelled with its callers
grp[, CALLER := fifelse(GROUP == "two_plus", sc$two_plus_label, sc$all_callers_label)]
plot_96context(grp, od, name = "snv_96context_groups")  # snv_96context_groups.pdf / .csv

# per caller, as in 01: one column per caller, rows all SNVs / shared (PASS in >= 1 other
# caller) / unique (this caller only)
add_sbs96(atom)
atom[, QC_SHARED := fifelse(N_CALLERS >= 2, "shared", "unique")]
pc <- atom[, .(CHROM, POS, REF, ALT1, SBS96, QC_SHARED, CALLER = disp(CALLER))]
plot_96context(pc, od, rowsplit = "QC_SHARED", name = "snv_96context_per_caller")
message("done: ", od)
