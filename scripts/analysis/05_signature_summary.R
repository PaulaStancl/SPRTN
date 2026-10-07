# ---------------------------------------------------------------------------
# 05_signature_summary.R - mutational signature attributions of several PASS sets side by side:
# SigProfilerAssignment (COSMIC, 01b) and FitMS (liver common + rare, 01c), one fit per set.
#   Rscript 05_signature_summary.R comparison    # muse, mutect2, strelka, sage, >= 2, all callers (03)
#   Rscript 05_signature_summary.R sarek         # mutect2, strelka, mutect2_strelka (01)
# Run after 01b (--pipeline <same>) and, for the FitMS part, 01c (<same>).
# Output: <OUT>/<pipeline>/signatures/summary/
#   attribution_<method>.pdf   stacked bars: fraction of each set's mutations per signature,
#                              signatures < 5% in every set grouped as "other"; n mutations on top
#   attribution_heatmap.pdf    set x signature, % of the set's mutations (all methods)
#   fit_quality.pdf / .csv     cosine similarity of each fit (and FitMS unassigned %, rare signature)
#   attribution_all.csv        everything in one long table: method, set, signature, mutations, fraction
# ---------------------------------------------------------------------------
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
PIPELINE <- commandArgs(TRUE)[1]; if (is.na(PIPELINE)) PIPELINE <- "comparison"
sig <- file.path(OUT, PIPELINE, "signatures")
od  <- file.path(sig, "summary"); dir.create(od, recursive = TRUE, showWarnings = FALSE)

# set order and display names
set_order <- c("muse", "mutect2", "strelka", "sage", "mutect2_strelka", "two_plus", "all_callers")
set_name  <- c(muse = "MuSE", mutect2 = "Mutect2", strelka = "Strelka2", sage = "SAGE/PURPLE",
               mutect2_strelka = "Mutect2 + Strelka2", two_plus = ">= 2 callers", all_callers = "all callers")
lab_set <- function(x) factor(ifelse(x %in% names(set_name), set_name[x], x),
                              levels = unique(c(set_name[set_order], sort(unique(x)))))

# ---- read the fits -----------------------------------------------------------------
att <- list(); fq <- list()
for (ctx in c("SBS96", "DBS78", "ID83")) {
  f <- file.path(sig, "sigprofiler", paste0("activities_", ctx, ".csv"))
  if (file.exists(f)) att[[paste("SigProfiler", ctx)]] <- fread(f)[, method := paste("SigProfiler COSMIC", ctx)]
  f <- file.path(sig, "sigprofiler", paste0("fit_stats_", ctx, ".csv"))
  if (file.exists(f)) { x <- fread(f)
    fq[[paste("SigProfiler", ctx)]] <- data.table(method = paste("SigProfiler COSMIC", ctx), set = x[[1]],
      n = x[["Total Mutations"]], cosine = x[["Cosine Similarity"]]) }
}
f <- file.path(sig, "fitms_liver", "exposures_refsig.csv")
if (file.exists(f)) att[["FitMS"]] <- fread(f)[, method := "FitMS liver (as reference signatures)"]
f <- file.path(sig, "fitms_liver", "fit_summary.csv")
if (file.exists(f)) { x <- fread(f)
  fq[["FitMS"]] <- x[, .(method = "FitMS liver (as reference signatures)", set, n = snvs, cosine = cosine_similarity,
                         unassigned_pct, rare_signature)] }
if (!length(att)) stop("no signature fits under ", sig, " - run 01b_sarek_signatures.py --pipeline ", PIPELINE,
                       " (and 01c_sarek_signatures_organ.R ", PIPELINE, ") first", call. = FALSE)
att <- rbindlist(att, use.names = TRUE, fill = TRUE)[mutations > 0]
att[, fraction := mutations / sum(mutations), by = .(method, set)]       # recompute: shares of the set's fit
fwrite(att[order(method, set, -mutations)], file.path(od, "attribution_all.csv"))

# ---- stacked bars per method ------------------------------------------------------------
for (m in unique(att$method)) {
  x <- att[method == m]
  keep <- x[, .(mx = max(fraction)), by = signature][mx >= 0.05, signature]
  x[, sig_lab := fifelse(signature %chin% keep, signature, "other (< 5% in every set)")]
  x <- x[, .(fraction = sum(fraction), mutations = sum(mutations)), by = .(set, sig_lab)]
  ord <- x[sig_lab != "other (< 5% in every set)", .(t = sum(fraction)), by = sig_lab][order(-t), sig_lab]
  x[, sig_lab := factor(sig_lab, levels = c(ord, "other (< 5% in every set)"))]
  x[, SET := lab_set(set)]
  tot <- x[, .(n = sum(mutations)), by = SET]
  save_plot(ggplot(x, aes(SET, fraction, fill = sig_lab)) + geom_col(width = 0.8) +
              geom_text(data = tot, aes(SET, 1.02, label = paste0("n=", format(round(n), big.mark = ","))),
                        inherit.aes = FALSE, size = 3, vjust = 0) +
              scale_y_continuous(labels = function(v) paste0(100 * v, "%"), expand = expansion(mult = c(0, 0.08))) +
              labs(x = NULL, y = "share of the set's mutations", fill = "signature",
                   title = paste("Signature attribution per PASS set -", m),
                   subtitle = "one fit per set; n = mutations attributed") +
              theme(axis.text.x = element_text(angle = 30, hjust = 1)),
            paste0("attribution_", sub("_$", "", gsub("[^A-Za-z0-9]+", "_", m))), od, w = 3 + 1.1 * uniqueN(x$SET), h = 6)
}

# ---- heatmap: all methods, % of each set's mutations ------------------------------------------
h <- att[, .(method, SET = lab_set(set), signature, pct = 100 * fraction)]
h <- h[h[, .(mx = max(pct)), by = .(method, signature)][mx >= 5], on = .(method, signature)]
save_plot(ggplot(h, aes(SET, signature, fill = pct)) + geom_tile(colour = "white") +
            geom_text(aes(label = sprintf("%.0f", pct)), size = 2.8) +
            facet_grid(method ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(18)) +
            scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 100)) +
            labs(x = NULL, y = NULL, fill = "% of set", title = "Signature attribution per PASS set (all methods)",
                 subtitle = "% of each set's attributed mutations; signatures >= 5% in at least one set") +
            theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text.y = element_text(angle = 0)),
          "attribution_heatmap", od, w = 3 + 1.1 * uniqueN(h$SET), h = 2 + 0.28 * nrow(unique(h[, .(method, signature)])))

# ---- fit quality ------------------------------------------------------------------------------
if (length(fq)) {
  fq <- rbindlist(fq, fill = TRUE)
  print(fq); fwrite(fq, file.path(od, "fit_quality.csv"))
  fq[, SET := lab_set(set)]
  save_plot(ggplot(fq, aes(SET, cosine, colour = method, group = method)) + geom_point(size = 2.5) + geom_line() +
              geom_hline(yintercept = 0.9, linetype = 2, colour = "grey50") +
              scale_y_continuous(limits = c(min(0.5, min(fq$cosine, na.rm = TRUE)), 1)) +
              labs(x = NULL, y = "cosine similarity (profile vs reconstruction)", colour = NULL,
                   title = "Signature fit quality per PASS set", subtitle = "dashed: 0.9, a good fit") +
              theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom"),
            "fit_quality", od, w = 3 + 1.1 * uniqueN(fq$SET), h = 5)
}
message("done: ", od)
