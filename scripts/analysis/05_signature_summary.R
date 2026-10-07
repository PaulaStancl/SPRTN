# ---------------------------------------------------------------------------
# 05_signature_summary.R - mutational signature attributions of several PASS sets side by side:
# SigProfilerAssignment (COSMIC, 01b) and FitMS (liver common + rare, 01c), one fit per set.
#   Rscript 05_signature_summary.R comparison    # muse, mutect2, strelka, sage, >= 2, all callers (03)
#   Rscript 05_signature_summary.R sarek         # mutect2, strelka, mutect2_strelka (01)
# Run after 01b (--pipeline <same>) and, for the FitMS part, 01c (<same>).
# Output: <OUT>/<pipeline>/signatures/summary/
#   attribution_<method>.pdf   stacked bars: fraction of each set's mutations per signature,
#                              signatures < 5% in every set grouped as "other"; n mutations on top.
#                              Colours as SigProfiler's activity plots (artefact signatures grey).
#                              FitMS twice: liver organ signatures + unassigned, and as RefSig.
#                              FitMS's own per-set plots and JSON: ../fitms_liver/<set>/ (01c)
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
# FitMS: liver signatures as fitted (incl. "unassigned" - mutations no signature explains)
# and converted to reference signatures (RefSig, comparable with COSMIC)
f <- file.path(sig, "fitms_liver", "exposures_organ.csv")
if (file.exists(f)) att[["FitMS organ"]] <- fread(f)[, method := "FitMS liver (organ signatures + unassigned)"]
f <- file.path(sig, "fitms_liver", "exposures_refsig.csv")
if (file.exists(f)) att[["FitMS"]] <- fread(f)[, method := "FitMS liver (as reference signatures)"]
if (!any(grepl("^FitMS", names(att))))
  message("no FitMS results under ", file.path(sig, "fitms_liver"), " - run 01c_sarek_signatures_organ.R ", PIPELINE)
f <- file.path(sig, "fitms_liver", "fit_summary.csv")
if (file.exists(f)) { x <- fread(f)
  # one FitMS fit per set, on the liver organ signatures; the RefSig view is a conversion of the
  # same exposures, so this cosine (catalogue vs reconstruction) belongs to both FitMS plots
  fq[["FitMS"]] <- x[, .(method = "FitMS liver (organ signatures; RefSig view = same fit)", set, n = snvs, cosine = cosine_similarity,
                         unassigned_pct, rare_signature)] }
if (!length(att)) stop("no signature fits under ", sig, " - run 01b_sarek_signatures.py --pipeline ", PIPELINE,
                       " (and 01c_sarek_signatures_organ.R ", PIPELINE, ") first", call. = FALSE)
att <- rbindlist(att, use.names = TRUE, fill = TRUE)[mutations > 0]
att[, fraction := mutations / sum(mutations), by = .(method, set)]       # recompute: shares of the set's fit
fwrite(att[order(method, set, -mutations)], file.path(od, "attribution_all.csv"))

# ---- colours: SigProfiler's own (sigProfilerPlotting plotActivity.py, v1.4.3) ---------------
# fixed colours for common signatures; COSMIC's artefact signatures in greys; every other
# signature takes the next colour of SigProfiler's list (matplotlib names, here as hex) in
# signature order. Returned in SigProfiler's stacking order (fixed, others, artefacts), then
# "other" (light grey, lighter than the artefact greys) and FitMS's "unassigned" (near white).
sp_fixed <- c(SBS1 = "#acf2d0", SBS5 = "#63d69e", SBS2 = "#f8b6b3", SBS13 = "#f17fb2", SBS3 = "#c4abc4",
  SBS4 = "#bcf2f5", SBS7a = "#b5d7f5", SBS7b = "#9ecef7", SBS7c = "#84bdf0", SBS7d = "#6cb2f0",
  SBS8 = "#dfc4f5", SBS9 = "#ebf5bc", SBS10a = "#f2aeae", SBS10b = "#f08080", SBS17a = "#d9f7b0",
  SBS17b = "#8cc63f", SBS40 = "#c4c4f5", SBS6 = "#faf1dc", SBS14 = "#faecca", SBS15 = "#fcebc2",
  SBS20 = "#fae4af", SBS21 = "#fae1a5", SBS26 = "#fcde97", SBS44 = "#fad682")
sp_fixed <- c(sp_fixed, SBS40a = "#c4c4f5", SBS40b = "#a8a8eb", SBS40c = "#8c8ce0")  # COSMIC v3.4 split SBS40
sp_artefact <- c(SBS27 = "#C8C8C8", SBS43 = "#C0C0C0", SBS45 = "#BEBEBE", SBS46 = "#B8B8B8", SBS47 = "#B0B0B0",
  SBS48 = "#A9A9A9", SBS49 = "#A8A8A8", SBS50 = "#A0A0A0", SBS51 = "#989898", SBS52 = "#909090",
  SBS53 = "#888888", SBS54 = "#808080", SBS55 = "#787878", SBS56 = "#707070", SBS57 = "#696969",
  SBS58 = "#686868", SBS59 = "#606060", SBS60 = "#585858")
sp_list <- c("#e377c2", "#ff7f0e", "#9467bd", "#bcbd22", "#8c564b", "#d62728", "#2ca02c", "#17becf",  # tab:*
  "#ff1493", "#ff4500", "#8a2be2", "#d2691e", "#006400", "#1e90ff", "#c71585", "#fa8072", "#ff00ff",
  "#f4a460", "#228b22", "#4169e1", "#da70d6", "#4b0082", "#8fbc8f", "#0000ff", "#db7093", "#483d8b",
  "#6b8e23", "#00ffff", "#ff69b4", "#663399", "#00ff00")
OTHER_LAB <- "other (< 5% in every set)"; OTHER_COL <- "#E3E3E3"
sig_colours <- function(sigs) {
  sigs <- setdiff(sigs, c(OTHER_LAB, "unassigned"))
  ref  <- sub("^.*_(common|rare)_", "", sigs)               # FitMS GEL-Liver_common_SBS1 -> SBS1
  num  <- suppressWarnings(as.numeric(sub(".*?([0-9]+)[^0-9]*$", "\\1", ref)))  # signature order: SBS2 < SBS10a
  fx   <- sigs[match(intersect(names(sp_fixed), ref), ref)]
  af   <- sigs[match(intersect(names(sp_artefact), ref), ref)]
  rest <- setdiff(sigs, c(fx, af)); rest <- rest[order(num[match(rest, sigs)], rest)]
  c(setNames(sp_fixed[ref[match(fx, sigs)]], fx), setNames(rep_len(sp_list, length(rest)), rest),
    setNames(sp_artefact[ref[match(af, sigs)]], af), setNames(OTHER_COL, OTHER_LAB), unassigned = "#F5F5F5")
}

# ---- stacked bars per method ------------------------------------------------------------
for (m in unique(att$method)) {
  x <- att[method == m]
  keep <- x[, .(mx = max(fraction)), by = signature][mx >= 0.05, signature]
  x[, sig_lab := fifelse(signature %chin% keep | signature == "unassigned", signature, OTHER_LAB)]
  x <- x[, .(fraction = sum(fraction), mutations = sum(mutations)), by = .(set, sig_lab)]
  cols <- sig_colours(unique(x$sig_lab))                       # SigProfiler order: SBS1 at the bottom
  x[, sig_lab := factor(sig_lab, levels = names(cols))]
  x[, SET := lab_set(set)]
  tot <- x[, .(n = sum(mutations)), by = SET]
  save_plot(ggplot(x, aes(SET, fraction, fill = sig_lab)) + geom_col(width = 0.8, colour = "white", linewidth = 0.3,
                                                                   position = position_stack(reverse = TRUE)) +
              scale_fill_manual(values = cols) +
              geom_text(data = tot, aes(SET, 1.02, label = paste0("n=", format(round(n), big.mark = ",", trim = TRUE))),
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
                   title = "Signature fit quality per PASS set",
                   subtitle = "dashed: 0.9, a good fit; FitMS: one fit per set (organ signatures), shown in 05 also as RefSig") +
              theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom"),
            "fit_quality", od, w = 3 + 1.1 * uniqueN(fq$SET), h = 5)
}
message("done: ", od)
