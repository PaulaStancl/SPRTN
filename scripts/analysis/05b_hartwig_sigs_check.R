# ---------------------------------------------------------------------------
# 05b_hartwig_sigs_check.R - SAGE/PURPLE set: Hartwig SIGS (oncoanalyser) vs our SigProfiler
# and FitMS fits.
#   Rscript 05b_hartwig_sigs_check.R
# Run after 01b_sarek_signatures.py --pipeline oncoanalyser (and 01c ... oncoanalyser for FitMS).
# Expected: identical SNV counts and contexts - both leave doublets (MNVs) out of SBS96.
#
# SIGS (hmftools): PURPLE's PASS single-base SNVs (MNVs skipped), context from PURPLE's
# annotation, least-squares fit of the 30 COSMIC v2 signatures (Sig1-Sig30); the unexplained
# part is reported as UNALLOC / MISALLOC.
#
# Output <OUT>/oncoanalyser/signatures/hartwig_check/
#   sbs96_check.csv        96 channels: SNVs in SIGS and in SigProfiler's matrix, difference
#   sbs96_check.pdf        the same, channel by channel (on the dashed line = identical)
#   signatures_compare.csv / .pdf   % of the SNVs per signature: SIGS, SigProfiler, FitMS
# ---------------------------------------------------------------------------
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
sig <- file.path(OUT, "oncoanalyser", "signatures")
od  <- file.path(sig, "hartwig_check"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
f_mg <- file.path(sig, "sigprofiler", "matrix_generator", "output", "SBS", paste0(PATIENT, ".SBS96.all"))
f_sp <- file.path(sig, "sigprofiler", "activities_SBS96.csv")
f_fm <- file.path(sig, "fitms_liver", "exposures_refsig.csv")
if (!file.exists(f_mg)) stop("missing ", f_mg, " - run 01b_sarek_signatures.py --pipeline oncoanalyser first", call. = FALSE)

# ---- 1. same SNVs and contexts? -------------------------------------------------------------
# SIGS bucket C>A_ACA = SigProfiler A[C>A]A
hmf <- fread(find_one(file.path(ONCO, "sigs"), "\\.sig\\.snv_counts\\.csv$"))
setnames(hmf, 1:2, c("bucket", "SIGS"))
hmf[, MutationType := paste0(substr(bucket, 5, 5), "[", substr(bucket, 1, 3), "]", substr(bucket, 7, 7))]
mg <- fread(f_mg)
chk <- merge(hmf[, .(MutationType, SIGS)], mg[, .(MutationType = get(names(mg)[1]), SigProfiler = sage)],
             by = "MutationType", all = TRUE)
chk[, diff := SigProfiler - SIGS]
fwrite(chk, file.path(od, "sbs96_check.csv"))
cos <- with(chk, sum(SIGS * SigProfiler) / sqrt(sum(SIGS^2) * sum(SigProfiler^2)))
txt <- sprintf("SNVs: SIGS %d, SigProfiler %d | channels identical: %d / 96 | cosine similarity %.4f",
               sum(chk$SIGS), sum(chk$SigProfiler), sum(chk$diff == 0), cos)
message(txt)
if (any(chk$diff != 0)) print(chk[diff != 0])
save_plot(ggplot(chk, aes(SIGS, SigProfiler, colour = substr(MutationType, 3, 5))) +
            geom_abline(linetype = 2, colour = "grey50") + geom_point(size = 2) +
            labs(colour = NULL, x = "Hartwig SIGS: SNVs per 96 channel", y = "SigProfiler matrix: SNVs per channel",
                 title = "SAGE/PURPLE PASS SNVs: SIGS vs SigProfiler", subtitle = txt),
          "sbs96_check", od, w = 7, h = 5.5)

# ---- 2. signatures side by side -----------------------------------------------------------
# SIGS uses COSMIC v2 (Sig1-30), SigProfiler COSMIC v3, FitMS reference signatures (RefSig):
# rows matched on the signature number (Sig3 ~ SBS3 ~ RefSig 3). v3 split some v2 signatures
# (Sig7 -> SBS7a-d, Sig17 -> SBS17a/b ...) and added many (SBS31 and up), which v2 cannot fit.
num_of <- function(s) fifelse(grepl("[0-9]", s), sub("^[^0-9]*([0-9]+).*$", "\\1", s), s)
ha <- fread(find_one(file.path(ONCO, "sigs"), "\\.sig\\.allocation\\.tsv$")); setnames(ha, tolower(names(ha)))
res <- list(`Hartwig SIGS (COSMIC v2)` = ha[, .(signature, mutations = allocation)])
if (file.exists(f_sp)) res[["SigProfiler (COSMIC v3)"]] <- fread(f_sp)[set == "sage", .(signature, mutations)]
if (file.exists(f_fm)) res[["FitMS liver (RefSig)"]] <- fread(f_fm)[set == "sage", .(signature, mutations)] else
  message("no FitMS result (", f_fm, ") - run 01c_sarek_signatures_organ.R oncoanalyser to add it")
res <- rbindlist(res, idcol = "method")
res[, `:=`(num = num_of(signature), pct = 100 * mutations / sum(mutations)), by = method]
res[, label := fifelse(grepl("^[0-9]+$", num), paste0("Sig/SBS ", num), num)]
wide <- dcast(res, label ~ method, value.var = "pct", fun.aggregate = function(v) round(sum(v), 1), fill = 0)
wide[res[, .(names = paste(unique(signature), collapse = " | ")), by = label], on = "label", names := i.names]
wide <- wide[order(-do.call(pmax, wide[, setdiff(names(wide), c("label", "names")), with = FALSE]))]
print(wide); fwrite(wide, file.path(od, "signatures_compare.csv"))

keep <- res[, .(mx = max(pct)), by = label][mx >= 2, label]
pl <- res[label %chin% keep, .(pct = sum(pct)), by = .(method, label)]
pl[, method := factor(method, levels = rev(unique(res$method)))]
pl[, label := factor(label, levels = pl[, .(m = max(pct)), by = label][order(m), label])]
save_plot(ggplot(pl, aes(pct, label, fill = method)) +
            geom_col(position = position_dodge(preserve = "single"), width = 0.8) +
            scale_fill_manual(values = c(`Hartwig SIGS (COSMIC v2)` = "#7f7f7f", `SigProfiler (COSMIC v3)` = "#63d69e",
                                         `FitMS liver (RefSig)` = "#9467bd"), breaks = unique(res$method)) +
            labs(x = "% of the SNVs", y = NULL, fill = NULL,
                 title = "SAGE/PURPLE set: signatures by method",
                 subtitle = "matched on signature number; signatures >= 2% in any method") +
            theme(legend.position = "bottom"),
          "signatures_compare", od, w = 8, h = 2.5 + 0.35 * length(keep))
message("done: ", od)
