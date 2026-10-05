# ---------------------------------------------------------------------------
# 01c_sarek_signatures_organ.R - liver-specific common + rare SBS signatures in the
# sarek PASS sets, with FitMS (signature.tools.lib; Degasperi et al. 2022, Science).
# Run after 01b_sarek_signatures.py, in the signature-tools env:
#
#   micromamba activate /common/WORK/pstancl/envs/signature-tools
#   Rscript 01c_sarek_signatures_organ.R
#
# Input  <OUT>/sarek/signatures/sigprofiler/matrix_generator/output/SBS/RJALS.SBS96.all
#        - SigProfilerMatrixGenerator's SBS96 matrix (01b); one column per set: mutect2,
#        strelka, mutect2_strelka.
# Output <OUT>/sarek/signatures/fitms_liver/
#   <set>/                         plotFitMS() plots + FitMS result as JSON, one folder per set
#   exposures_organ.csv            set, signature (GEL-Liver_common_* or a rare one), mutations, fraction
#   exposures_refsig.csv           the same in reference signatures (organ-specific -> RefSig)
#   fit_summary.csv                per set: mutations, cosine similarity, unassigned %, the
#                                  rare signature chosen (if any) and the other candidates
#
# FitMS fits the liver *common* signatures first, then tests the *rare* signatures one at a
# time and keeps one only if it reduces the fit error enough (errorReduction, >= 15%).
# Each set is fitted on its own (one FitMS call per set, as in 01b). Settings follow the
# package README / help: common tier T1 (organ-specific), rare tier T2 (advised), KLD,
# bootstrap with the Gini-scaled exposure filter.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(signature.tools.lib) })

PATIENT <- "RJALS"
ORGAN   <- "Liver"
RESULTS <- Sys.getenv("SPRTN_RESULTS", "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/results/wgs")
OUT     <- Sys.getenv("SPRTN_OUT", file.path(RESULTS, "analysis"))
if (grepl("OneDrive|CloudStorage|Dropbox|iCloud|Google Drive", normalizePath(OUT, mustWork = FALSE), ignore.case = TRUE))
  stop("Refusing to use a cloud-synced folder for patient data: ", OUT, call. = FALSE)
NCPU    <- max(1L, as.integer(Sys.getenv("NCPUS", "4")))
NBOOT   <- as.integer(Sys.getenv("SPRTN_NBOOT", "200"))              # package default
sig_dir <- file.path(OUT, "sarek", "signatures")
mat_f   <- file.path(sig_dir, "sigprofiler", "matrix_generator", "output", "SBS", paste0(PATIENT, ".SBS96.all"))
od      <- file.path(sig_dir, "fitms_liver"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(mat_f)) stop("missing ", mat_f, " - run 01b_sarek_signatures.py first", call. = FALSE)

# SigProfiler's SBS96 matrix -> catalogue matrix (channels as rows, sets as columns), rows in
# the order of the liver signatures (SigProfiler sorts A[C>A]A, A[C>A]C, A[C>G]A ...;
# signature.tools.lib groups by substitution: A[C>A]A ... T[C>A]T, A[C>G]A ...)
m   <- fread(mat_f)
cat <- as.matrix(m[, -1]); rownames(cat) <- m[[1]]
sigs <- getSignaturesForFitting(organ = ORGAN, typemut = "subs", commontier = "T1", raretier = "T2", verbose = FALSE)
stopifnot(setequal(rownames(cat), rownames(sigs$common)))
cat <- cat[rownames(sigs$common), , drop = FALSE]
message("liver signatures: ", ncol(sigs$common), " common (", paste(colnames(sigs$common), collapse = ", "),
        "), ", ncol(sigs$rare), " rare candidates")
fwrite(data.table(signature = c(colnames(sigs$common), colnames(sigs$rare)),
                  type = rep(c("common", "rare"), c(ncol(sigs$common), ncol(sigs$rare)))),
       file.path(od, "liver_signatures_used.csv"))

exp_org <- list(); exp_ref <- list(); summ <- list()
for (st in colnames(cat)) {
  x <- cat[, st, drop = FALSE]
  if (sum(x) == 0) { message(st, ": no SNVs - skipped"); next }
  message("FitMS ", st, " (", sum(x), " SNVs)")
  res <- FitMS(catalogues = x, organ = ORGAN, commonSignatureTier = "T1", rareSignatureTier = "T2",
               exposureFilterType = "giniScaledThreshold", useBootstrap = TRUE, nboot = NBOOT,
               nparallel = NCPU, randomSeed = 1, verbose = FALSE)
  sd <- file.path(od, st); unlink(sd, recursive = TRUE); dir.create(sd)
  plotFitMS(res, outdir = paste0(sd, "/"))
  writeFitResultsToJSON(fitObj = res, filename = file.path(sd, paste0(st, "_FitMS.json")))   # gzipped

  e  <- res$exposures[st, ]                                       # signatures + "unassigned"
  exp_org[[st]] <- data.table(set = st, signature = names(e), mutations = as.numeric(e),
                              fraction = as.numeric(e) / sum(x))[mutations > 0]
  # organ-specific -> reference signatures (rare signatures are already reference ones)
  org <- setdiff(names(e)[e > 0 & names(e) %in% colnames(sigs$common)], "unassigned")
  conv <- if (length(org)) convertExposuresFromOrganToRefSigs(expMatrix = matrix(e[org], ncol = 1, dimnames = list(org, st)),
                                                              typemut = "subs") else NULL
  rest <- setdiff(names(e)[e > 0], c(org, "unassigned"))
  ref  <- c(if (!is.null(conv)) setNames(conv[, 1], rownames(conv)), setNames(e[rest], rest))
  ref  <- tapply(ref, names(ref), sum)
  exp_ref[[st]] <- data.table(set = st, signature = names(ref), mutations = as.numeric(ref),
                              fraction = as.numeric(ref) / sum(x))[mutations > 0]
  summ[[st]] <- data.table(set = st, snvs = sum(x),
                           cosine_similarity = round(as.numeric(res$cossim_catalogueVSreconstructed[1]), 3),
                           unassigned_pct = round(100 * e[["unassigned"]] / sum(x), 1),
                           rare_signature = if (!is.null(res$rareSigChoice[[st]])) res$rareSigChoice[[st]] else "none",
                           rare_candidates = paste(res$candidateRareSigs[[st]], collapse = ";"))
}
exp_org <- rbindlist(exp_org)[order(set, -mutations)]
exp_ref <- rbindlist(exp_ref)[order(set, -mutations)]
summ    <- rbindlist(summ)
fwrite(exp_org, file.path(od, "exposures_organ.csv"))
fwrite(exp_ref, file.path(od, "exposures_refsig.csv"))
fwrite(summ,    file.path(od, "fit_summary.csv"))
print(summ); print(exp_ref)
message("done: ", od)
