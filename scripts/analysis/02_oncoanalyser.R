# ---------------------------------------------------------------------------
# 02_oncoanalyser.R - oncoanalyser (Hartwig WiGiTS) output for RJALS, tumour vs normal.
#   SNV / indel: SAGE + PAVE + PURPLE    SV: ESVEE + LINX    copy number + purity: PURPLE
# Run:  Rscript 02_oncoanalyser.R   (or line by line in R)
# ---------------------------------------------------------------------------
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
od <- file.path(OUT, "oncoanalyser"); dir.create(od, recursive = TRUE, showWarnings = FALSE)

# ---- 1. input files ---------------------------------------------------------
pd <- file.path(ONCO, "purple")
files <- c(
  purple_somatic = find_one(pd, "\\.purple\\.somatic\\.vcf\\.gz$", required = FALSE),
  esvee          = find_one(file.path(ONCO, "esvee"), "\\.esvee\\.somatic\\.vcf\\.gz$", required = FALSE),
  linx_svs       = find_one(file.path(ONCO, "linx"),  "\\.linx\\.svs\\.tsv$",           required = FALSE),
  purple_cnv     = find_one(pd, "\\.purple\\.cnv\\.somatic\\.tsv$"),
  purple_pur     = find_one(pd, "\\.purple\\.purity\\.tsv$")
)
files <- files[!is.na(files)]
writeLines(paste(names(files), files, sep = "\t"), file.path(od, "inputs_used.tsv"))

# ---- 2. SNV / indel: PURPLE-annotated SAGE calls -----------------------------
calls <- load_calls(files, "purple_somatic", "sage_purple")
if (!is.null(calls)) calls_pass <- snv_indel_summary(calls, od)

# ---- 3. structural variants: ESVEE and LINX ------------------------------------
if ("esvee" %in% names(files)) {
  sv <- read_sv_vcf(files[["esvee"]])[, CALLER := "esvee"]
  sv_summary(sv, od)
}
# TODO: match breakpoints between ESVEE and Manta (StructuralVariantAnnotation) in a comparison script.
if ("linx_svs" %in% names(files)) {
  linx <- fread(files[["linx_svs"]])
  if ("type" %in% names(linx)) print(linx[, .N, by = type])    # LINX's classification of each SV
}

# ---- 4. copy number and purity: PURPLE ---------------------------------------
seg <- purple_segments(files[["purple_cnv"]])
cn_summary(seg, od)
pp <- fread(files[["purple_pur"]])
purity <- data.table(field = names(pp), value = unlist(lapply(pp[1], as.character)))   # transposed: one row per field
print(purity); fwrite(purity, file.path(od, "purple_purity.csv"))

# ---- 5. qcVCF -----------------------------------------------------------------
# TODO(Paula): call your qcVCF functions here and tell me what they take.
if (HAVE_QCVCF) {
  library(qcVCF)
  # qc <- qcVCF::<function>(files[["purple_somatic"]])
}
message("done: ", od)
