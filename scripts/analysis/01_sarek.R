# ---------------------------------------------------------------------------
# 01_sarek.R - sarek output for RJALS, tumour vs normal.
#   SNV / indel: Mutect2 + Strelka2     SV: Manta     copy number + purity: ASCAT
# Run:  Rscript 01_sarek.R   (or line by line in R - any of the folders above works)
# ---------------------------------------------------------------------------
# Finds 00_setup.R from this folder, from scripts/, or from the project root.
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
od <- file.path(OUT, "sarek"); dir.create(od, recursive = TRUE, showWarnings = FALSE)

# ---- 1. input files ---------------------------------------------------------
vc <- file.path(SAREK, "variant_calling")
files <- c(
  mutect2       = find_one(file.path(vc, "mutect2", PAIR), "\\.mutect2\\.filtered\\.vcf\\.gz$"),
  strelka_snv   = find_one(file.path(vc, "strelka", PAIR), "somatic_snvs\\.vcf\\.gz$",   required = FALSE),
  strelka_indel = find_one(file.path(vc, "strelka", PAIR), "somatic_indels\\.vcf\\.gz$", required = FALSE),
  manta         = find_one(file.path(vc, "manta",   PAIR), "\\.manta\\.somatic_sv\\.vcf\\.gz$", required = FALSE),
  ascat_seg     = find_one(file.path(vc, "ascat",   PAIR), "\\.segments\\.txt$"),
  ascat_pp      = find_one(file.path(vc, "ascat",   PAIR), "\\.purityploidy\\.txt$")
)
files <- files[!is.na(files)]
writeLines(paste(names(files), files, sep = "\t"), file.path(od, "inputs_used.tsv"))

# ---- 2. SNV / indel: Mutect2 and Strelka2 -------------------------------------
calls <- rbindlist(list(
  load_calls(files, "mutect2", "mutect2"),
  load_calls(files, c("strelka_snv", "strelka_indel"), "strelka")      # NULL (skipped) if not found
), fill = TRUE)                                                         # the callers' metric columns differ
calls_pass <- snv_indel_summary(calls, od)

# ---- 3. structural variants: Manta ---------------------------------------------
if ("manta" %in% names(files)) {
  sv <- read_sv_vcf(files[["manta"]])[, caller := "manta"]
  sv_summary(sv, od)
}

# ---- 4. copy number and purity: ASCAT ------------------------------------------
seg <- ascat_segments(files[["ascat_seg"]])
cn_summary(seg, od)
purity <- fread(files[["ascat_pp"]])
print(purity); fwrite(purity, file.path(od, "ascat_purity_ploidy.csv"))

# ---- 5. qcVCF -----------------------------------------------------------------
# TODO(Paula): call your qcVCF functions here and tell me what they take (a VCF
# path? a data frame?) so I can wire them in.
if (HAVE_QCVCF) {
  library(qcVCF)
  # qc <- qcVCF::<function>(files[["mutect2"]])
}
message("done: ", od)
