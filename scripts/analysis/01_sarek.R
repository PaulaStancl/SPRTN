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
  sv <- read_sv_vcf(files[["manta"]])[, CALLER := "manta"]
  sv_summary(sv, od)
}

# ---- 4. copy number and purity: ASCAT ------------------------------------------
seg <- ascat_segments(files[["ascat_seg"]])
cn_summary(seg, od)
purity <- fread(files[["ascat_pp"]])
print(purity); fwrite(purity, file.path(od, "ascat_purity_ploidy.csv"))

# ---- 5. qcVCF on the PASS calls ----------------------------------------------
# qcVCF's CHROM / POS / REF / ALT match these tables; ALT is ALT1 (the first allele - PASS
# calls have one). mutType uses qcVCF's own palette classes (SNV, MNV, INS, DEL) rather than
# add_mutational_type(), which only knows SNV / INDEL and would call Mutect2's MNVs indels.
if (HAVE_QCVCF) {
  library(qcVCF)
  qc <- calls_pass[, .(CHROM, POS, REF, ALT = ALT1, subjectID = PATIENT, tool = CALLER,
                       mutType = fcase(TYPE == "SNV", "SNV", TYPE == "MNV", "MNV",
                                       TYPE == "INDEL" & nchar(ALT1) > nchar(REF), "INS",
                                       TYPE == "INDEL", "DEL", default = "INDEL"))]

  # PASS calls per caller and mutation type
  save_plot(plot_mutation_counts(qc, grp_col = "tool", cohort = PATIENT),
            "qcVCF_mutation_counts", od, w = 8, h = 6)

  # % of each caller's mutations also found by the other caller
  pw <- plot_pairwise_shared_mutations(qc, sample_col = "tool", cohort = PATIENT)
  fwrite(pw$data, file.path(od, "qcVCF_pairwise_shared.csv"))
  save_plot(pw$plot, "qcVCF_pairwise_shared", od, w = 7, h = 6)

  # Per mutation: shared (found by >1 caller) or unique. mode "any" - "within_subject"
  # calls stringr::str_remove_all(), which qcVCF does not import, so it fails without stringr.
  ov <- overlap_and_annotate_shared_mutations(qc, sample_col = "tool", mode = "any")
  shared <- unique(ov[, .(CHROM, POS, REF, ALT, tool = tool.x, shared_status)])
  calls_pass[shared, on = .(CHROM, POS, REF, ALT1 = ALT, CALLER = tool), QC_SHARED := i.shared_status]
  print(calls_pass[, .N, by = .(CALLER, TYPE, QC_SHARED)])
  fwrite(calls_pass, file.path(od, "snv_indel_pass.csv"))      # now with QC_SHARED

  # Not used: plot96_matrix() needs each SNV's trinucleotide context (NC_3), which sarek's
  # VCFs do not carry - SAGE's TNC field has it, so it fits 02_oncoanalyser.R instead.
}
message("done: ", od)
