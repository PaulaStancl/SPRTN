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
  sage_somatic   = find_one(file.path(ONCO, "sage", "somatic"), "\\.sage\\.somatic\\.vcf\\.gz$", required = FALSE),
  purple_somatic = find_one(pd, "\\.purple\\.somatic\\.vcf\\.gz$", required = FALSE),
  hmf_sigs       = find_one(file.path(ONCO, "sigs"), "\\.sig\\.allocation\\.tsv$", required = FALSE),
  esvee          = find_one(file.path(ONCO, "esvee"), "\\.esvee\\.somatic\\.vcf\\.gz$", required = FALSE),
  linx_svs       = find_one(file.path(ONCO, "linx"),  "\\.linx\\.svs\\.tsv$",           required = FALSE),
  purple_cnv     = find_one(pd, "\\.purple\\.cnv\\.somatic\\.tsv$"),
  purple_pur     = find_one(pd, "\\.purple\\.purity\\.tsv$")
)
files <- files[!is.na(files)]
writeLines(paste(names(files), files, sep = "\t"), file.path(od, "inputs_used.tsv"))

# ---- 2. SNV / indel: SAGE (+ PAVE, PURPLE) ---------------------------------------
# qc_vcf/raw:  every SAGE record, PASS and filtered (sage/somatic - PURPLE's VCF is read for
#              PASS only) - FILTER outcome, reasons, depth / VAF / SAGE's QUAL
# qc_vcf/pass: PASS calls from PURPLE's somatic VCF (SAGE's calls with PAVE + PURPLE annotation,
#              e.g. info_PURPLE_AF purity-adjusted VAF, info_PURPLE_CN, info_SUBCL) - VAF,
#              spectrum, qcVCF plots (section 5), MNVs split as for sarek (overlaps);
#              96-context and spectrum leave doublets out (a doublet is one DBS event)
raw_dir  <- file.path(od, "qc_vcf", "raw");  dir.create(raw_dir,  recursive = TRUE, showWarnings = FALSE)
pass_dir <- file.path(od, "qc_vcf", "pass"); dir.create(pass_dir, recursive = TRUE, showWarnings = FALSE)
calls_raw <- load_calls(files, if ("sage_somatic" %in% names(files)) "sage_somatic" else "purple_somatic", "sage")
calls     <- load_calls(files, "purple_somatic", "sage")
if (!is.null(calls_raw)) raw_qc(calls_raw, raw_dir)
if (!is.null(calls)) {
  pass <- snv_indel_summary(calls, pass_dir)
  calls_pass <- pass$pass      # PASS records as SAGE wrote them (MNVs whole)
  atom       <- pass$atom      # the same with MNVs split into SNVs
  # trinucleotide context from the GATK GRCh38 FASTA - same coordinates as Hartwig's GRCh38;
  # add_sbs96() warns if a REF does not match
  add_sbs96(atom)
}

# ---- 3. structural variants: ESVEE and LINX -> sv/ ------------------------------------
sv_dir <- file.path(od, "sv"); dir.create(sv_dir, showWarnings = FALSE)
if ("esvee" %in% names(files)) {
  sv <- read_sv_vcf(files[["esvee"]])[, CALLER := "esvee"]
  sv_summary(sv, sv_dir)
}
# Manta vs ESVEE breakpoints: 04_summary_table.R (within 100 bp)
if ("linx_svs" %in% names(files)) {
  linx <- fread(files[["linx_svs"]])
  if ("type" %in% names(linx)) {                                 # LINX's classification of each SV
    lt <- linx[, .N, by = type][order(-N)]
    print(lt); fwrite(lt, file.path(sv_dir, "linx_sv_types.csv"))
  }
}

# ---- 4. copy number and purity: PURPLE -> cnv/ ---------------------------------------
cn_dir <- file.path(od, "cnv"); dir.create(cn_dir, showWarnings = FALSE)
seg <- purple_segments(files[["purple_cnv"]])
cn_summary(seg, cn_dir)
pp <- fread(files[["purple_pur"]])
purity <- data.table(field = names(pp), value = unlist(lapply(pp[1], as.character)))   # transposed: one row per field
print(purity); fwrite(purity, file.path(cn_dir, "purple_purity.csv"))
# outputs of the earlier flat layout, now in sv/ and cnv/
unlink(file.path(od, c("sv_counts.csv", "sv_types.pdf", "cn_segments.csv", "cn_genome.pdf",
                       "cn_state_fraction.csv", "cn_loh_fraction.csv", "purple_purity.csv")))

# ---- 5. qcVCF on the PASS calls ----------------------------------------------
# One caller here, so no pairwise sharing / shared-unique split (that is the sarek comparison).
if (HAVE_QCVCF && exists("calls_pass")) {
  library(qcVCF)
  save_plot(plot_mutation_counts(calls_pass[, .(CHROM, POS, REF, ALT = ALT1, subjectID = PATIENT, tool = CALLER,
                                                mutType = MUTTYPE)], grp_col = "tool", cohort = PATIENT),
            "mutation_counts", pass_dir, w = 6, h = 6)
  plot_96context(atom, pass_dir)                                    # 96-context profile of the PASS SNVs
}

# ---- 6. mutational signatures: PASS calls for SigProfiler and FitMS ------------------
# Same as sarek's section 7, one set: sage (all PASS calls, MNVs split). Then
#   python 01b_sarek_signatures.py --pipeline oncoanalyser
#   Rscript --no-environ 01c_sarek_signatures_organ.R oncoanalyser
# Hartwig's own signature fit (sigs/, on PURPLE's SNVs) is copied next to them to compare.
if (exists("atom")) {
  sig_in <- file.path(od, "signatures", "input")
  write_signature_sets(rbindlist(list(sage = atom), idcol = "SET"), sig_in, "02_oncoanalyser.R")
  if ("hmf_sigs" %in% names(files))
    fwrite(fread(files[["hmf_sigs"]]), file.path(od, "signatures", "hartwig_sigs_allocation.csv"))
  message("signature inputs: ", sig_in,
          "  ->  python 01b_sarek_signatures.py --pipeline oncoanalyser; Rscript 01c_sarek_signatures_organ.R oncoanalyser")
}
message("done: ", od)
