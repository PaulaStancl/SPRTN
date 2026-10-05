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
# SNV / indel VCFs left-aligned by scripts/wgs/06_normalize_vcfs.sh (bcftools norm), used
# instead of sarek's raw ones once they exist - otherwise identical indels written at
# different places in a repeat count as two calls.
nd <- file.path(SAREK, "normalized_bcftools", PAIR)
norm <- if (dir.exists(nd)) c(
  mutect2       = find_one(nd, "\\.mutect2\\.filtered\\.norm\\.vcf\\.gz$", required = FALSE),
  strelka_snv   = find_one(nd, "somatic_snvs\\.norm\\.vcf\\.gz$",           required = FALSE),
  strelka_indel = find_one(nd, "somatic_indels\\.norm\\.vcf\\.gz$",         required = FALSE))
if (is.null(norm)) message("no normalised VCFs (", nd, ") - using sarek's raw ones; run scripts/wgs/06_normalize_vcfs.sh")
files <- c(
  mutect2       = find_one(file.path(vc, "mutect2", PAIR), "\\.mutect2\\.filtered\\.vcf\\.gz$"),
  strelka_snv   = find_one(file.path(vc, "strelka", PAIR), "somatic_snvs\\.vcf\\.gz$",   required = FALSE),
  strelka_indel = find_one(file.path(vc, "strelka", PAIR), "somatic_indels\\.vcf\\.gz$", required = FALSE),
  manta         = find_one(file.path(vc, "manta",   PAIR), "\\.manta\\.somatic_sv\\.vcf\\.gz$", required = FALSE),
  ascat_seg     = find_one(file.path(vc, "ascat",   PAIR), "\\.segments\\.txt$"),
  ascat_pp      = find_one(file.path(vc, "ascat",   PAIR), "\\.purityploidy\\.txt$")
)
files <- files[!is.na(files)]
if (!is.null(norm)) files[names(norm)[!is.na(norm)]] <- norm[!is.na(norm)]
writeLines(paste(names(files), files, sep = "\t"), file.path(od, "inputs_used.tsv"))

# ---- 2. SNV / indel: Mutect2 and Strelka2 -------------------------------------
calls <- rbindlist(list(
  load_calls(files, "mutect2", "mutect2"),
  load_calls(files, c("strelka_snv", "strelka_indel"), "strelka")      # NULL (skipped) if not found
), fill = TRUE)                                                         # the callers' metric columns differ
# qc_vcf/raw:  every record, PASS and filtered - FILTER outcome, reasons, depth / VAF / score
# qc_vcf/pass: PASS calls only - caller overlap, VAF, spectrum, qcVCF plots (section 5)
raw_dir  <- file.path(od, "qc_vcf", "raw");  dir.create(raw_dir,  recursive = TRUE, showWarnings = FALSE)
pass_dir <- file.path(od, "qc_vcf", "pass"); dir.create(pass_dir, recursive = TRUE, showWarnings = FALSE)
raw_qc(calls, raw_dir)
pass <- snv_indel_summary(calls, pass_dir)
calls_pass <- pass$pass      # PASS records as the callers wrote them (Mutect2 MNVs whole)
atom       <- pass$atom      # the same with MNVs split into SNVs - for comparing callers

# trinucleotide context of every PASS SNV (split MNVs included) from the reference FASTA:
# NC_3 on the + strand, SBS96 the pyrimidine-strand class (A[C>T]G) used by 96-context
# plots and SigProfiler
atom[MUTTYPE == "SNV", NC_3 := trinuc_context(.SD)]
bad <- atom[!is.na(NC_3) & substr(NC_3, 2, 2) != REF, .N]      # middle base must be REF
if (bad) warning(bad, " SNVs whose REF is not the reference base - is FASTA the genome sarek used?")
atom[MUTTYPE == "SNV" & substr(NC_3, 2, 2) == REF, SBS96 := sbs96(REF, ALT1, NC_3)]

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
# calls have one). mutType is MUTTYPE (SNV, MNV, INS, DEL - qcVCF's own palette classes)
# rather than add_mutational_type(), which only knows SNV / INDEL.
# Counts use the records as called; the caller comparisons use `atom` (MNVs split into
# SNVs), since Strelka writes a CC>TT as two C>T and would otherwise never match Mutect2.
if (HAVE_QCVCF) {
  library(qcVCF)
  to_qc <- function(x) x[, .(CHROM, POS, REF, ALT = ALT1, subjectID = PATIENT, tool = CALLER, mutType = MUTTYPE)]

  # PASS calls per caller and mutation type
  save_plot(plot_mutation_counts(to_qc(calls_pass), grp_col = "tool", cohort = PATIENT),
            "mutation_counts", pass_dir, w = 8, h = 6)

  # % of each caller's mutations also found by the other caller, separately for SNVs (split
  # MNVs included) and indels (insertions + deletions).
  qc <- to_qc(atom)
  qc[, varClass := fifelse(mutType %chin% c("INS", "DEL"), "INDEL", mutType)]
  # qcVCF groups with by = get(...), which data.table rejects when there are only a handful
  # of mutations (<= ~3 per caller), so a tiny class is skipped rather than stopping the script.
  pw_all <- rbindlist(lapply(unique(qc$varClass), function(cl) {
    pw <- tryCatch(plot_pairwise_shared_mutations(qc[varClass == cl], sample_col = "tool",
                                                  cohort = paste(PATIENT, cl)),
                   error = function(e) { message("qcVCF pairwise ", cl, " skipped: ", conditionMessage(e)); NULL })
    if (is.null(pw)) return(NULL)
    save_plot(pw$plot, paste0("pairwise_shared_", cl), pass_dir, w = 7, h = 6)
    pw$data[, varClass := cl]
  }))
  print(pw_all); fwrite(pw_all, file.path(pass_dir, "pairwise_shared.csv"))

  # Per mutation: shared (found by >1 caller) or unique. mode "any" - "within_subject"
  # calls stringr::str_remove_all(), which qcVCF does not import, so it fails without stringr.
  ov <- overlap_and_annotate_shared_mutations(qc, sample_col = "tool", mode = "any")
  shared <- unique(ov[, .(CHROM, POS, REF, ALT, tool = tool.x, shared_status)])
  atom[shared, on = .(CHROM, POS, REF, ALT1 = ALT, CALLER = tool), QC_SHARED := i.shared_status]
  # back onto the records as called: an MNV is shared / unique if all its bases are,
  # "partial" if only some are
  calls_pass[atom[FROM_MNV == FALSE], on = .(CHROM, POS, REF, ALT1, CALLER), QC_SHARED := i.QC_SHARED]
  if (any(atom$FROM_MNV)) {
    mnv_sh <- atom[FROM_MNV == TRUE, .(QC_SHARED = if (uniqueN(QC_SHARED) == 1) QC_SHARED[1] else "partial"),
                   by = .(MNV_KEY, CALLER)]
    calls_pass[, MNV_KEY := fifelse(MUTTYPE == "MNV", paste0(CHROM, ":", POS, ":", REF, ">", ALT1), NA_character_)]
    calls_pass[mnv_sh, on = .(MNV_KEY, CALLER), QC_SHARED := i.QC_SHARED][, MNV_KEY := NULL]
  }
  print(calls_pass[, .N, by = .(CALLER, MUTTYPE, QC_SHARED)])
  fwrite(calls_pass, file.path(pass_dir, "snv_indel_pass.csv"))            # now with QC_SHARED
  fwrite(atom,       file.path(pass_dir, "snv_indel_pass_atomized.csv"))

  # 96-context profile of the PASS SNVs (split MNVs included), one column per caller, rows
  # all / shared / unique. plot96_matrix() wants NC_3 on the pyrimidine strand (as palimpsest
  # writes it) and does not flip it itself, so it gets SBS96's bases, not the + strand NC_3.
  # It returns the figure(s) without saving; orderplots / showperc / dropempty / dontshowall
  # must be given as single values (their defaults are vectors, which its if() checks reject).
  snv96 <- atom[!is.na(SBS96), .(CHROM, POS, REF, ALT = ALT1, tool = CALLER, QC_SHARED,
                                 NC_3 = paste0(substr(SBS96, 1, 1), substr(SBS96, 3, 3), substr(SBS96, 7, 7)))]
  # plot96_matrix() also needs these, which qcVCF does not declare - skip the plot, not the script
  miss96 <- Filter(function(pk) !requireNamespace(pk, quietly = TRUE), c("cowplot", "stringr", "stringi", "ggtext"))
  if (length(miss96)) message("96-context plot skipped - install: ", paste(miss96, collapse = ", "))
  if (nrow(snv96)) fwrite(snv96, file.path(pass_dir, "snv_96context.csv"))
  if (nrow(snv96) && !length(miss96)) {
    fig96 <- plot96_matrix(snv96, rowsplit = "QC_SHARED", plotsplitcol = "tool",
                           orderplots = "no", showperc = "yes", dropempty = "no", dontshowall = "no")
    for (k in seq_along(fig96))
      save_plot(fig96[[k]], paste0("snv_96context", if (k > 1) paste0("_", k)), pass_dir,
                w = 9 * uniqueN(snv96$tool), h = 2 + 1.7 * (uniqueN(snv96$QC_SHARED) + 1))
  }
}
# ---- 6. mutational signatures: PASS sets for SigProfiler ------------------------
# Each PASS set as a minimal VCF in signatures/input/vcf/ (one "sample" per file), for
# 01b_sarek_signatures.py: SigProfilerMatrixGenerator builds SBS96 / DBS78 / ID83 from them
# (the only matrices the fits use)
# and SigProfilerAssignment fits COSMIC signatures (cosmic_fit, which replaced
# SigProfilerSingleSample). MNVs are written split: the matrix generator rejoins adjacent
# SNVs into doublets itself, the same way for both callers. Sets (keys CHROM:POS:REF:ALT):
#   mutect2 / strelka                 all PASS calls of that caller
#   mutect2_strelka                   PASS in both callers (intersection)
sig_in <- file.path(od, "signatures", "input")
dir.create(file.path(sig_in, "vcf"), recursive = TRUE, showWarnings = FALSE)
unlink(list.files(file.path(sig_in, "vcf"), "\\.vcf$", full.names = TRUE))   # no sets left over from an earlier run
atom[, CALLERS := paste(sort(unique(CALLER)), collapse = "+"), by = .(CHROM, POS, REF, ALT1)]
sets <- rbindlist(list(
  mutect2         = atom[CALLER == "mutect2"],
  strelka         = atom[CALLER == "strelka"],
  mutect2_strelka = atom[CALLER == "mutect2" & CALLERS == "mutect2+strelka"]
), idcol = "SET")[, .(SET, CHROM, POS, REF, ALT = ALT1, MUTTYPE, FROM_MNV, CALLERS, NC_3, SBS96)]
sets <- unique(sets, by = c("SET", "CHROM", "POS", "REF", "ALT"))
sets <- sets[order(SET, match(CHROM, STD_CHR), POS)]
set_counts <- dcast(sets, SET ~ MUTTYPE, fun.aggregate = length, value.var = "POS")
print(set_counts); fwrite(set_counts, file.path(sig_in, "pass_set_counts.csv"))
fwrite(sets, file.path(sig_in, "pass_sets.csv"))
for (st in unique(sets$SET)) {
  f <- file.path(sig_in, "vcf", paste0(st, ".vcf"))
  writeLines(c("##fileformat=VCFv4.2", paste0("##source=01_sarek.R PASS set ", st),
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), f)
  fwrite(sets[SET == st, .(CHROM, POS, ID = ".", REF, ALT, QUAL = ".", FILTER = "PASS", INFO = ".")],
         f, sep = "\t", append = TRUE, col.names = FALSE)
}
# The same SBS96 counts computed here (+ strand context from the GATK FASTA -> sbs96()), in
# SigProfiler's matrix format - NOT used for fitting, only to compare with the matrix
# generator's (01b writes sbs96_compare.csv, class by class).
m96 <- dcast(sets[!is.na(SBS96)], SBS96 ~ SET, fun.aggregate = length, value.var = "POS")
m96 <- m96[data.table(SBS96 = SBS96_TYPES), on = "SBS96"]
setnames(m96, "SBS96", "MutationType")
for (cl in setdiff(names(m96), "MutationType")) set(m96, which(is.na(m96[[cl]])), cl, 0L)
fwrite(m96, file.path(sig_in, paste0(PATIENT, ".SBS96.from_R.txt")), sep = "\t")
message("signature inputs: ", sig_in, "  ->  python 01b_sarek_signatures.py")

message("done: ", od)
