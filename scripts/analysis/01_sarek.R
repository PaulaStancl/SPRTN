# ---------------------------------------------------------------------------
# 01_sarek.R - sarek output for RJALS, tumour vs normal.
#   SNV / indel: Mutect2 + Strelka2 (+ MuSE, extra run RJALS_vc)     SV: Manta     copy number + purity: ASCAT
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
  strelka_indel = find_one(nd, "somatic_indels\\.norm\\.vcf\\.gz$",         required = FALSE),
  muse          = find_one(nd, "\\.muse\\.norm\\.vcf\\.gz$",                  required = FALSE))
if (is.null(norm)) message("no normalised VCFs (", nd, ") - using sarek's raw ones; run scripts/wgs/06_normalize_vcfs.sh")
files <- c(
  mutect2       = find_one(file.path(vc, "mutect2", PAIR), "\\.mutect2\\.filtered\\.vcf\\.gz$"),
  strelka_snv   = find_one(file.path(vc, "strelka", PAIR), "somatic_snvs\\.vcf\\.gz$",   required = FALSE),
  strelka_indel = find_one(file.path(vc, "strelka", PAIR), "somatic_indels\\.vcf\\.gz$", required = FALSE),
  # MuSE: run later on the same recalibrated CRAMs, in its own sarek run (results/wgs/sarek/RJALS_vc)
  muse          = find_one(file.path(RESULTS, "sarek", paste0(PATIENT, "_vc"), "variant_calling", "muse", PAIR),
                           "\\.muse\\.vcf\\.gz$", required = FALSE),
  manta         = find_one(file.path(vc, "manta",   PAIR), "\\.manta\\.somatic_sv\\.vcf\\.gz$", required = FALSE),
  ascat_seg     = find_one(file.path(vc, "ascat",   PAIR), "\\.segments\\.txt$"),
  ascat_pp      = find_one(file.path(vc, "ascat",   PAIR), "\\.purityploidy\\.txt$")
)
files <- files[!is.na(files)]
if (!is.null(norm)) files[names(norm)[!is.na(norm)]] <- norm[!is.na(norm)]
writeLines(paste(names(files), files, sep = "\t"), file.path(od, "inputs_used.tsv"))

# ---- 2. SNV / indel: Mutect2, Strelka2 and MuSE (SNVs only) ----------------------
calls <- rbindlist(list(
  load_calls(files, "mutect2", "mutect2"),
  load_calls(files, c("strelka_snv", "strelka_indel"), "strelka"),     # NULL (skipped) if not found
  load_calls(files, "muse", "muse")                                     # MuSE: SNVs only, FILTER PASS / Tier2-5
), fill = TRUE)                                                         # the callers' metric columns differ
# qc_vcf/raw:  every record, PASS and filtered - FILTER outcome, reasons, depth / VAF / score
# qc_vcf/pass: PASS calls only - caller overlap, VAF, spectrum, qcVCF plots (section 5)
raw_dir  <- file.path(od, "qc_vcf", "raw");  dir.create(raw_dir,  recursive = TRUE, showWarnings = FALSE)
pass_dir <- file.path(od, "qc_vcf", "pass"); dir.create(pass_dir, recursive = TRUE, showWarnings = FALSE)
raw_qc(calls, raw_dir)
pass <- snv_indel_summary(calls, pass_dir)
calls_pass <- pass$pass      # PASS records as the callers wrote them (Mutect2 MNVs whole)
atom       <- pass$atom      # the same with MNVs split into SNVs - for comparing callers

# trinucleotide context of every PASS SNV (NC_3 + strand, SBS96 class) from the FASTA
add_sbs96(atom)
# which callers have each PASS mutation (key CHROM:POS:REF:ALT, MNVs split), and the record
# each row of atom comes from (REC_KEY - an MNV's split bases share their MNV's key)
atom[, CALLERS := paste(sort(unique(CALLER)), collapse = "+"), by = .(CHROM, POS, REF, ALT1)]
atom[, REC_KEY := fifelse(FROM_MNV, MNV_KEY, paste0(CHROM, ":", POS, ":", REF, ">", ALT1))]
in_strelka <- function(callers) grepl("(^|\\+)strelka(\\+|$)", callers)   # CALLERS contains strelka (any other callers too)

# ---- 3. structural variants: Manta -> sv/ ------------------------------------------
sv_dir <- file.path(od, "sv"); dir.create(sv_dir, showWarnings = FALSE)
if ("manta" %in% names(files)) {
  sv <- read_sv_vcf(files[["manta"]])[, CALLER := "manta"]
  sv_summary(sv, sv_dir)
}

# ---- 4. copy number and purity: ASCAT -> cnv/ ------------------------------------------
cn_dir <- file.path(od, "cnv"); dir.create(cn_dir, showWarnings = FALSE)
seg <- ascat_segments(files[["ascat_seg"]])
cn_summary(seg, cn_dir)
purity <- fread(files[["ascat_pp"]])
print(purity); fwrite(purity, file.path(cn_dir, "ascat_purity_ploidy.csv"))
# outputs of the earlier flat layout, now in sv/ and cnv/
unlink(file.path(od, c("sv_counts.csv", "sv_types.pdf", "cn_segments.csv", "cn_genome.pdf",
                       "cn_state_fraction.csv", "cn_loh_fraction.csv", "ascat_purity_ploidy.csv")))

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

  # 96-context profile, one column per caller, rows all / by number of callers
  # rows by how many SNV callers found each SNV: unique / 2 of 3 / all 3 (the last is the same
  # mutations in every column; "shared" alone would mean "with any other caller", different per caller)
  n_snv_callers <- uniqueN(atom[MUTTYPE == "SNV", CALLER])
  atom[, SHARING := sharing_label(lengths(strsplit(CALLERS, "+", fixed = TRUE)), n_snv_callers)]
  plot_96context(atom, pass_dir, rowsplit = "SHARING", roworder = sharing_levels(n_snv_callers))
  atom[, SHARING := NULL]
}

# ---- 5b. Mutect2 read-orientation artefacts: profile by ROQ ---------------------------------------
# ROQ (Mutect2 INFO): Phred-scaled quality that the ALT allele is NOT a read-orientation artefact
# (LearnReadOrientationModel: e.g. 8-oxoG oxidation -> C>A, deamination -> C>T, seen on one read
# orientation only). FilterMutectCalls has no fixed ROQ cut-off, so low-ROQ calls can PASS.
# If they are artefacts, the low-ROQ group has excess C>A (or C>T) and is rarely confirmed by
# another caller. SPRTN_ROQ_CUT (default 20 = 1% artefact probability) sets the split.
# Output (qc_vcf/pass/): snv_96context_mutect2_by_ROQ.pdf/.csv, mutect2_roq_spectrum.pdf,
#   mutect2_roq_spectrum.csv (class share per group, Fisher test per class, BH), mutect2_roq_confirmed.csv
if ("info_ROQ" %in% names(atom) && atom[CALLER == "mutect2" & MUTTYPE == "SNV" & !is.na(info_ROQ), .N] > 0) {
  ROQ_CUT <- as.numeric(Sys.getenv("SPRTN_ROQ_CUT", "20"))
  lv  <- c(sprintf("ROQ < %g", ROQ_CUT), sprintf("ROQ >= %g", ROQ_CUT))
  m2  <- atom[CALLER == "mutect2" & MUTTYPE == "SNV" & !is.na(SBS96) & !is.na(info_ROQ)]
  m2[, ROQ_GROUP := factor(fifelse(as.numeric(info_ROQ) < ROQ_CUT, lv[1], lv[2]), levels = lv)]
  print(m2[, .N, by = ROQ_GROUP])
  plot_96context(m2[, .(CHROM, POS, REF, ALT1, SBS96, ROQ_GROUP = as.character(ROQ_GROUP),
                        CALLER = "Mutect2 PASS SNVs by read-orientation quality (ROQ)")],
                 pass_dir, rowsplit = "ROQ_GROUP", name = "snv_96context_mutect2_by_ROQ", show_all = FALSE, roworder = lv)

  # 6-class spectrum per group; Fisher test per class: share of this class, low vs high ROQ
  m2[, class := substr(SBS96, 3, 5)]
  sp <- dcast(m2[, .N, by = .(class, ROQ_GROUP)], class ~ ROQ_GROUP, value.var = "N", fill = 0)
  setnames(sp, lv, c("n_low", "n_high"))
  tot <- sp[, .(L = sum(n_low), H = sum(n_high))]
  sp[, `:=`(pct_low = round(100 * n_low / tot$L, 1), pct_high = round(100 * n_high / tot$H, 1))]
  sp[, p := mapply(function(a, b) fisher.test(matrix(c(a, tot$L - a, b, tot$H - b), 2))$p.value, n_low, n_high)]
  sp[, p_BH := signif(p.adjust(p, "BH"), 3)][, p := signif(p, 3)]
  print(sp); fwrite(sp, file.path(pass_dir, "mutect2_roq_spectrum.csv"))
  sl <- melt(sp, id.vars = c("class", "p_BH"), measure.vars = c("pct_low", "pct_high"), variable.name = "grp", value.name = "pct")
  sl[, grp := factor(fifelse(grp == "pct_low", sprintf("%s (n=%d)", lv[1], tot$L), sprintf("%s (n=%d)", lv[2], tot$H)))]
  lab <- sp[, .(class, y = pmax(pct_low, pct_high) + 2, txt = fifelse(p_BH < 0.001, "p<0.001", sprintf("p=%.2g", p_BH)))]
  save_plot(ggplot(sl, aes(class, pct, fill = grp)) + geom_col(position = position_dodge(preserve = "single"), width = 0.75) +
              geom_text(data = lab, aes(class, y, label = txt), inherit.aes = FALSE, size = 3) +
              scale_fill_manual(values = c("#d62728", "grey60")) +
              labs(x = NULL, y = "% of the group's SNVs", fill = NULL,
                   title = "Mutect2 PASS SNVs: substitution spectrum by ROQ",
                   subtitle = "artefacts: excess C>A (8-oxoG) or C>T (deamination) at low ROQ\nFisher test per class (share low vs high ROQ), BH-adjusted") +
              theme(legend.position = "bottom"),
            "mutect2_roq_spectrum", pass_dir, w = 7, h = 4.5)

  # confirmed by another caller? (artefacts are mostly Mutect2-only)
  m2[, confirmed := lengths(strsplit(CALLERS, "+", fixed = TRUE)) > 1]
  cf <- m2[, .(snvs = .N, confirmed = sum(confirmed), confirmed_pct = round(100 * mean(confirmed), 1),
               C_to_A_pct = round(100 * mean(class == "C>A"), 1)), by = ROQ_GROUP][order(ROQ_GROUP)]
  print(cf); fwrite(cf, file.path(pass_dir, "mutect2_roq_confirmed.csv"))
} else message("no Mutect2 ROQ values - section 5b skipped")
# ---- 6. clonal vs subclonal: PyClone-VI clusters (tumourevo) on the Mutect2 calls ----
# tumourevo ran PyClone-VI on sarek's Mutect2 PASS calls (autosomes; copy number and purity
# from ASCAT). Its clusters are joined back to the Mutect2 records here (mutation_id is
# RJALS:<chrom>:<VCF POS>:<ALT>), and the clusters are compared on what separates a real
# subclone from a neutral tail or low-VAF false positives: support by the second caller
# (Strelka), allele fraction, alt reads, depth, Mutect2's TLOD, PyClone's assignment
# certainty, karyotype. PyClone-VI's result is the
# same with or without tumourevo's CNAqc filter (same 4,928 mutations), so the unfiltered
# run is read; set SPRTN_TEVO_RUN=RJALS_cnaqcPASS for the other.
cd_dir <- file.path(od, "clonality"); dir.create(cd_dir, showWarnings = FALSE)
py_dir <- file.path(RESULTS, "tumourevo", Sys.getenv("SPRTN_TEVO_RUN", PATIENT), "subclonal_deconvolution", "pyclonevi")
py_fit <- if (dir.exists(py_dir)) find_one(py_dir, "_best_fit\\.txt$", required = FALSE) else NA_character_
if (!is.na(py_fit)) {
  fit <- fread(py_fit)
  py_in <- find_one(py_dir, "_pyclone_input\\.tsv$", required = FALSE)   # not *_all_samples.tsv
  if (!is.na(py_in))
    fit <- unique(fread(py_in)[, .(mutation_id, major_cn, minor_cn)])[fit, on = "mutation_id"]
  # clusters named by their mean cellular prevalence (CP): >= 0.9 clonal, below that subclonal
  cl <- fit[, .(n = .N, CP = mean(cellular_prevalence)), by = cluster_id][order(-CP)]
  cl[, CLONE := fifelse(CP >= 0.9, "clonal", sprintf("subclonal_CP%02.0f", 100 * CP))]
  cl[, CLONE := make.unique(CLONE, sep = "_")]
  print(cl); fwrite(cl, file.path(cd_dir, "pyclone_clusters.csv"))
  fit <- cl[, .(cluster_id, CLONE)][fit, on = "cluster_id"]
  ids <- tstrsplit(fit$mutation_id, ":", fixed = TRUE)
  fit[, `:=`(CHROM = paste0("chr", ids[[2]]), POS = as.integer(ids[[3]]), ALT1 = ids[[4]])]

  # Mutect2 PASS records + Strelka support of each (an MNV: all / some / none of its bases)
  m2 <- calls_pass[CALLER == "mutect2"]
  m2[, REC_KEY := paste0(CHROM, ":", POS, ":", REF, ">", ALT1)]
  sup <- atom[CALLER == "mutect2", .(STRELKA = if (all(in_strelka(CALLERS))) "PASS in Strelka"
                                               else if (any(in_strelka(CALLERS))) "partly (MNV)"
                                               else "not PASS in Strelka"), by = REC_KEY]
  m2 <- sup[m2, on = "REC_KEY"]
  keep <- intersect(c("CHROM", "POS", "REF", "ALT1", "MUTTYPE", "REC_KEY", "STRELKA", "VAF", "N_VAF", "T_ALT",
                      "t_DP", "n_DP", "info_TLOD", "info_ECNT", "info_MBQ", "info_MMQ", "info_MPOS", "info_POPAF"), names(m2))
  clon <- m2[, ..keep][fit, on = .(CHROM, POS, ALT1), nomatch = NULL]
  message(sprintf("PyClone-VI: %d mutations, %d matched to Mutect2 PASS records", nrow(fit), nrow(clon)))
  if (nrow(clon) < nrow(fit))   # e.g. indels shifted by bcftools norm (tumourevo read the raw VCF)
    print(fit[!clon, on = "mutation_id", .N, by = .(type = fifelse(nchar(ALT1) == 1, "SNV-like ALT", "longer ALT"))])
  num1 <- function(x) suppressWarnings(as.numeric(sub(",.*", "", x)))
  for (m in intersect(c("info_TLOD", "info_ECNT", "info_MBQ", "info_MMQ", "info_MPOS", "info_POPAF", "t_DP", "n_DP"), names(clon)))
    set(clon, j = m, value = num1(clon[[m]]))     # Number=A / R fields: first (tumour / ALT) value
  clon[, KARYOTYPE := if ("major_cn" %in% names(clon)) paste0(major_cn, ":", minor_cn) else NA_character_]
  clon[, CLONE := factor(CLONE, levels = cl$CLONE)]
  fwrite(clon, file.path(cd_dir, "clonality_mutations.csv"))

  # per cluster: size, caller agreement, read evidence
  summ <- clon[, .(n = .N, pct_PASS_in_Strelka = round(100 * mean(STRELKA == "PASS in Strelka"), 1),
                   median_VAF = median(VAF, na.rm = TRUE), median_alt_reads = median(T_ALT, na.rm = TRUE),
                   pct_alt_reads_le5 = round(100 * mean(T_ALT <= 5, na.rm = TRUE), 1),
                   median_tumour_DP = median(t_DP, na.rm = TRUE),
                   median_TLOD = if ("info_TLOD" %in% names(clon)) median(info_TLOD, na.rm = TRUE) else NA_real_,
                   median_assignment_prob = median(cluster_assignment_prob, na.rm = TRUE),
                   pct_SNV = round(100 * mean(MUTTYPE == "SNV"), 1)), by = CLONE][order(CLONE)]
  print(summ); fwrite(summ, file.path(cd_dir, "clonality_summary.csv"))
  kt <- clon[, .N, by = .(CLONE, KARYOTYPE)][, pct := round(100 * N / sum(N), 1), by = CLONE][order(CLONE, -N)]
  fwrite(kt, file.path(cd_dir, "clonality_karyotypes.csv"))

  sv <- clon[, .N, by = .(CLONE, STRELKA)]
  save_plot(ggplot(sv, aes(CLONE, N, fill = STRELKA)) + geom_col(position = "fill") +
              geom_text(aes(label = N), position = position_fill(vjust = 0.5), size = 3) +
              scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
              labs(x = NULL, y = "share of the cluster's Mutect2 PASS calls", fill = NULL,
                   title = "Caller agreement per PyClone-VI cluster"), "clonality_strelka_support", cd_dir, w = 7, h = 5)
  save_plot(ggplot(clon, aes(VAF, fill = CLONE)) + geom_histogram(bins = 60, position = "identity", alpha = 0.6) +
              facet_wrap(~KARYOTYPE, scales = "free_y") +
              labs(x = "tumour VAF (raw)", y = "mutations", title = "VAF per cluster, by karyotype (major:minor)"),
            "clonality_vaf", cd_dir, w = 10, h = 6)
  save_plot(ggplot(clon, aes(cellular_prevalence, fill = CLONE)) + geom_histogram(bins = 50) +
              labs(x = "cellular prevalence (PyClone-VI)", y = "mutations", title = "PyClone-VI cellular prevalence"),
            "clonality_cp", cd_dir, w = 7, h = 4)
  ev <- melt(clon, id.vars = c("CLONE", "STRELKA"), na.rm = TRUE, variable.name = "metric", variable.factor = FALSE,
             measure.vars = intersect(c("T_ALT", "t_DP", "info_TLOD", "cluster_assignment_prob", "info_MPOS"), names(clon)))
  save_plot(ggplot(ev, aes(CLONE, value, fill = STRELKA)) + geom_boxplot(outlier.size = 0.3) +
              facet_wrap(~metric, scales = "free_y") +
              labs(x = NULL, y = NULL, fill = NULL, title = "Read evidence per cluster and caller agreement"),
            "clonality_evidence", cd_dir, w = 11, h = 7)
} else message("no PyClone-VI best fit under ", py_dir, " - section 6 skipped")

# ---- 7. mutational signatures: PASS sets for SigProfiler ------------------------
# Each PASS set as a minimal VCF in signatures/input/vcf/ (one "sample" per file), for
# 01b_sarek_signatures.py: SigProfilerMatrixGenerator builds SBS96 / DBS78 / ID83 from them
# (the only matrices the fits use)
# and SigProfilerAssignment fits COSMIC signatures (cosmic_fit, which replaced
# SigProfilerSingleSample). MNVs are written split: the matrix generator rejoins adjacent
# SNVs into doublets (DBS78) itself, the same way for both callers; SBS96 is built from the
# same sets without the doublet halves (vcf_sbs/, see write_signature_sets()), so a doublet
# counts once. Sets (keys CHROM:POS:REF:ALT):
#   mutect2 / strelka                 all PASS calls of that caller
#   mutect2_strelka                   PASS in both callers (intersection)
sig_in <- file.path(od, "signatures", "input")
write_signature_sets(rbindlist(list(
  mutect2         = atom[CALLER == "mutect2"],
  strelka         = atom[CALLER == "strelka"],
  mutect2_strelka = atom[CALLER == "mutect2" & in_strelka(CALLERS)]   # PASS in both (MuSE may have it too)
), idcol = "SET"), sig_in, "01_sarek.R")
message("signature inputs: ", sig_in, "  ->  python 01b_sarek_signatures.py; Rscript 01c_sarek_signatures_organ.R")

message("done: ", od)
