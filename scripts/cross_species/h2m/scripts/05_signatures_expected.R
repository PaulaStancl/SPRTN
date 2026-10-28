# ---------------------------------------------------------------------------
# 05_signatures_expected.R - human SBS96 spectrum -> expected spectrum in mouse sequence
# opportunity, and the SIGS signature exposures.
#
#   Rscript scripts/05_signatures_expected.R [--config ...] [--release ...]
#
# SBS96 from PURPLE's PASS single-base SNVs (MNV records left out, as in analysis/01b: a doublet
# counts in DBS78 only), context from human_fasta (Rsamtools). Built here rather than with
# SigProfilerMatrixGenerator (its own multi-GB GRCh38 install); if config sbs96_reference (the
# SigProfiler matrix of the same SAGE set from analysis/01b) exists, the two are compared.
# Trinucleotide counts of the autosomes (A/C/G/T triplets only, Biostrings) of GRCh38 and each mouse
# build are cached in <resources_dir>/trinuc/. Expected mouse spectrum: each channel x (mouse / human
# frequency of its pyrimidine-centred trinucleotide), scaled back to the human total. With config
# signature_definitions (COSMIC-format 96 x K), each signature's exposure is rescaled the same way
# (exposure x sum_c profile_c x ratio_c, renormalised).
#
# Input   <oncoanalyser>/purple/*.purple.somatic.vcf.gz, <oncoanalyser>/sigs/*.sig.allocation.tsv
# Output  human_predictions/
#   expected_mouse_sbs96.tsv       build MutationType trinucleotide human_count human_prop human_tri_freq
#                                  mouse_tri_freq ratio_mouse_human expected_mouse_count expected_mouse_prop
#   expected_signatures.tsv        signature human_allocation human_pct [<build>_expected_pct]
#   trinucleotide_frequencies.tsv  trinucleotide, count / frequency per genome
#   sbs96_crosscheck.tsv           per channel: ours vs sbs96_reference (if present)
#   summary/signatures.tsv
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("05_signatures_expected", P)
  vcf <- find_one(file.path(cfg$oncoanalyser_dir, "purple"), "\\.purple\\.somatic\\.vcf\\.gz$", "PURPLE somatic VCF", log)
  f_sig <- find_one(file.path(cfg$oncoanalyser_dir, "sigs"), "\\.sig\\.allocation\\.tsv$", "SIGS allocation", log, required = FALSE)

  v <- read_purple_somatic(vcf, log)[vtype == "SNV"]
  fa <- Rsamtools::FaFile(cfg$human_fasta)
  v[, ctx := toupper(as.character(Rsamtools::scanFa(fa, GenomicRanges::GRanges(chrom, IRanges::IRanges(pos - 1L, pos + 1L)))))]
  v[, ok := substr(ctx, 2, 2) == ref & !grepl("[^ACGT]", ctx)]
  if (v[ok == FALSE, .N]) log("%d SNVs with REF != FASTA base or N in context - left out", v[ok == FALSE, .N], level = "WARNING")
  lab <- v[ok == TRUE, sbs96_label(ctx, ref, alt)]
  sbs <- setNames(as.numeric(table(factor(lab, levels = SBS96))), SBS96)
  log("SBS96: %d SNVs", sum(sbs))

  ref_f <- cfg$sbs96_reference
  cross <- data.table()
  if (!is.null(ref_f) && file.exists(ref_f)) {
    r <- fread(ref_f)
    col <- setdiff(names(r), "MutationType")[1]
    cross <- data.table(MutationType = SBS96, ours = unname(sbs), reference = r[[col]][match(SBS96, r$MutationType)])
    cross[is.na(reference), reference := 0]
    cross[, difference := ours - reference]
    log("SBS96 vs %s: total %d vs %d, channels differing %d", basename(ref_f), sum(cross$ours), sum(cross$reference), cross[difference != 0, .N])
  } else if (!is.null(ref_f)) log("sbs96_reference not found: %s", ref_f, level = "WARNING")
  write_tsv(cross, file.path(P$hp, "sbs96_crosscheck.tsv"), log)

  tri_h <- genome_trinuc(cfg$human_fasta, 22, file.path(P$trinuc_dir, "GRCh38.autosomes.tsv"), log)
  tri <- list(GRCh38 = tri_h)
  spectra <- list()
  for (build in names(P$mouse)) {
    tri[[build]] <- genome_trinuc(P$mouse[[build]]$fasta, 19, file.path(P$trinuc_dir, paste0(build, ".autosomes.tsv")), log)
    spectra[[build]] <- rescale_spectrum(sbs, tri_h, tri[[build]])[, build := build]
  }
  sp_all <- rbindlist(spectra)
  setcolorder(sp_all, "build")
  f1 <- write_tsv(sp_all, file.path(P$hp, "expected_mouse_sbs96.tsv"), log)
  tf <- data.table(trinucleotide = TRI32)
  for (k in names(tri)) tf[, (paste0(k, c("_count", "_freq"))) := .(unname(tri[[k]][TRI32]), unname(tri[[k]][TRI32] / sum(tri[[k]])))]
  f2 <- write_tsv(tf, file.path(P$hp, "trinucleotide_frequencies.tsv"), log)

  sig <- data.table(signature = character(), human_allocation = numeric(), human_pct = numeric())
  if (!is.null(f_sig)) {
    s <- fread(f_sig)
    setnames(s, tolower(names(s)))
    need_cols(s, c("signature", "allocation"), "sig.allocation.tsv")
    sig <- s[, .(signature, human_allocation = allocation, human_pct = 100 * allocation / sum(allocation))]
    defs_f <- cfg$signature_definitions
    if (!is.null(defs_f) && file.exists(defs_f)) {
      d <- fread(defs_f)
      key <- intersect(c("Type", "MutationType"), names(d))[1]
      d <- d[match(SBS96, d[[key]])]
      miss <- setdiff(sig$signature, names(d))
      if (length(miss)) log("signatures without a definition: %s", paste(miss, collapse = ","), level = "WARNING")
      for (build in names(P$mouse)) {
        ratio <- spectra[[build]]$ratio_mouse_human
        w <- vapply(sig$signature, function(k) if (k %chin% names(d)) sum(d[[k]] * ratio) else NA_real_, 0)
        e <- w * sig$human_allocation
        sig[, (paste0(build, "_expected_pct")) := 100 * e / sum(e, na.rm = TRUE)]
      }
    } else if (!is.null(defs_f)) log("signature_definitions not found: %s", defs_f, level = "WARNING")
    setorder(sig, -human_allocation)
  }
  f3 <- write_tsv(sig, file.path(P$hp, "expected_signatures.tsv"), log)

  sp <- spectra[[cfg$h2m_mouse_build]]
  cosv <- if (sum(sbs) > 0) sum(sp$human_prop * sp$expected_mouse_prop) / sqrt(sum(sp$human_prop^2) * sum(sp$expected_mouse_prop^2)) else NA_real_
  log("cosine human vs expected mouse spectrum (%s): %.4f", cfg$h2m_mouse_build, cosv)
  write_summary(P, "signatures", rbind(srow("snvs_in_sbs96", sum(sbs), "yes"), srow("snvs_left_out", v[ok == FALSE, .N], "no"),
                                       srow("signatures_reported", nrow(sig), "info")), log)
  write_provenance(P, cfg, "05_signatures_expected", inputs = list(purple_somatic_vcf = vcf, sigs_allocation = f_sig),
                   resources = list(human_fasta_fai = paste0(cfg$human_fasta, ".fai")), outputs = c(f1, f2, f3),
                   extra = list(cosine_human_vs_mouse_expected = cosv, sbs96_reference = ref_f))
  log("done")
}

main()
