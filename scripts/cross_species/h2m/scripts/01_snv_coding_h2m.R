# ---------------------------------------------------------------------------
# 01_snv_coding_h2m.R - coding SNVs/indels (PURPLE PASS, PAVE canonical coding effect in
# snv.coding_effects) modelled in mouse with H2M (bioh2m 1.0.3, direction h2m).
#
#   Rscript scripts/01_snv_coding_h2m.R [--config ...] [--release ...] [--skip-validation]
#
# H2M builds: human GRCh38 (ver=38) against mouse GRCm39 - H2M's mouse side is GRCm39 only (its
# bundled transcripts are GENCODE vM33). GRCm38 coordinates come from a liftOver of the GRCm39
# result (mm39ToMm10, liftover.min_match_mouse). Orthologs: Ensembl Compara (config release),
# not H2M's own name table, so everything stays keyed on Ensembl gene IDs. Transcripts: PAVE's
# transcript if it has a CDS in the release GTF, else the Ensembl canonical one; mouse: the
# Ensembl canonical transcript of each ortholog. H2M itself is Python - this script writes its
# input table and calls scripts/h2m_model.py (config `python`); everything else is R.
#
# H2M validation first: config snv.h2m_validation (KRAS G12D, TP53 R175H, CTNNB1 S45F, PIK3CA
# H1047R, BRAF V600E) goes through the same call; REF must match human_fasta and H2M's mouse
# HGVSp must equal expected_mouse. Any FAIL stops the script unless --skip-validation.
#
# Input   <oncoanalyser>/purple/*.purple.somatic.vcf.gz, *.purple.driver.catalog.somatic.tsv
# Output  human_predictions/
#   expected_mouse_snv.tsv   one row per human variant x mouse ortholog:
#     variant_id chrom pos ref alt vtype gene_symbol_h gene_id_h transcript_h hgvsc_h hgvsp_h
#     pave_effect pave_coding_effect reported hotspot purple_driver_gene purple_af biallelic
#     homology_type n_mouse gene_id_m symbol_m perc_id_h_to_m perc_id_m_to_h tx_h tx_m
#     h2m_class h2m_statement n_alternatives
#     same_dna_<build>_{chrom,start,end,ref,alt} same_dna_hgvsc_m same_dna_hgvsp_m   (H2M "_ori")
#     same_protein_<build>_{chrom,start,end,ref,alt} same_protein_hgvsc_m same_protein_hgvsp_m
#     alternatives_hgvsc_m alternatives_hgvsp_m (all H2M routes, ';'; class 1 = same protein via other codons;
#     at most snv.h2m_max_alternative) uncertain_effect uncertain_reason prediction_status
#     equivalence_level functional_equivalence mouse_aligned_residue mouse_hypothesis counterpart
# equivalence_level (H2M's class is kept as it is; this is a separate interpretation):
#   exact       H2M class 0 - same nucleotide change                         functional_equivalence = exact
#   functional  H2M class 1 - same protein change via other codons             supported
#               H2M class 2 - modelled, effect may differ                      predicted
#   gene_level  H2M could not model (class 3-6 / failed) a truncating or splice variant of a gene with a 1:1
#               ortholog: the mouse gene can be disrupted (near mouse_aligned_residue), but equivalence to this
#               specific human variant is not shown                            predicted_not_validated
#   none        everything else (missense not modelable, no ortholog, ...)    not_applicable
# counterpart: yes (exact / functional class 1), gene_level, uncertain, no.
#   h2m_raw.tsv, h2m_failed.tsv   H2M output as returned (all alternatives) / inputs it rejected
#   h2m_validation.tsv            the validation variants with H2M's answer and PASS/FAIL
#   summary/snv_coding.tsv
# prediction_status: modelled_original (H2M class 0), modelled_alternative (1),
#   modelled_effect_uncertain (2), not_modelable (3-6), h2m_failed, no_mouse_transcript,
#   not_modelled_homology, no_ortholog, not_in_release.
# uncertain_effect: frameshift, splice, or H2M class 2.
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

H2M_STATUS <- c(`0` = "modelled_original", `1` = "modelled_alternative", `2` = "modelled_effect_uncertain",
                `3` = "not_modelable", `4` = "not_modelable", `5` = "not_modelable", `6` = "not_modelable")
COUNTERPART <- c(modelled_original = "yes", modelled_alternative = "yes", modelled_effect_uncertain = "uncertain",
                 not_modelable = "no", h2m_failed = "uncertain", no_mouse_transcript = "uncertain",
                 not_modelled_homology = "uncertain", no_ortholog = "no", not_in_release = "no")

# gene -> transcript with a CDS: preferred (version-stripped match) if valid, else Ensembl canonical
pick_tx <- function(gene_ids, tx, preferred = NULL) {
  t <- tx[has_cds == TRUE]
  t[, tid := strip_version(transcript_id)]
  canon <- t[canonical == TRUE][!duplicated(gene_id)]
  out <- canon$transcript_id[match(gene_ids, canon$gene_id)]
  if (!is.null(preferred)) {
    pref <- strip_version(fcoalesce(as.character(preferred), ""))
    out <- fifelse(pref %chin% t$tid, pref, out)
  }
  out
}

h2m_input <- function(d) cbind(d[, .(index, tx_id_h = tx_h, tx_id_m = tx_m)], vcf_to_maf(d$pos, d$ref, d$alt))

main <- function() {
  a <- parse_cli("skip-validation")
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("01_snv_coding_h2m", P)
  onco <- cfg$oncoanalyser_dir
  hb <- cfg$h2m_mouse_build
  M <- P$mouse[[hb]]

  vcf <- find_one(file.path(onco, "purple"), "\\.purple\\.somatic\\.vcf\\.gz$", "PURPLE somatic VCF", log)
  drv_f <- find_one(file.path(onco, "purple"), "\\.purple\\.driver\\.catalog\\.somatic\\.tsv$", "PURPLE driver catalog", log, required = FALSE)
  drv_genes <- if (is.null(drv_f)) character() else unique(fread(drv_f)$gene)

  v <- read_purple_somatic(vcf, log)[pave_coding_effect %chin% cfg$snv$coding_effects]
  log("coding (%s): %d variants", paste(cfg$snv$coding_effects, collapse = ","), nrow(v))

  hg <- load_genes(P$human_genes); htx <- load_genes(P$human_tx)
  mg <- load_genes(M$genes); mtx <- load_genes(M$tx)
  orth <- read_orthologs(P$orthologs)

  v[, gene_id_h := gene_ids_for(gene_symbol_h, hg, htx, transcript_h)]
  v[, purple_driver_gene := gene_symbol_h %chin% drv_genes]
  v[, uncertain_reason := fcase(pave_coding_effect == "SPLICE", "splice",
                                pave_coding_effect == "NONSENSE_OR_FRAMESHIFT" & nchar(ref) != nchar(alt) & is_frameshift(ref, alt), "frameshift",
                                default = "")]
  log("no human gene ID: %d", v[is.na(gene_id_h), .N])

  om <- map_orthologs(v$gene_id_h, orth)[, .(gene_id_h, gene_id_m, symbol_m, homology_type, n_mouse, perc_id_h_to_m, perc_id_m_to_h)]
  x <- om[v, on = "gene_id_h", allow.cartesian = TRUE]
  x[is.na(homology_type), homology_type := "not_in_release"]
  x[is.na(n_mouse), n_mouse := 0L]
  x[, tx_h := pick_tx(gene_id_h, htx, transcript_h)]
  x[, tx_m := pick_tx(gene_id_m, mtx)]
  x[, index := seq_len(.N) - 1L]
  run <- x[!is.na(gene_id_m) & homology_type %chin% cfg$snv$model_homology & !is.na(tx_h) & !is.na(tx_m)]
  log("H2M input: %d variant x ortholog pairs (%d variants)", nrow(run), uniqueN(run$variant_id))

  # ---- validation set ------------------------------------------------------------------------
  val <- rbindlist(lapply(cfg$snv$h2m_validation, as.data.table))
  val[, pos := as.integer(pos)]
  val[, gene_id_h := gene_ids_for(gene, hg, htx)]
  val[orth[homology_type == "ortholog_one2one"], on = "gene_id_h", gene_id_m := i.gene_id_m]
  val[, `:=`(tx_h = pick_tx(gene_id_h, htx), tx_m = pick_tx(gene_id_m, mtx), index = -seq_len(.N))]
  fa <- Rsamtools::FaFile(cfg$human_fasta)
  val[, fasta_ref := as.character(Rsamtools::scanFa(fa, GenomicRanges::GRanges(chrom, IRanges::IRanges(pos, pos + nchar(ref) - 1L))))]
  val[, ref_in_fasta := toupper(fasta_ref) == ref]

  # ---- one H2M call (Python helper) ---------------------------------------------------------------
  inp <- rbind(h2m_input(val[!is.na(tx_h) & !is.na(tx_m)]), h2m_input(run))
  f_in <- file.path(P$tmp, "h2m_input.tsv"); f_raw <- file.path(P$hp, "h2m_raw.tsv"); f_fail <- file.path(P$hp, "h2m_failed.tsv")
  write_tsv(inp, f_in, log)
  py <- if (is.null(cfg$python)) "python" else cfg$python
  args <- c(file.path(H2M_SCRIPTS, "h2m_model.py"), "--input", f_in, "--out", f_raw, "--failed", f_fail,
            "--human-fasta", cfg$human_fasta, "--mouse-fasta", M$fasta, "--human-gtf", P$human_gtf, "--mouse-gtf", M$gtf,
            "--workdir", file.path(P$tmp, "h2m"), "--max-alternative", cfg$snv$h2m_max_alternative,
            if (!isTRUE(cfg$snv$h2m_search_alternative)) "--no-alternative")
  log("running %s %s", py, paste(args, collapse = " "))
  if (system2(py, shQuote(args)) != 0) stop("h2m_model.py failed", call. = FALSE)
  raw <- fread(f_raw, na.strings = c("", "NA", "None", "nan"))
  fail <- fread(f_fail, na.strings = c("", "NA", "None", "nan"))
  failed <- if (nrow(fail)) fail$index else integer()
  log("H2M: %d result rows, %d inputs failed", nrow(raw), length(failed))

  first <- raw[order(index)][, .SD[1], by = index]
  # each H2M row of a variant is one modelling route: for class 1 the alternatives give the same protein
  # change through different codons, so their nucleotide changes (HGVSc_m) are what differs
  alts <- raw[, .(n_alternatives = .N, alternatives_hgvsc_m = paste(unique(na.omit(HGVSc_m)), collapse = ";"),
                  alternatives_hgvsp_m = paste(unique(na.omit(HGVSp_m)), collapse = ";")), by = index]

  # ---- validation verdict ---------------------------------------------------------------------
  val[first, on = "index", `:=`(h2m_class = i.class, hgvsp_h_h2m = i.HGVSp_h, hgvsp_m_h2m = i.HGVSp_m)]
  if (!"hgvsp_m_h2m" %chin% names(val)) val[, hgvsp_m_h2m := NA_character_]
  val[, pass := ref_in_fasta & same_protein_change(hgvsp_m_h2m, expected_mouse)]
  write_tsv(val[, !"index"], file.path(P$hp, "h2m_validation.tsv"), log)
  for (i in seq_len(nrow(val)))
    log("H2M validation %-6s %-7s -> %-16s expected %-7s REF %s  %s", val$gene[i], val$protein[i], val$hgvsp_m_h2m[i],
        val$expected_mouse[i], if (val$ref_in_fasta[i]) "ok" else "MISMATCH", if (val$pass[i]) "PASS" else "FAIL")
  if (!all(val$pass)) {
    msg <- sprintf("H2M validation: %d of %d known variants failed - see h2m_validation.tsv", sum(!val$pass), nrow(val))
    if (isTRUE(a$skip_validation)) log("%s (continuing: --skip-validation)", msg, level = "WARNING")
    else { log("%s; rerun with --skip-validation to continue anyway", msg, level = "ERROR"); quit(status = 1) }
  }

  # ---- one row per input: original (same DNA change) + first alternative (same protein change)
  first <- first[index >= 0]
  first[, chr_m := fifelse(is.na(chr_m), NA_character_, norm_chrom(chr_m))]
  x[first, on = "index", `:=`(h2m_class = i.class, h2m_statement = i.statement,
                              same_dna_chrom = i.chr_m, same_dna_start = i.start_m_ori, same_dna_end = i.end_m_ori,
                              same_dna_ref = i.ref_seq_m_ori, same_dna_alt = i.alt_seq_m_ori,
                              same_dna_hgvsc_m = i.HGVSc_m_ori, same_dna_hgvsp_m = i.HGVSp_m_ori,
                              same_protein_chrom = i.chr_m, same_protein_start = i.start_m, same_protein_end = i.end_m,
                              same_protein_ref = i.ref_seq_m, same_protein_alt = i.alt_seq_m,
                              same_protein_hgvsc_m = i.HGVSc_m, same_protein_hgvsp_m = i.HGVSp_m)]
  x[alts, on = "index", `:=`(n_alternatives = i.n_alternatives, alternatives_hgvsc_m = i.alternatives_hgvsc_m,
                                 alternatives_hgvsp_m = i.alternatives_hgvsp_m)]
  for (tag in c("same_dna", "same_protein")) {
    old <- paste0(tag, "_", c("chrom", "start", "end", "ref", "alt"))
    if (all(old %chin% names(x))) setnames(x, old, paste0(tag, "_", hb, "_", c("chrom", "start", "end", "ref", "alt")))
  }

  # ---- other builds: liftOver of the H2M build's coordinates -------------------------------------
  for (build in setdiff(names(P$mouse), hb)) {
    if (!(hb == "GRCm39" && build == "GRCm38")) { log("no chain for %s -> %s, coordinates left empty", hb, build, level = "WARNING"); next }
    for (tag in c("same_dna", "same_protein")) {
      cc <- paste0(tag, "_", hb, "_", c("chrom", "start", "end"))
      if (!all(cc %chin% names(x))) next
      bed <- x[!is.na(get(cc[1])) & !is.na(get(cc[2])) & !is.na(get(cc[3])),
               .(chrom = get(cc[1]), start = pmin(as.integer(get(cc[2])), as.integer(get(cc[3]))) - 1L,
                 end = pmax(as.integer(get(cc[2])), as.integer(get(cc[3]))), name = as.character(index))]
      if (!nrow(bed)) next
      lo <- run_liftover(bed, P$mm39_to_mm10, cfg$liftover$min_match_mouse, P$tmp, paste0("snv_", tag, "_mm10"))
      mp <- lo$mapped[!duplicated(name)]
      x[, key := as.character(index)]
      x[mp, on = c(key = "name"), (paste0(tag, "_", build, "_", c("chrom", "start", "end"))) := .(i.chrom, i.start + 1L, i.end)]
      x[, key := NULL]
    }
  }

  if (!"h2m_class" %chin% names(x)) x[, `:=`(h2m_class = NA_integer_, h2m_statement = NA_character_)]
  x[, prediction_status := fcase(
    homology_type %chin% c("no_ortholog", "not_in_release"), homology_type,
    !homology_type %chin% cfg$snv$model_homology, "not_modelled_homology",
    is.na(tx_m) | is.na(tx_h), "no_mouse_transcript",
    index %in% failed | is.na(h2m_class), "h2m_failed",
    default = unname(H2M_STATUS[as.character(h2m_class)]))]
  x[, counterpart := fcoalesce(unname(COUNTERPART[prediction_status]), "uncertain")]
  x[, uncertain_effect := uncertain_reason != "" | h2m_class %in% 2L]
  x[uncertain_reason == "" & h2m_class %in% 2L, uncertain_reason := "h2m_class2"]
  x[mg, on = c(gene_id_m = "gene_id"), symbol_m := fcoalesce(symbol_m, i.gene_name)]

  # ---- equivalence level, next to (never replacing) H2M's class -----------------------------------------
  # exact (class 0) / functional (class 1 supported, class 2 predicted) / gene_level (H2M could not model a
  # truncating or splice variant of a gene with a 1:1 ortholog: the mouse gene can be disrupted, equivalence to
  # this specific human variant is not shown) / none
  x[, equivalence_level := fcase(
    h2m_class %in% 0L, "exact", h2m_class %in% 1:2, "functional",
    prediction_status %chin% c("not_modelable", "h2m_failed") & pave_coding_effect %chin% c("NONSENSE_OR_FRAMESHIFT", "SPLICE") &
      homology_type == "ortholog_one2one" & !is.na(tx_m), "gene_level",
    default = "none")]
  x[, functional_equivalence := fcase(equivalence_level == "exact", "exact", h2m_class %in% 1L, "supported",
                                      h2m_class %in% 2L, "predicted", equivalence_level == "gene_level", "predicted_not_validated",
                                      default = "not_applicable")]
  x[equivalence_level == "gene_level", counterpart := "gene_level"]
  # gene-level rows: mouse residue aligned to the human HGVSp position (BLOSUM62 protein alignment)
  x[, mouse_aligned_residue := NA_integer_]
  gl <- x[equivalence_level == "gene_level"]
  if (nrow(gl)) {
    hprot <- proteins_from(attach_seq(read_tx_features(P$human_gtf, tx_ids = unique(gl$tx_h)), cfg$human_fasta))
    mprot <- proteins_from(attach_seq(read_tx_features(M$gtf, tx_ids = unique(gl$tx_m)), M$fasta))
    hp <- parse_hgvsp(gl$hgvsp_h)
    ar <- vapply(seq_len(nrow(gl)), function(i) {
      h <- hprot[gl$tx_h[i]]; m <- mprot[gl$tx_m[i]]
      if (is.na(h) || is.na(m) || is.na(hp$pos[i])) return(NA_integer_)
      al <- align_proteins(h, m)
      fcoalesce(map_aligned_pos(al$p, al$s, hp$pos[i])$mouse_pos, NA_integer_)
    }, 0L)
    x[gl[, .(index)][, ar := ar], on = "index", mouse_aligned_residue := i.ar]
    log("gene-level equivalence: %d variant x ortholog rows", nrow(gl))
  }
  eff <- fifelse(x$pave_coding_effect == "SPLICE", "splice", "truncating")
  x[, mouse_hypothesis := fcase(
    equivalence_level %chin% c("exact", "functional"),
    paste0(symbol_m, " ", fcoalesce(fifelse(h2m_class %in% 1L, same_protein_hgvsp_m, same_dna_hgvsp_m), same_dna_hgvsc_m, ""),
           " (H2M class ", h2m_class, ")"),
    equivalence_level == "gene_level",
    paste0("loss of function of ", symbol_m, " (", eff, ")",
           fifelse(is.na(mouse_aligned_residue), "", paste0(" near aligned residue ", mouse_aligned_residue)),
           "; gene-level - the exact variant is not modelable (H2M class ", fcoalesce(as.character(h2m_class), "NA"), ")"),
    default = NA_character_)]

  lead <- c("variant_id", "chrom", "pos", "ref", "alt", "vtype", "gene_symbol_h", "gene_id_h", "transcript_h", "hgvsc_h",
            "hgvsp_h", "pave_effect", "pave_coding_effect", "reported", "hotspot", "purple_driver_gene", "purple_af",
            "biallelic", "homology_type", "n_mouse", "gene_id_m", "symbol_m", "perc_id_h_to_m", "perc_id_m_to_h",
            "tx_h", "tx_m", "h2m_class", "h2m_statement", "n_alternatives")
  mid <- sort(grep("^same_(dna|protein)_", names(x), value = TRUE))
  tail <- c("alternatives_hgvsc_m", "alternatives_hgvsp_m", "uncertain_effect", "uncertain_reason", "prediction_status",
            "equivalence_level", "functional_equivalence", "mouse_aligned_residue", "mouse_hypothesis", "counterpart")
  out <- x[order(chrom_order(chrom), pos, variant_id, gene_id_m), intersect(c(lead, mid, tail), names(x)), with = FALSE]
  f <- write_tsv(out, file.path(P$hp, "expected_mouse_snv.tsv"), log)

  # per human variant: the best counterpart over its orthologs
  best <- out[order(match(counterpart, c("yes", "gene_level", "uncertain", "no")))][!duplicated(variant_id)]
  rows <- rbind(srow("variants_total", nrow(best), "info"),
                # H2M's own result per variant (unchanged by the equivalence level)
                best[, .(n = .N), by = prediction_status][order(prediction_status),
                     .(category = paste0("h2m_", prediction_status), n, counterpart = "info")],
                # counterpart per variant: exact / functional = yes, gene_level = gene-level only, rest uncertain / no
                best[, .(n = .N), by = .(counterpart, equivalence_level)][order(match(counterpart, c("yes", "gene_level", "uncertain", "no")), equivalence_level),
                     .(category = paste0("equivalence_", equivalence_level), n, counterpart)],
                srow("uncertain_effect", best[uncertain_effect == TRUE, .N], "info"),
                srow("h2m_validation_pass", sum(val$pass), "info"), srow("h2m_validation_fail", sum(!val$pass), "info"))
  write_summary(P, "snv_coding", rows, log)
  write_provenance(P, cfg, "01_snv_coding_h2m", inputs = list(purple_somatic_vcf = vcf, purple_driver_catalog = drv_f),
                   resources = list(orthologs = P$orthologs, human_gtf = P$human_gtf, mouse_gtf = M$gtf, mm39_to_mm10 = P$mm39_to_mm10),
                   outputs = c(f, f_raw, f_fail),
                   extra = list(h2m = list(human_ver = 38, mouse_build = hb, transcripts = "PAVE else Ensembl canonical",
                                           orthologs = paste("Ensembl Compara", cfg$ensembl_release),
                                           validation_pass = sum(val$pass), validation_n = nrow(val))))
  log("done")
}

main()
