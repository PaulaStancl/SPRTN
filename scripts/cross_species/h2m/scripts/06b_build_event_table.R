# ---------------------------------------------------------------------------
# 06b_build_event_table.R - event-level master table: one row per human alteration, or one row
# per mouse fragment when its mapping splits (several orthologs, several mouse intervals).
# Built from the same per-layer outputs as 06 (and 05b for pathways).
#
#   Rscript scripts/06b_build_event_table.R [--config ...] [--release ...]
#
# Rows per layer
#   snv_coding     per variant x mouse ortholog (01). Method H2M. Mouse locus = same protein change
#                  if H2M had to use an alternative (class 1), else same DNA change.
#                  confidence: high (class 0), partial (class 1 / 2), ambiguous (1:many), none.
#   snv_noncoding  TERT promoter / promoter / UTR / splice-region variants (02). Method liftOver.
#                  confidence: high (single hit, REF conserved, in the ortholog's promoter or gene),
#                  partial (single hit otherwise), ambiguous (multi-mapped), none (unmapped).
#   cn / loh       per non-neutral (or LOH) human segment x mouse interval of primary_mouse_build (03).
#                  Method ortholog + synteny block. mapped_fraction = fraction of the segment's genes
#                  placed in mouse; confidence from cn_segment_mapping.tsv.
#   sv             per disrupted gene (04). Method ortholog. confidence high (1:1 placed), ambiguous
#                  (1:many), none.
#   fusion         per LINX fusion: a feasibility flag only - hypothesis says no mouse fusion is
#                  predicted; confidence none.
#   signature      one row: SBS96 spectrum rescaled to mouse opportunity (05).
# Every row: event_id layer event human_locus human_genes human_gene_ids mouse_orthologs
# mouse_gene_ids homology mouse_locus_<build> mapping_method mapped_fraction mapping_confidence
# hypothesis expected_direction (of the event: type + CGC role) cgc_genes reactome_top_level mouse_1 mouse_2 mouse_3 (= "pending").
#
# Output  human_predictions/event_table.tsv, summary/event_table.tsv
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

colv <- function(dt, name) if (name %chin% names(dt)) dt[[name]] else rep(NA, nrow(dt))
with_loci <- function(dt, base, builds) cbind(base, dt[, paste0("mouse_locus_", builds), with = FALSE])

locus <- function(chrom, start, end = start) fifelse(is.na(chrom) | is.na(start), NA_character_,
                                                     fifelse(start == end, paste0(chrom, ":", start), paste0(chrom, ":", start, "-", end)))

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("06b_build_event_table", P)
  hb <- cfg$primary_mouse_build      # main build for counterpart calls, event table, report
  builds <- names(P$mouse)
  other <- setdiff(builds, hb)
  rows <- list()

  # ---- coding SNVs -------------------------------------------------------------------------------
  snv <- read_layer(P, "expected_mouse_snv.tsv")
  if (!is.null(snv) && nrow(snv)) {
    s <- copy(snv)
    s[, use := fifelse(h2m_class %in% 1L, "same_protein", "same_dna")]
    for (b in builds) {
      cc <- paste0(c("same_dna_", "same_protein_"), b)
      pick <- function(k) { p <- colv(s, paste0(cc[2], "_", k)); d <- colv(s, paste0(cc[1], "_", k))
                            if (k == "chrom") fifelse(s$use == "same_protein", as.character(p), as.character(d))
                            else fifelse(s$use == "same_protein", as.numeric(p), as.numeric(d)) }
      set(s, j = paste0("mouse_locus_", b), value = locus(pick("chrom"), pick("start"), pick("end")))
    }
    s[, mapping_confidence := fcase(prediction_status == "modelled_original" & homology_type == "ortholog_one2one", "high",
                                    prediction_status %chin% c("modelled_alternative", "modelled_effect_uncertain") & homology_type == "ortholog_one2one", "partial",
                                    prediction_status %like% "^modelled" , "ambiguous", default = "none")]
    # modelled in GRCm39 but the position did not lift over to the main build -> partial
    s[mapping_confidence == "high" & is.na(get(paste0("mouse_locus_", hb))), mapping_confidence := "partial"]
    mprot <- fcoalesce(fifelse(s$use == "same_protein", s$same_protein_hgvsp_m, s$same_dna_hgvsp_m), s$same_dna_hgvsc_m)
    s[, hypothesis := fcase(
      prediction_status %chin% c("no_ortholog", "not_in_release"), paste0("no mouse counterpart: ", gene_symbol_h, " has no mouse ortholog"),
      prediction_status == "not_modelable", paste0("not modelable in mouse (H2M class ", h2m_class, "): ", gene_symbol_h, " ", fcoalesce(hgvsp_h, hgvsc_h)),
      prediction_status %like% "^modelled", paste0(fifelse(uncertain_effect %in% TRUE, "possible ", ""), symbol_m, " ", fcoalesce(mprot, ""),
                                                 " (human ", gene_symbol_h, " ", fcoalesce(hgvsp_h, hgvsc_h, ""), ")"),
      default = paste0(gene_symbol_h, " ", fcoalesce(hgvsp_h, hgvsc_h, ""), ": ", prediction_status))]
    rows$snv <- with_loci(s, s[, .(layer = "snv_coding", event = pave_coding_effect, event_key = variant_id,
                           human_locus = paste0(chrom, ":", pos, ":", ref, ">", alt), human_genes = gene_symbol_h,
                           human_gene_ids = gene_id_h, mouse_orthologs = symbol_m, mouse_gene_ids = gene_id_m,
                           homology = homology_type, mapping_method = "H2M",
                           mapped_fraction = fifelse(prediction_status %like% "^modelled", 1, 0),
                           mapping_confidence = mapping_confidence, hypothesis = hypothesis)], builds)
  }

  # ---- non-coding (gene-associated classes) ---------------------------------------------------------
  nc <- read_layer(P, "expected_mouse_noncoding.tsv")
  if (!is.null(nc) && nrow(nc)) {
    n <- nc[class %chin% NONCODING_GENE_CLASSES]
    stc <- paste0(hb, "_status")
    n[, mapping_confidence := fcase(
      get(stc) == "mapped" & get(paste0(hb, "_ref_conserved")) %in% TRUE &
        (get(paste0(hb, "_in_ortholog_promoter")) %in% TRUE | get(paste0(hb, "_in_ortholog_gene")) %in% TRUE), "high",
      get(stc) == "mapped", "partial", get(stc) == "multi_mapped", "ambiguous", default = "none")]
    for (b in builds) set(n, j = paste0("mouse_locus_", b),
                          value = fifelse(n[[paste0(b, "_status")]] == "multi_mapped", paste0("multi_mapped(", n[[paste0(b, "_n_hits")]], ")"),
                                          locus(n[[paste0(b, "_chrom")]], n[[paste0(b, "_pos")]])))
    n[, hypothesis := fifelse(mapping_confidence == "none", paste0("no mouse counterpart: ", class, " variant near ", fcoalesce(gene_symbol_h, "?"), " does not lift over"),
                              paste0(class, " variant of ", fcoalesce(gene_symbol_h, "?"), " -> ", get(paste0("mouse_locus_", hb)),
                                     fifelse(get(paste0(hb, "_ref_conserved")) %in% TRUE, " (REF conserved)", " (REF differs)")))]
    rows$nc <- with_loci(n, n[, .(layer = "snv_noncoding", event = class, event_key = variant_id,
                          human_locus = paste0(chrom, ":", pos, ":", ref, ">", alt), human_genes = gene_symbol_h, human_gene_ids = gene_id_h,
                          mouse_orthologs = NA_character_, mouse_gene_ids = gene_id_m, homology = homology_type,
                          mapping_method = "liftOver", mapped_fraction = fifelse(get(stc) == "unmapped", 0, 1),
                          mapping_confidence = mapping_confidence, hypothesis = hypothesis)], builds)
  }

  # ---- copy number and LOH: human segment x mouse interval ---------------------------------------------
  segm <- read_layer(P, "cn_segment_mapping.tsv")
  if (!is.null(segm) && nrow(segm)) {
    bed_of <- function(b) { f <- file.path(P$hp, sprintf("expected_mouse_cn.%s.bed", b)); if (file.exists(f)) fread(f) else NULL }
    beds <- lapply(setNames(builds, builds), bed_of)
    explode <- function(bed) bed[, .(segment = unlist(strsplit(source_human_segments, ",", fixed = TRUE))),
                                 by = .(chrom = `#chrom`, start, end, name, direction, genes, human_genes, confidence)]
    hs <- segm[build == hb & (direction %chin% c("gain", "loss") | truthy(loh))]
    ex <- explode(beds[[hb]])
    for (d in c("gain", "loss", "LOH")) {
      sg <- if (d == "LOH") hs[truthy(loh)] else hs[direction == d]
      if (!nrow(sg)) next
      fr <- ex[direction == d][sg[, .(segment, state, length_bp, n_genes, gene_frac_mapped, seg_conf = confidence)], on = "segment"]
      oth <- lapply(other, function(b) { e <- explode(beds[[b]])[direction == d]; e[, .(l = paste(sprintf("%s:%d-%d", chrom, start + 1L, end), collapse = ";")), by = segment] })
      fr[, mouse_genes_n := lengths(strsplit(fcoalesce(genes, ""), ",", fixed = TRUE)) * !is.na(genes)]
      fr[, hypothesis := fifelse(is.na(name), paste0("no mouse counterpart: ", d, " segment ", segment, " has no gene placed in mouse"),
                                 paste0(fcase(d == "gain", "gain", d == "loss", "loss", default = "LOH"), " of ",
                                        vapply(strsplit(genes, ",", fixed = TRUE), function(x) paste(head(x, 6), collapse = ", "), ""),
                                        fifelse(mouse_genes_n > 6, paste0(" (+", mouse_genes_n - 6, ")"), ""), " on mouse ", chrom,
                                        fifelse(d == "LOH", " - not assessable in inbred mice unless heterozygous regions exist", "")))]
      r <- fr[, .(layer = if (d == "LOH") "loh" else "cn", event = if (d == "LOH") "LOH" else state, event_key = paste(d, segment, sep = ":"),
                  human_locus = segment, human_genes, human_gene_ids = NA_character_, mouse_orthologs = genes, mouse_gene_ids = NA_character_,
                  homology = "ortholog_one2one", mapping_method = "ortholog + synteny block",
                  mapped_fraction = gene_frac_mapped, mapping_confidence = fifelse(is.na(name), "none", fcoalesce(seg_conf, confidence)),
                  hypothesis, mouse_interval = name, segment)]
      r[, (paste0("mouse_locus_", hb)) := fr[, fifelse(is.na(name), NA_character_, sprintf("%s:%d-%d", chrom, start + 1L, end))]]
      for (k in seq_along(other)) r[oth[[k]], on = "segment", (paste0("mouse_locus_", other[k])) := i.l]
      r[, event_key := paste(event_key, fcoalesce(mouse_interval, "none"), sep = "|")][, c("mouse_interval", "segment") := NULL]
      rows[[paste0("cn_", d)]] <- r
    }
  }

  # ---- SV-disrupted genes and fusions ------------------------------------------------------------
  sv <- read_layer(P, "expected_mouse_sv_genes.tsv")
  if (!is.null(sv) && nrow(sv)) {
    d <- sv[event %chin% c("HOM_DISRUPTION", "DISRUPTION")]
    d[, mapping_confidence := fcase(prediction == "disruption_of_ortholog", "high", prediction == "one2many", "ambiguous", default = "none")]
    d[, hypothesis := fcase(prediction == "disruption_of_ortholog", paste0(fifelse(event == "HOM_DISRUPTION", "homozygous ", ""), "disruption of ", get(paste0(hb, "_symbol"))),
                            prediction == "one2many", paste0("disruption of one of the mouse orthologs of ", gene),
                            default = paste0("no mouse counterpart: ", gene, " (", prediction, ")"))]
    for (b in builds) set(d, j = paste0("mouse_locus_", b), value = locus(d[[paste0(b, "_chrom")]], d[[paste0(b, "_start")]], d[[paste0(b, "_end")]]))
    rows$sv <- with_loci(d, d[, .(layer = "sv", event = event, event_key = paste0("SV:", sv_ids, ":", gene), human_locus = paste0("SV ", sv_ids),
                          human_genes = gene, human_gene_ids = gene_id_h, mouse_orthologs = get(paste0(hb, "_symbol")), mouse_gene_ids = fcoalesce(gene_id_m, mouse_targets),
                          homology = homology, mapping_method = "ortholog", mapped_fraction = fifelse(mapping_confidence == "none", 0, 1),
                          mapping_confidence = mapping_confidence, hypothesis = hypothesis)], builds)
  }
  fus <- read_layer(P, "expected_mouse_fusions.tsv")
  if (!is.null(fus) && nrow(fus)) {
    fcol <- paste0(hb, "_fusion_feasibility")
    rows$fusion <- fus[, .(layer = "fusion", event = "FUSION", event_key = paste0("FUSION:", name), human_locus = name,
                           human_genes = paste(gene5, gene3, sep = ","), human_gene_ids = paste(gene_id5_h, gene_id3_h, sep = ","),
                           mouse_orthologs = NA_character_, mouse_gene_ids = paste(gene_id5_m, gene_id3_m, sep = ","),
                           homology = NA_character_, mapping_method = "ortholog (feasibility flag)", mapped_fraction = NA_real_,
                           mapping_confidence = "none",
                           hypothesis = paste0("no mouse fusion predicted (feasibility flag: ", get(fcol), "; human ", human_sv_class, ")"))]
  }

  # ---- signature spectrum ------------------------------------------------------------------------
  sb <- read_layer(P, "expected_mouse_sbs96.tsv")
  if (!is.null(sb) && nrow(sb)) {
    x <- sb[build == hb]
    cosv <- x[, sum(human_prop * expected_mouse_prop) / sqrt(sum(human_prop^2) * sum(expected_mouse_prop^2))]
    sig <- read_layer(P, "expected_signatures.tsv")
    top <- if (!is.null(sig) && nrow(sig)) paste(head(sprintf("%s %.0f%%", sig$signature, sig$human_pct), 3), collapse = ", ") else "n/a"
    rows$sig <- data.table(layer = "signature", event = "SBS96", event_key = "SBS96", human_locus = "genome",
                           human_genes = NA_character_, human_gene_ids = NA_character_, mouse_orthologs = NA_character_,
                           mouse_gene_ids = NA_character_, homology = NA_character_, mapping_method = "trinucleotide rescaling",
                           mapped_fraction = 1, mapping_confidence = "high",
                           hypothesis = sprintf("mouse SBS96 like the human one (%d SNVs; cosine to rescaled %.3f); top signatures %s",
                                                as.integer(sum(x$human_count)), cosv, top))
  }

  et <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(et)) stop("no events - run 01-05 first", call. = FALSE)

  # ---- pathway columns (05b) --------------------------------------------------------------------
  gp <- read_layer(P, "gene_pathways.tsv")
  if (!is.null(gp) && nrow(gp)) {
    g1 <- unique(et[, .(event_key, human_genes, human_gene_ids)], by = "event_key")[, {
      sy <- strsplit(fcoalesce(human_genes, ""), ",", fixed = TRUE)[[1]]
      gi <- if (is.na(human_gene_ids)) character() else strsplit(human_gene_ids, ",", fixed = TRUE)[[1]]
      if (length(gi) != length(sy)) gi <- rep(NA_character_, length(sy))
      .(gid = gi, sym = sy)
    }, by = event_key]
    g1[gp, on = c(gid = "gene_id_h"), `:=`(dir = i.expected_direction, role = i.cgc_role, top = i.reactome_top_level)]
    g1[gp, on = c(sym = "symbol_h"), `:=`(dir = fcoalesce(dir, i.expected_direction), role = fcoalesce(role, i.cgc_role),
                                          top = fcoalesce(top, i.reactome_top_level))]
    agg <- g1[, .(role1 = role[1],
                  cgc_genes = paste(unique(na.omit(fifelse(is.na(role), NA_character_, paste0(sym, "(", role, ")")))), collapse = ";"),
                  reactome_top_level = paste(head(sort(unique(unlist(strsplit(na.omit(top), ";", fixed = TRUE)))), 8), collapse = ";")),
              by = event_key]
    et[agg, on = "event_key", `:=`(role1 = i.role1, cgc_genes = i.cgc_genes, reactome_top_level = i.reactome_top_level)]
  } else log("gene_pathways.tsv missing - run 05b for pathway columns", level = "WARNING")
  # direction of this event (not of the gene): event type, plus the CGC role for missense / promoter
  if (!"role1" %chin% names(et)) et[, role1 := NA_character_]
  et[, expected_direction := mapply(event_direction, event, role1, USE.NAMES = FALSE)]
  et[layer %chin% c("signature", "fusion"), expected_direction := NA_character_]
  et[, role1 := NULL]

  et[, event_id := sprintf("E%05d", seq_len(.N))]
  et[, c("mouse_1", "mouse_2", "mouse_3") := "pending"]
  lead <- c("event_id", "layer", "event", "human_locus", "human_genes", "human_gene_ids", "mouse_orthologs", "mouse_gene_ids", "homology",
            paste0("mouse_locus_", builds), "mapping_method", "mapped_fraction", "mapping_confidence", "hypothesis",
            "expected_direction", "cgc_genes", "reactome_top_level", "event_key", "mouse_1", "mouse_2", "mouse_3")
  for (c in setdiff(lead, names(et))) set(et, j = c, value = NA_character_)
  f <- write_tsv(et[, lead, with = FALSE], file.path(P$hp, "event_table.tsv"), log)

  s <- et[, .N, by = .(layer, mapping_confidence)][order(layer, mapping_confidence)]
  write_summary(P, "event_table", s[, .(category = paste0(layer, "_", mapping_confidence), n = N,
                                        counterpart = fcase(mapping_confidence == "high", "yes", mapping_confidence == "none", "no", default = "uncertain"))], log)
  write_provenance(P, cfg, "06b_build_event_table", outputs = f)
  log("done")
}

main()
