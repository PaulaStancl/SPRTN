# ---------------------------------------------------------------------------
# 06_build_matrix.R - gene x sample matrix: the human column filled now, mouse_1..3 empty.
#
#   Rscript scripts/06_build_matrix.R [--config ...] [--release ...]
#
# Rows: every gene PURPLE assessed for copy number (purple.cnv.gene.tsv, via 03) plus every gene
# with an event in 01-04. Cells hold several events, ';'-separated, from the codes SNV, INDEL,
# AMP, GAIN, DEL, LOH, SV, FUSION. A cell without events is
#   neutral          the gene was assessed (in PURPLE's gene table, CN assessable) and has none
#   not_assessable   it was not assessed
# expected_mouse: the human codes when the gene has a 1:1 ortholog placed in the H2M build;
# "one2many:<codes>" for 1:many / many:many; no_ortholog / not_in_release / not_in_build when the
# event has no possible mouse counterpart (its own category, never a discordance); neutral /
# not_assessable as for the human column. LOH is predicted as LOH, but a later mouse call of LOH
# is only possible where the mouse genome is heterozygous (inbred: usually not_assessable).
# freq_HCC / freq_HB: from config tumour_frequency_table (gene_symbol, freq_HCC, freq_HB) when set.
#
# Output  human_predictions/
#   gene_matrix.tsv   gene_id_h symbol_h chrom_h homology gene_id_m symbol_m human human_confidence
#                     human_detail expected_mouse mouse_1 mouse_2 mouse_3 freq_HCC freq_HB
#   summary/matrix.tsv
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("06_build_matrix", P)
  hb <- cfg$h2m_mouse_build

  ev <- collect_gene_events(P)[!is.na(gene_id_h)]
  cn <- read_layer(P, "expected_mouse_cn_genes.tsv")
  if (is.null(cn)) stop("expected_mouse_cn_genes.tsv missing - run 03 first", call. = FALSE)
  hg <- load_genes(P$human_genes)
  M <- P$mouse[[hb]]
  mg <- load_genes(M$genes)

  genes <- unique(rbind(cn[!is.na(gene_id_h), .(gene_id_h, assessed = state != "not_assessable")],
                        ev[, .(gene_id_h, assessed = TRUE)]), by = "gene_id_h")
  genes[cn[!is.na(gene_id_h)], on = "gene_id_h", assessed := i.state != "not_assessable"]
  cls <- classify_orthologs(read_orthologs(P$orthologs), genes$gene_id_h)
  genes[cls, on = "gene_id_h", `:=`(homology = i.homology, gene_id_m = i.gene_id_m)]
  genes[hg, on = c(gene_id_h = "gene_id"), `:=`(symbol_h = i.gene_name, chrom_h = i.chrom)]
  genes[mg, on = c(gene_id_m = "gene_id"), symbol_m := i.gene_name]
  genes[, in_build := !is.na(gene_id_m) & gene_id_m %chin% mg$gene_id]

  per_gene <- ev[, .(human = paste(sort(unique(matrix_code)), collapse = ";"), human_confidence = best_confidence(confidence),
                     human_detail = paste(unique(paste0(layer, ":", event, " ", detail)), collapse = " | "),
                     ev_symbol = symbol_h[1]), by = gene_id_h]
  genes[per_gene, on = "gene_id_h", `:=`(human = i.human, human_confidence = i.human_confidence,
                                         human_detail = i.human_detail, symbol_h = fcoalesce(symbol_h, i.ev_symbol))]
  genes[is.na(human), human := fifelse(assessed, "neutral", "not_assessable")]
  has_ev <- !genes$human %chin% c("neutral", "not_assessable")
  genes[, expected_mouse := fcase(
    !has_ev, human,
    homology %chin% c("no_ortholog", "not_in_release"), homology,
    homology %chin% c("ortholog_one2many", "ortholog_many2many"), paste0("one2many:", human),
    !in_build, "not_in_build",
    default = human)]
  genes[, c("mouse_1", "mouse_2", "mouse_3") := NA_character_]

  ff <- cfg$tumour_frequency_table
  genes[, c("freq_HCC", "freq_HB") := NA_real_]
  if (!is.null(ff)) {
    fq <- fread(ff)
    need_cols(fq, c("gene_symbol", "freq_HCC", "freq_HB"), "tumour_frequency_table")
    genes[fq, on = c(symbol_h = "gene_symbol"), `:=`(freq_HCC = i.freq_HCC, freq_HB = i.freq_HB)]
  }
  out <- genes[order(!has_ev, chrom_order(fcoalesce(chrom_h, "chrUn")), symbol_h),
               .(gene_id_h, symbol_h, chrom_h, homology, gene_id_m, symbol_m, human, human_confidence, human_detail,
                 expected_mouse, mouse_1, mouse_2, mouse_3, freq_HCC, freq_HB)]
  f <- write_tsv(out, file.path(P$hp, "gene_matrix.tsv"), log)

  e <- out[!human %chin% c("neutral", "not_assessable")]
  cp <- fcase(e$expected_mouse %chin% c("no_ortholog", "not_in_release", "not_in_build"), "no",
              startsWith(e$expected_mouse, "one2many:"), "uncertain", default = "yes")
  rows <- rbind(srow("genes_in_matrix", nrow(out), "info"), srow("genes_neutral", out[human == "neutral", .N], "info"),
                srow("genes_not_assessable", out[human == "not_assessable", .N], "info"),
                srow("genes_with_event", nrow(e), "info"),
                data.table(category = paste0("event_gene_", fcase(cp == "yes", "counterpart", cp == "uncertain", "one2many",
                                                                  default = e$expected_mouse)), counterpart = cp)[, .(n = .N), by = .(category, counterpart)][, .(category, n, counterpart)])
  write_summary(P, "matrix", rows, log)
  write_provenance(P, cfg, "06_build_matrix", inputs = list(tumour_frequency_table = ff), outputs = f)
  log("done")
}

main()
