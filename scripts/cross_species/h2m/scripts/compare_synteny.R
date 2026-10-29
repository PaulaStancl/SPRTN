# ---------------------------------------------------------------------------
# compare_synteny.R - our gene-order synteny blocks (03) vs UCSC's human-mouse alignment chains.
# Run after 03.
#
#   Rscript scripts/compare_synteny.R [--config ...] [--release ...] [--build GRCm38] [--min-chain-kb 100]
#                                     [--plot-chroms chr9,chr1]
#
# Reference: the UCSC liftOver chain hg38To<Mm39|Mm10>.over.chain.gz (00_resources), i.e. UCSC's net-derived
# human (target) - mouse (query) whole-genome alignment chains. Chains with a human span >= --min-chain-kb
# (default 100 kb) are the reference syntenic segments. Ensembl's synteny regions are not compared: they are
# only in the Compara MySQL database, not in a download file.
#
# Comparisons
#   1. block support    per block: fraction of its human span covered by reference chains on the same mouse
#                       chromosome whose mouse interval overlaps the block's mouse interval
#   2. gene agreement   every 1:1 ortholog's human gene midpoint is lifted with the chain (liftover.min_match);
#                       agrees if it lands on the ortholog's mouse chromosome within the gene +/- 50 kb
#   3. boundaries       per internal block boundary on a human chromosome: distance to the nearest reference
#                       chain boundary (start or end)
#   4. orientation      block orientation (+ / -) vs the strand of the reference chain covering most of it
#
# Input   human_predictions/synteny_blocks.<build>.tsv (03), orthologs, gene tables, the chain file
# Output  human_predictions/synteny_compare.<build>.tsv          per block: the four measures
#         human_predictions/synteny_compare_genes.<build>.tsv    per 1:1 ortholog: lifted position, agreement
#         human_predictions/synteny_compare_summary.<build>.tsv  overall numbers
#         human_predictions/synteny_compare.<build>.html         summary + dot plots (our blocks over UCSC chains)
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))
suppressPackageStartupMessages(library(ggplot2))

# chain headers of a UCSC chain file: human (t) and mouse (q) intervals, mouse coordinates on the + strand
read_chains <- function(f) {
  ch <- fread(cmd = sprintf("gzip -dc %s | grep '^chain'", shQuote(f)), header = FALSE, sep = " ",
              col.names = c("x", "score", "tName", "tSize", "tStrand", "tStart", "tEnd", "qName", "qSize", "qStrand", "qStart", "qEnd", "id"))
  ch[, `:=`(h_chrom = tName, h_start = tStart, h_end = tEnd, m_chrom = qName, strand = qStrand,
            m_start = fifelse(qStrand == "-", qSize - qEnd, qStart), m_end = fifelse(qStrand == "-", qSize - qStart, qEnd))]
  ch[, .(id, score, h_chrom, h_start, h_end, m_chrom, m_start, m_end, strand, h_span = h_end - h_start)]
}

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("compare_synteny", P)
  build <- if (is.null(a$build)) cfg$primary_mouse_build else a$build
  M <- P$mouse[[build]]
  min_kb <- as.numeric(if (is.null(a$min_chain_kb)) 100 else a$min_chain_kb)
  plot_chroms <- if (is.null(a$plot_chroms)) c("chr1", "chr9") else strsplit(a$plot_chroms, ",")[[1]]

  sbf <- file.path(P$hp, sprintf("synteny_blocks.%s.tsv", build))
  if (!file.exists(sbf)) stop("missing ", sbf, " - run 03 first", call. = FALSE)
  sb <- fread(sbf)
  ch <- read_chains(M$chain)
  ref <- ch[h_span >= min_kb * 1e3 & h_chrom %chin% STD_CHROMS]
  log("%s: %d blocks; %d chains, %d with human span >= %g kb (%.0f%% of chained human bp)", build, nrow(sb), nrow(ch), nrow(ref),
      min_kb, 100 * sum(ref$h_span) / sum(ch$h_span))

  # ---- 1 + 4: block support and orientation ------------------------------------------------------
  b <- sb[, .(block, h_chrom, s = h_start, e = h_end, m_chrom, m_start, m_end, n_genes, orientation)]
  r <- ref[, .(id, h_chrom, s = h_start + 1, e = h_end, rm_chrom = m_chrom, rm_start = m_start + 1, rm_end = m_end, strand)]
  setkey(r, h_chrom, s, e)
  ov <- foverlaps(b, r, by.x = c("h_chrom", "s", "e"), nomatch = NULL)
  ov[, ov_bp := pmax(0, pmin(e, i.e) - pmax(s, i.s) + 1)]
  ov[, consistent := rm_chrom == m_chrom & rm_start <= m_end & rm_end >= m_start]
  sup <- ov[, .(supported_bp = sum(ov_bp[consistent]), any_chain_bp = sum(ov_bp),
                main_strand = if (any(consistent)) strand[consistent][which.max(ov_bp[consistent])] else NA_character_), by = block]
  res <- sup[sb, on = "block"]
  res[, h_span := h_end - h_start + 1]
  res[is.na(supported_bp), `:=`(supported_bp = 0, any_chain_bp = 0)]
  res[, support_frac := pmin(1, supported_bp / h_span)]
  res[, orientation_agrees := fifelse(is.na(main_strand) | n_genes < 2, NA, main_strand == orientation)]

  # ---- 3: boundaries -----------------------------------------------------------------------------
  bd <- sb[order(h_chrom, h_start)][, .(pos = c(h_start, h_end), block = block), by = h_chrom]
  bd <- bd[, .SD[!(pos == min(pos) | pos == max(pos))], by = h_chrom]      # internal boundaries only
  rb <- ref[, .(pos = c(h_start, h_end)), by = h_chrom]
  bd[, dist_to_chain_boundary := vapply(seq_len(.N), function(i) {
    p <- rb$pos[rb$h_chrom == h_chrom[i]]; if (!length(p)) NA_real_ else min(abs(p - pos[i]))
  }, 0)]
  if (nrow(bd)) res[bd[, .(d = if (all(is.na(dist_to_chain_boundary))) NA_real_ else min(dist_to_chain_boundary, na.rm = TRUE)), by = block],
                    on = "block", boundary_dist_bp := i.d]
  if (!"boundary_dist_bp" %chin% names(res)) res[, boundary_dist_bp := NA_real_]

  # ---- 2: gene agreement -------------------------------------------------------------------------
  hg <- load_genes(P$human_genes); mg <- load_genes(M$genes)
  o1 <- classify_orthologs(read_orthologs(P$orthologs))[homology == "ortholog_one2one", .(gene_id_h, gene_id_m)]
  o1[hg, on = c(gene_id_h = "gene_id"), `:=`(h_chrom = i.chrom, mid = (i.start + i.end) %/% 2)]
  o1[mg, on = c(gene_id_m = "gene_id"), `:=`(m_chrom = i.chrom, m_start = i.start, m_end = i.end)]
  o1 <- o1[h_chrom %chin% STD_CHROMS & !is.na(m_chrom)]
  lo <- run_liftover(o1[, .(chrom = h_chrom, start = mid - 1L, end = mid, name = gene_id_h)], M$chain, cfg$liftover$min_match,
                     P$tmp, paste0("synteny_genes_", build))
  lm <- lo$mapped[!duplicated(name)]
  o1[lm, on = c(gene_id_h = "name"), `:=`(lift_chrom = i.chrom, lift_pos = i.end)]
  o1[, agrees := fifelse(is.na(lift_chrom), NA, lift_chrom == m_chrom & lift_pos >= m_start - 5e4 & lift_pos <= m_end + 5e4)]
  o1[, status := fcase(is.na(agrees), "not_lifted", agrees, "agrees", lift_chrom != m_chrom, "other_chromosome", default = "same_chromosome_far")]
  # block of each gene from 03 (synteny_blocks has extents only): genes whose mouse position lies in a block
  gk <- o1[, .(gene_id_h, m_chrom, s = m_start, e = m_start)]
  sk <- sb[, .(block, m_chrom, s = m_start, e = m_end)]; setkey(sk, m_chrom, s, e)
  gb <- foverlaps(gk, sk, by.x = c("m_chrom", "s", "e"), nomatch = NULL)[!duplicated(gene_id_h), .(gene_id_h, block)]
  o1[gb, on = "gene_id_h", block := i.block]
  res[o1[!is.na(block), .(gene_agree_frac = mean(agrees, na.rm = TRUE), genes_lifted = sum(!is.na(agrees))), by = block],
      on = "block", `:=`(gene_agree_frac = i.gene_agree_frac, genes_lifted = i.genes_lifted)]

  res <- res[order(chrom_order(h_chrom), h_start)]
  f1 <- write_tsv(res[, .(block, h_chrom, h_start, h_end, m_chrom, m_start, m_end, n_genes, orientation, support_frac,
                          main_strand, orientation_agrees, boundary_dist_bp, gene_agree_frac, genes_lifted)],
                  file.path(P$hp, sprintf("synteny_compare.%s.tsv", build)), log)
  f2 <- write_tsv(o1[, .(gene_id_h, gene_id_m, h_chrom, mid, m_chrom, m_start, m_end, lift_chrom, lift_pos, status, block)],
                  file.path(P$hp, sprintf("synteny_compare_genes.%s.tsv", build)), log)

  big <- res[n_genes >= 5]
  sm <- data.table(
    measure = c("blocks", "blocks with >= 5 genes", "human bp in blocks", "reference chains (>= min kb)",
                "block support: median fraction of human span on a consistent chain (blocks >= 5 genes)",
                "block support: blocks >= 5 genes with support >= 0.5",
                "orientation agrees with main chain (blocks >= 5 genes)",
                "internal boundaries: median distance to a chain boundary (kb)",
                "internal boundaries within 100 kb of a chain boundary",
                "1:1 orthologs lifted", "lifted orthologs landing on their mouse ortholog (+/- 50 kb)",
                "lifted orthologs on another mouse chromosome"),
    value = c(nrow(res), nrow(big), sum(res$h_span), nrow(ref),
              round(median(big$support_frac), 3), sprintf("%d of %d", sum(big$support_frac >= 0.5), nrow(big)),
              sprintf("%d of %d", sum(big$orientation_agrees, na.rm = TRUE), sum(!is.na(big$orientation_agrees))),
              round(median(bd$dist_to_chain_boundary, na.rm = TRUE) / 1e3, 1),
              sprintf("%d of %d", sum(bd$dist_to_chain_boundary <= 1e5, na.rm = TRUE), sum(!is.na(bd$dist_to_chain_boundary))),
              sum(!is.na(o1$agrees)), sprintf("%.1f%%", 100 * mean(o1$agrees, na.rm = TRUE)),
              sprintf("%.1f%%", 100 * mean(o1$status == "other_chromosome"))))
  sm[, value := as.character(value)]
  f3 <- write_tsv(sm, file.path(P$hp, sprintf("synteny_compare_summary.%s.tsv", build)), log)
  for (i in seq_len(nrow(sm))) log("  %s: %s", sm$measure[i], sm$value[i])

  # ---- dot plots: UCSC chains (thin, coloured by mouse chromosome) with our blocks (black) on top -----------
  imgs <- vapply(intersect(plot_chroms, unique(sb$h_chrom)), function(hc) {
    rc <- ref[h_chrom == hc]; ob <- sb[h_chrom == hc & n_genes >= 2]
    seg <- function(d) d[, .(x = h_start / 1e6, xend = h_end / 1e6, y = fifelse(strand == "-", m_end, m_start) / 1e6,
                             yend = fifelse(strand == "-", m_start, m_end) / 1e6, m_chrom)]
    ob[, strand := orientation]
    p <- ggplot() +
      geom_segment(data = seg(rc), aes(x = x, xend = xend, y = y, yend = yend, colour = m_chrom), linewidth = 1.2, alpha = 0.6) +
      geom_segment(data = seg(ob), aes(x = x, xend = xend, y = y, yend = yend), colour = "black", linewidth = 0.5) +
      facet_wrap(~ m_chrom, scales = "free_y") +
      labs(x = paste(hc, "(human, Mb)"), y = "mouse position (Mb)", colour = "UCSC chain on",
           title = sprintf("%s: UCSC chains >= %g kb (colour) and our synteny blocks (black), mouse %s", hc, min_kb, build)) +
      theme_minimal(base_size = 9) + theme(legend.position = "none")
    f <- tempfile(fileext = ".png"); ggsave(f, p, width = 10, height = 6, dpi = 110, bg = "white")
    sprintf("<img style='max-width:100%%' src='%s'/>", base64enc::dataURI(file = f, mime = "image/png"))
  }, "")
  tab <- paste0("<table>", paste0("<tr><td>", sm$measure, "</td><td><b>", sm$value, "</b></td></tr>", collapse = ""), "</table>")
  fh <- file.path(P$hp, sprintf("synteny_compare.%s.html", build))
  writeLines(c("<!doctype html><html><head><meta charset='utf-8'><title>Synteny blocks vs UCSC chains</title>",
               "<style>body{font:14px system-ui,sans-serif;max-width:1100px;margin:auto;padding:16px}td{border:1px solid #ccc;padding:3px 8px}table{border-collapse:collapse}</style></head><body>",
               sprintf("<h1>Our synteny blocks vs UCSC human-mouse chains (%s)</h1>", build),
               "<p>Ours: runs of 1:1 orthologs in mouse order whose human gene rank moves by &le; synteny.max_rank_jump (03). ",
               sprintf("Reference: UCSC liftOver chains (net-derived whole-genome alignment) with a human span &ge; %g kb.</p>", min_kb),
               tab, "<h2>Dot plots</h2><p>One panel per mouse chromosome. A black segment lying on a coloured one = our block agrees with the alignment.</p>",
               imgs, "</body></html>"), fh)
  log("wrote %s", fh)
  write_provenance(P, cfg, "compare_synteny", resources = list(chain = M$chain), outputs = c(f1, f2, f3, fh),
                   extra = list(min_chain_kb = min_kb, build = build))
  log("done")
}

main()
