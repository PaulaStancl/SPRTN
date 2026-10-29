# ---------------------------------------------------------------------------
# 07_report.R - one self-contained HTML report of the human -> mouse expectation map.
#
#   Rscript scripts/07_report.R [--config ...] [--release ...]
#
# Sections: per layer, how many human events have a predicted mouse counterpart, how many are
# uncertain (1:many, multi-mapped, LOH in inbred mice, H2M class 2) and how many cannot have one
# (no ortholog, not modelable, unmapped) - from human_predictions/summary/*.tsv; H2M validation;
# human CN map (PURPLE segments) and expected mouse CN map per build (expected_mouse_cn.<build>.bed);
# circos plot human -> mouse (synteny blocks as grey ribbons, CN regions as coloured ribbons; needs the
# circlize package, skipped with a warning if it fails); CN liftOver
# cross-check and per-segment mapping confidence; the human events with no possible mouse
# counterpart (event_table.tsv, mapping_confidence none); ortholog differences between Ensembl
# releases; software and reference versions (provenance.json).
#
# Output  human_predictions/report.html
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))
suppressPackageStartupMessages(library(ggplot2))

LAYERS <- c(snv_coding = "Coding SNVs / indels (H2M)", snv_noncoding = "Non-coding variants (liftOver)",
            cn = "Copy number (gain / loss genes)", loh = "LOH genes", sv = "SV-disrupted genes",
            signatures = "SBS96 spectrum")

esc <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE) }

html_table <- function(dt, max_rows = 200) {
  if (is.null(dt) || !nrow(dt)) return("<p class='muted'>none</p>")
  more <- if (nrow(dt) > max_rows) sprintf("<p class='muted'>%d of %d rows shown - full table in the TSV</p>", max_rows, nrow(dt)) else ""
  d <- head(dt, max_rows)
  num <- vapply(d, is.numeric, TRUE)
  head_html <- paste0("<tr>", paste0("<th", ifelse(num, " class='n'", ""), ">", esc(names(d)), "</th>", collapse = ""), "</tr>")
  body <- vapply(seq_len(nrow(d)), function(i) paste0("<tr>", paste0("<td", ifelse(num, " class='n'", ""), ">",
                 vapply(d, function(col) { v <- col[i]; if (is.numeric(v) && !is.na(v) && v != round(v)) esc(signif(v, 3)) else esc(v) }, ""),
                 "</td>", collapse = ""), "</tr>"), "")
  paste0("<div class='tw'><table>", head_html, paste(body, collapse = ""), "</table></div>", more)
}

png_uri <- function(p, w, h) {
  f <- tempfile(fileext = ".png")
  ggsave(f, p, width = w, height = h, dpi = 110, bg = "white")
  base64enc::dataURI(file = f, mime = "image/png")
}

chrom_lengths <- function(M) {
  fai <- paste0(M$fasta, ".fai")
  if (file.exists(fai)) {
    d <- fread(fai, header = FALSE, select = 1:2, col.names = c("chrom", "len"))
    d[, chrom := norm_chrom(chrom)]
  } else {
    d <- load_genes(M$genes)[, .(len = max(end)), by = chrom]
  }
  d[chrom %chin% c(paste0("chr", 1:19), "chrX", "chrY")]
}

CN_COLS <- c(gain = "#c0392b", loss = "#2471a3", LOH = "#b7950b")

human_lengths <- function(cfg, P) {
  fai <- paste0(cfg$human_fasta, ".fai")
  d <- if (file.exists(fai)) fread(fai, header = FALSE, select = 1:2, col.names = c("chrom", "len"))
       else load_genes(P$human_genes)[, .(len = max(end)), by = chrom]
  d[, chrom := norm_chrom(chrom)][chrom %chin% c(paste0("chr", 1:22), "chrX", "chrY")]
}

# d: chrom start end direction (gain / loss / LOH); LOH drawn as a second, thinner track
ideogram <- function(d, lens, title) {
  lens <- copy(lens)[, chrom := factor(chrom, levels = rev(lens$chrom[order(chrom_order(lens$chrom))]))]
  b <- copy(d)[chrom %chin% as.character(lens$chrom)]
  b[, chrom := factor(chrom, levels = levels(lens$chrom))]
  b[, `:=`(y0 = as.integer(chrom) + fifelse(direction == "LOH", 0.22, -0.18), y1 = as.integer(chrom) + fifelse(direction == "LOH", 0.36, 0.18))]
  ggplot() +
    geom_rect(data = lens, aes(xmin = 0, xmax = len / 1e6, ymin = as.integer(chrom) - 0.18, ymax = as.integer(chrom) + 0.18),
              fill = "grey92", colour = "grey70", linewidth = 0.2) +
    geom_rect(data = b, aes(xmin = start / 1e6, xmax = pmax(end, start + 2e5) / 1e6, ymin = y0, ymax = y1, fill = direction)) +
    scale_fill_manual(values = CN_COLS, drop = FALSE, name = NULL) +
    scale_y_continuous(breaks = seq_along(levels(lens$chrom)), labels = levels(lens$chrom), expand = c(0.01, 0.01)) +
    labs(x = "Mb", y = NULL, title = title) +
    theme_minimal(base_size = 10) + theme(panel.grid.major.y = element_blank(), legend.position = "top")
}

# human PURPLE segments (non-neutral + LOH) from cn_segment_mapping.tsv
human_segments <- function(segm) {
  s <- unique(segm, by = "segment")
  s[, c("chrom", "start", "end") := tstrsplit(segment, "[:-]", type.convert = TRUE)]
  rbind(s[direction %chin% c("gain", "loss"), .(chrom, start, end, direction)],
        s[truthy(loh), .(chrom, start, end, direction = "LOH")])
}

# circos: human chr1-22,X (clockwise) facing mouse chr1-19,X (mirrored); outer ring chromosomes, inner ring
# CN (human PURPLE segments / expected mouse intervals; LOH inner half); grey ribbons = synteny blocks with
# >= min_genes 1:1 orthologs, coloured ribbons = human CN region -> expected mouse interval
circos_png <- function(f, hseg, bed, sb, hlens, mlens, build, min_genes) {
  hl <- hlens[chrom %chin% c(paste0("chr", 1:22), "chrX")][order(chrom_order(chrom))]
  ml <- mlens[chrom %chin% c(paste0("chr", 1:19), "chrX")][order(-chrom_order(chrom))]
  hs <- function(ch) paste0("h", sub("^chr", "", ch)); ms <- function(ch) paste0("m", sub("^chr", "", ch))
  sec <- data.frame(chr = factor(c(hs(hl$chrom), ms(ml$chrom)), levels = c(hs(hl$chrom), ms(ml$chrom))),
                    start = 0, end = c(hl$len, ml$len))
  args <- list(f, width = 2200, height = 2200, res = 250)
  if (capabilities("cairo")) args$type <- "cairo"
  do.call(png, args)
  on.exit({ circlize::circos.clear(); dev.off() }, add = TRUE)
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.after = c(rep(0.6, nrow(hl) - 1), 6, rep(0.6, nrow(ml) - 1), 6),
                       cell.padding = c(0, 0, 0, 0), track.margin = c(0.004, 0.004))
  circlize::circos.genomicInitialize(sec, plotType = NULL)
  circlize::circos.track(ylim = c(0, 1), track.height = 0.06, bg.border = NA, panel.fun = function(x, y) {
    s <- circlize::CELL_META$sector.index; xl <- circlize::CELL_META$xlim
    circlize::circos.rect(xl[1], 0, xl[2], 1, col = if (startsWith(s, "h")) "grey75" else "#d8c3a0", border = NA)
    circlize::circos.text(mean(xl), 0.5, sub("^[hm]", "", s), cex = 0.5, facing = "clockwise", niceFacing = TRUE)
  })
  cn <- rbind(hseg[, .(sector = hs(chrom), start, end, direction)], bed[, .(sector = ms(chrom), start, end, direction)])
  cn <- cn[sector %chin% levels(sec$chr)]
  circlize::circos.track(ylim = c(0, 2), track.height = 0.07, bg.border = "grey85", panel.fun = function(x, y) {
    d <- cn[sector == circlize::CELL_META$sector.index]
    if (nrow(d)) circlize::circos.rect(d$start, fifelse(d$direction == "LOH", 0, 1), pmax(d$end, d$start + 3e5),
                                       fifelse(d$direction == "LOH", 1, 2), col = CN_COLS[d$direction], border = NA)
  })
  bg <- sb[n_genes >= min_genes & h_chrom %chin% hl$chrom & m_chrom %chin% ml$chrom]
  if (nrow(bg)) circlize::circos.genomicLink(bg[, .(hs(h_chrom), h_start, h_end)], bg[, .(ms(m_chrom), m_start, m_end)],
                                             col = adjustcolor("grey55", 0.18), border = NA)
  lk <- bed[!is.na(human_locus)]
  lk[, c("hc", "hs0", "he0") := tstrsplit(human_locus, "[:-]", type.convert = TRUE)]
  lk <- lk[hc %chin% hl$chrom & chrom %chin% ml$chrom]
  if (nrow(lk)) circlize::circos.genomicLink(lk[, .(hs(hc), hs0, he0)], lk[, .(ms(chrom), start, end)],
                                             col = adjustcolor(CN_COLS[lk$direction], 0.7), border = NA)
  title(sprintf("Human (grey, GRCh38) -> mouse (tan, %s)", build), cex.main = 0.9)
  legend("bottomleft", fill = c(CN_COLS, adjustcolor("grey55", 0.4)), border = NA, bty = "n", cex = 0.75,
         legend = c("gain", "loss", "LOH", sprintf("synteny block (>= %d genes)", min_genes)))
  invisible(f)
}

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("07_report", P)
  hb <- cfg$primary_mouse_build      # main build for counterpart calls, event table, report

  sm <- rbindlist(lapply(list.files(P$summary, "\\.tsv$", full.names = TRUE), fread), use.names = TRUE, fill = TRUE)
  ov <- sm[layer %chin% names(LAYERS) & counterpart %chin% c("yes", "uncertain", "no"),
           .(n = sum(n)), by = .(layer, counterpart)]
  ov <- dcast(ov, layer ~ counterpart, value.var = "n", fill = 0)
  for (c in c("yes", "uncertain", "no")) if (!c %chin% names(ov)) ov[, (c) := 0L]
  ov <- ov[match(intersect(names(LAYERS), layer), layer)]
  ov[, `:=`(Layer = LAYERS[layer], total = yes + uncertain + no)]
  ov <- ov[, .(Layer, `human events` = total, `predicted mouse counterpart` = yes, `uncertain / not assessable` = uncertain,
               `no possible counterpart` = no)]

  sec <- list()
  sec$overview <- paste0("<h2>Mouse counterpart per layer</h2>",
    "<p>Counts are human events (variants, genes, SNVs in the spectrum). <b>No possible counterpart</b> = no mouse ortholog, ",
    "not modelable by H2M, outside the mouse build or not lifting over; it is its own category, not a discordance. ",
    "<b>Uncertain</b> = 1:many orthologs, multi-mapped, H2M class 2, or LOH (not assessable in inbred mice).</p>", html_table(ov))
  sec$layers <- paste0("<h2>Per-layer detail</h2>", paste(vapply(split(sm, by = "layer"), function(d)
    paste0("<h3>", esc(d$layer[1]), "</h3>", html_table(d[, .(category, n, counterpart)])), ""), collapse = ""))

  val <- read_layer(P, "h2m_validation.tsv")
  sec$val <- paste0("<h2>H2M validation (known hotspots, run before the patient's variants)</h2>",
                    html_table(if (is.null(val)) NULL else val[, intersect(c("gene", "protein", "expected_mouse", "hgvsp_m_h2m", "h2m_class", "ref_in_fasta", "pass"), names(val)), with = FALSE]))

  segm <- read_layer(P, "cn_segment_mapping.tsv")
  hseg <- if (is.null(segm)) data.table(chrom = character(), start = numeric(), end = numeric(), direction = character()) else human_segments(segm[build == hb])
  hlens <- human_lengths(cfg, P)
  plots <- sprintf("<img alt='human CN' src='%s'/>", png_uri(ideogram(hseg, hlens, "Human tumour CN (PURPLE segments, GRCh38)"), 9, 6.5))
  beds <- list()
  for (b in names(P$mouse)) {
    f <- file.path(P$hp, sprintf("expected_mouse_cn.%s.bed", b))
    if (!file.exists(f)) next
    bed <- fread(f)
    setnames(bed, "#chrom", "chrom")
    beds[[b]] <- bed
    plots <- c(plots, sprintf("<img alt='expected mouse CN %s' src='%s'/>", b, png_uri(ideogram(bed, chrom_lengths(P$mouse[[b]]), paste("Expected mouse CN -", b)), 9, 6)),
               html_table(bed[, .(interval = paste0(chrom, ":", start + 1, "-", end), direction, states, n_genes, confidence,
                                 human_locus = if ("human_locus" %chin% names(bed)) human_locus else NA_character_,
                                 genes = substr(genes, 1, 120))], 60))
  }
  circ <- "<p class='muted'>circos plot not drawn</p>"
  sbf <- file.path(P$hp, sprintf("synteny_blocks.%s.tsv", hb))
  if (!is.null(beds[[hb]]) && file.exists(sbf) && "human_locus" %chin% names(beds[[hb]])) {
    circ <- tryCatch({
      f <- tempfile(fileext = ".png")
      circos_png(f, hseg, beds[[hb]], fread(sbf), hlens, chrom_lengths(P$mouse[[hb]]), hb, cfg$report$circos_min_block_genes)
      sprintf("<img alt='circos human to mouse' src='%s'/>", base64enc::dataURI(file = f, mime = "image/png"))
    }, error = function(e) { log("circos plot failed: %s", conditionMessage(e), level = "WARNING"); sprintf("<p class='muted'>circos plot failed: %s</p>", esc(conditionMessage(e))) })
  } else log("circos plot skipped: rerun 03 (needs synteny_blocks.%s.tsv and human_locus in the BED)", hb, level = "WARNING")
  sec$cn <- paste0("<h2>Copy number: human tumour and expected mouse</h2>",
                   "<p>Human: PURPLE segments, gain / loss relative to ploidy, LOH as the upper thin track. ",
                   "Mouse: the same events moved through 1:1 orthologs and merged into intervals per build.</p>",
                   paste(plots, collapse = ""),
                   "<h3>Human -> mouse (", esc(hb), ")</h3>",
                   "<p>Outer ring: human (grey) and mouse (tan) chromosomes. Inner ring: CN (outer half gain / loss, inner half LOH). ",
                   "Grey ribbons: conserved synteny blocks. Coloured ribbons: where each human CN region is expected in mouse.</p>", circ,
                   "<h3>Direct liftOver cross-check</h3>", html_table(read_layer(P, "cn_liftover_crosscheck.tsv")),
                   "<h3>Mapping confidence of human segments (", esc(hb), ")</h3>",
                   html_table(if (is.null(segm)) NULL else segm[build == hb, .N, by = .(direction, confidence)][order(direction, confidence)]),
                   "<h3>Human segments split across mouse chromosomes</h3>",
                   html_table(if (is.null(segm)) NULL else segm[build == hb & split_across_chromosomes == TRUE & direction != "neutral",
                                                                 .(segment, state, n_genes, mouse_chroms, n_blocks, confidence)]))

  et <- read_layer(P, "event_table.tsv")
  sec$none <- paste0("<h2>Human events with no possible mouse counterpart</h2>",
                     html_table(if (is.null(et)) NULL else et[mapping_confidence == "none", .(event_id, layer, event, human_locus, human_genes, hypothesis)], 300))
  sec$fus <- paste0("<h2>Fusions (feasibility flag only - no mouse fusion is predicted)</h2>", html_table(read_layer(P, "expected_mouse_fusions.tsv")))
  sig <- read_layer(P, "expected_signatures.tsv")
  sec$sig <- paste0("<h2>Signatures (SIGS) and expected mouse exposure</h2>", html_table(sig))
  dfiles <- list.files(P$hp, "^ortholog_release_diff_.*\\.tsv$", full.names = TRUE)
  sec$diff <- paste0("<h2>Ortholog differences between Ensembl releases</h2>",
                     if (length(dfiles)) paste(vapply(dfiles, function(f) paste0("<h3>", esc(basename(f)), "</h3>", html_table(fread(f)[, .N, by = change])), ""), collapse = "")
                     else "<p class='muted'>not run (ensembl_release_check null)</p>")

  prov <- if (file.exists(P$provenance)) fromJSON(P$provenance, simplifyVector = FALSE) else list()
  sw <- if (length(prov$scripts)) prov$scripts[[length(prov$scripts)]]$software else software_versions()
  swt <- data.table(software = names(sw), version = vapply(sw, function(x) if (is.null(x)) "" else as.character(x), ""))
  rv <- prov$scripts$`00_resources`$resource_versions
  rvt <- if (length(rv)) data.table(file = names(rv), url = vapply(rv, function(x) fcoalesce(x$url, ""), ""),
                                    downloaded = vapply(rv, function(x) fcoalesce(x$downloaded, ""), "")) else NULL
  pw <- prov$scripts$`05b_pathways`
  sec$ver <- paste0("<h2>Software and references</h2>",
                    sprintf("<p>Ensembl release <b>%s</b> (%s); H2M mouse build %s; Reactome %s; COSMIC CGC %s; seed %s.</p>",
                            cfg$ensembl_release, esc(cfg$ensembl_release_tag), esc(hb), esc(fcoalesce(as.character(pw$reactome_release), "n/a")),
                            esc(fcoalesce(as.character(pw$cosmic_cgc_version), "not set")), esc(cfg$seed)),
                    html_table(swt), "<h3>Reference files</h3>", html_table(rvt))

  css <- "
:root{--fg:#1f2328;--bg:#ffffff;--mut:#59636e;--line:#d1d9e0;--head:#f6f8fa;color-scheme:light}
@media (prefers-color-scheme: dark){:root{--fg:#e6edf3;--bg:#0d1117;--mut:#9198a1;--line:#3d444d;--head:#151b23;color-scheme:dark}}
body{font:14px/1.45 system-ui,sans-serif;color:var(--fg);background:var(--bg);margin:0;padding:24px 16px;max-width:1200px;margin-inline:auto}
h1{font-size:22px}h2{font-size:18px;margin-top:32px;border-bottom:1px solid var(--line);padding-bottom:4px}h3{font-size:15px}
.tw{overflow-x:auto}table{border-collapse:collapse;font-size:12.5px;margin:8px 0}th,td{border:1px solid var(--line);padding:3px 7px;text-align:left;vertical-align:top}
th{background:var(--head)}td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}.muted{color:var(--mut)}img{max-width:100%;background:#fff}"
  html <- paste0("<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
                 "<title>", esc(cfg$patient), " human-to-mouse expectation map</title><style>", css, "</style></head><body>",
                 "<h1>", esc(cfg$patient), " (", esc(cfg$tumour_id), ") - expected mouse alterations</h1>",
                 sprintf("<p class='muted'>Generated %s from %s. Predictions are fixed before mouse data; mouse columns are pending.</p>",
                         format(Sys.time(), "%Y-%m-%d %H:%M"), esc(cfg$oncoanalyser_dir)),
                 paste(unlist(sec), collapse = "\n"), "</body></html>")
  f <- file.path(P$hp, "report.html")
  writeLines(html, f)
  log("wrote %s", f)
  write_provenance(P, cfg, "07_report", outputs = f)
  log("done")
}

main()
