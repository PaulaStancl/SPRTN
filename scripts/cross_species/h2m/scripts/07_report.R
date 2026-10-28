# ---------------------------------------------------------------------------
# 07_report.R - one self-contained HTML report of the human -> mouse expectation map.
#
#   Rscript scripts/07_report.R [--config ...] [--release ...]
#
# Sections: per layer, how many human events have a predicted mouse counterpart, how many are
# uncertain (1:many, multi-mapped, LOH in inbred mice, H2M class 2) and how many cannot have one
# (no ortholog, not modelable, unmapped) - from human_predictions/summary/*.tsv; H2M validation;
# expected mouse CN map (ideogram per build, from expected_mouse_cn.<build>.bed); CN liftOver
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

ideogram <- function(bed, lens, build) {
  lens <- copy(lens)[, chrom := factor(chrom, levels = rev(lens$chrom[order(chrom_order(lens$chrom))]))]
  b <- copy(bed)
  setnames(b, "#chrom", "chrom")
  b <- b[chrom %chin% as.character(lens$chrom)]
  b[, chrom := factor(chrom, levels = levels(lens$chrom))]
  b[, track := fifelse(direction == "LOH", 0.28, 0)]
  ggplot() +
    geom_rect(data = lens, aes(xmin = 0, xmax = len / 1e6, ymin = as.integer(chrom) - 0.18, ymax = as.integer(chrom) + 0.18),
              fill = "grey92", colour = "grey70", linewidth = 0.2) +
    geom_rect(data = b, aes(xmin = start / 1e6, xmax = pmax(end, start + 2e5) / 1e6,
                            ymin = as.integer(chrom) - 0.18 + track, ymax = as.integer(chrom) + 0.18 + track, fill = direction)) +
    scale_fill_manual(values = c(gain = "#c0392b", loss = "#2471a3", LOH = "#7d6608"), drop = FALSE, name = NULL) +
    scale_y_continuous(breaks = seq_along(levels(lens$chrom)), labels = levels(lens$chrom), expand = c(0.01, 0.01)) +
    labs(x = "Mb", y = NULL, title = paste("Expected mouse CN -", build)) +
    theme_minimal(base_size = 10) + theme(panel.grid.major.y = element_blank(), legend.position = "top")
}

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("07_report", P)
  hb <- cfg$h2m_mouse_build

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

  plots <- character()
  for (b in names(P$mouse)) {
    f <- file.path(P$hp, sprintf("expected_mouse_cn.%s.bed", b))
    if (!file.exists(f)) next
    bed <- fread(f)
    lens <- chrom_lengths(P$mouse[[b]])
    plots <- c(plots, sprintf("<img alt='expected mouse CN %s' src='%s'/>", b, png_uri(ideogram(bed, lens, b), 9, 6)),
               html_table(bed[, .(interval = paste0(`#chrom`, ":", start + 1, "-", end), direction, states, n_genes, confidence,
                                 genes = substr(genes, 1, 120))], 60))
  }
  sec$cn <- paste0("<h2>Expected mouse CN map</h2><p>Merged 1:1-ortholog intervals; LOH drawn as a second track.</p>",
                   paste(plots, collapse = ""),
                   "<h3>Direct liftOver cross-check</h3>", html_table(read_layer(P, "cn_liftover_crosscheck.tsv")),
                   "<h3>Mapping confidence of human segments (", esc(hb), ")</h3>",
                   html_table({ s <- read_layer(P, "cn_segment_mapping.tsv"); if (is.null(s)) NULL else s[build == hb, .N, by = .(direction, confidence)][order(direction, confidence)] }),
                   "<h3>Human segments split across mouse chromosomes</h3>",
                   html_table({ s <- read_layer(P, "cn_segment_mapping.tsv"); if (is.null(s)) NULL else s[build == hb & split_across_chromosomes == TRUE & direction != "neutral",
                                                                                                   .(segment, state, n_genes, mouse_chroms, n_blocks, confidence)] }))

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
