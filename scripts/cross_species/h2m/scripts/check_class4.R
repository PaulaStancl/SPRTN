# ---------------------------------------------------------------------------
# check_class4.R - are H2M class 4 calls ("mutated sequences are not identical") real, or an
# artefact of H2M's alignment / transcript choice? Run after 01.
#
#   Rscript scripts/check_class4.R [--config ...] [--release ...] [--all-classes] [--h2m-plots] [--build GRCm38]
#
# H2M calls class 4 when the mutated residues / bases are not identical between human and mouse. For every
# class 4 variant this script works out why, from the transcripts themselves (release GTF exons + CDS, FASTAs):
#   1. the human change is applied to the human mRNA (coding strand, so minus-strand genes such as AXIN1 are
#      reverse-complemented): first and last affected codon from the coding phase, not from HGVSp;
#   2. human protein vs every protein-coding isoform of the mouse ortholog, global BLOSUM62 alignment (gap open
#      10, extension 0.5); the affected human codons are mapped to mouse residues and codons;
#   3. span verdict (H2M's mouse transcript first): span_identical (residues and codons identical - class 4 not
#      explained by the span), identical_in_other_isoform, aa_differs, nt_differs_synonymous (same residues,
#      different codons), gap_in_mouse; change_outside_exons / change_outside_cds for splice / UTR changes;
#   4. consequence: the human change translated to its new stop (premature_stop / extended / same_length /
#      no_stop_in_transcript; NMD = stop ends > 50 nt before the last exon-exon junction). For truncating
#      variants a hypothetical analogous mouse change (same deletion length / inserted bases at the aligned codon
#      and phase, or a stop codon for nonsense) gets the same consequence;
#   5. functional_equivalence: predicted (both premature stop with the same NMD outcome), unresolved, or
#      not_applicable (missense). A prediction, not a validation.
# --build GRCm38 repeats the protein check against the GRCm38 / Ensembl 102 isoforms (H2M itself models in
# GRCm39 only): class4_check.GRCm38.tsv / class4_alignments.GRCm38.html. Same verdicts = the class 4 calls hold
# for GRCm38. Transcript IDs missing from Ensembl 102 count as "no canonical"; the other isoforms are still used.
# --all-classes checks every modelled variant, not only class 4 (a sanity run of the method).
# --h2m-plots adds H2M's own result plot (h2m_visualize.py: model(show_sequence=True) + visualization(), as in
# the H2M tutorial) for the checked variants and for every modelled PURPLE-reported / driver-gene / hotspot
# variant. H2M can only draw variants it modelled (class 0-2); for class 4 the page says why not.
#
# Input   human_predictions/expected_mouse_snv.tsv, release GTFs (00), human_fasta, GRCm39 FASTA
# Output  human_predictions/class4_check.tsv   variant_id gene_symbol_h hgvsp_h pave_coding_effect h2m_class
#           tx_h human_strand first_codon last_codon human_aa protein_pos gene_id_m symbol_m tx_m n_mouse_isoforms
#           span_aa_h span_aa_m span_codons_h span_codons_m mouse_first_codon local_identity_canonical best_isoform
#           human_outcome human_wt_aa human_mut_aa human_nmd mouse_outcome mouse_wt_aa mouse_mut_aa mouse_nmd
#           verdict functional_equivalence
#         human_predictions/class4_alignments.html  per variant: consequences in both species, protein views
#           (H2M-style and BLOSUM62, +/-20 aa), coding-DNA codons human vs mouse (+/-6 codons, deleted bases struck
#           through), the isoform with an identical span if any, and H2M's own plot (--h2m-plots)
#         logs/check_class4.log (counts per verdict)
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

# analogous change in mouse: at mRNA index idx delete n_del bases and insert ins (coding strand); for a
# nonsense SNV/MNV (n_del = 0, stop = TRUE) the codon at idx becomes TAA
analog_mutation <- function(tm, idx, n_del, ins, stop = FALSE) {
  b <- strsplit(tm$mrna, "", fixed = TRUE)[[1]]
  if (stop) { b[idx:(idx + 2L)] <- c("T", "A", "A"); return(list(mrna = paste(b, collapse = ""), first = idx, delta = 0L)) }
  if (n_del > 0) b <- b[-(idx:(idx + n_del - 1L))]
  if (nzchar(ins)) b <- append(b, strsplit(ins, "")[[1]], after = idx - 1L)
  list(mrna = paste(b, collapse = ""), first = idx, delta = nchar(ins) - n_del)
}

# HTML: codons of human and mouse along the protein alignment, k columns either side of human residue hpos;
# deleted human bases (mRNA indices del_idx) struck through
codon_html <- function(al, hpos, tmh, tmm, del_idx, label, k = 6L) {
  a <- strsplit(al$p, "", fixed = TRUE)[[1]]; s <- strsplit(al$s, "", fixed = TRUE)[[1]]
  col <- which(cumsum(a != "-") == hpos & a != "-")[1]
  if (is.na(col)) return("")
  w <- max(1L, col - k):min(length(a), col + k)
  hi <- cumsum(a != "-"); mi <- cumsum(s != "-")
  cod <- function(tm, r) substr(tm$mrna, tm$cds_start + 3L * (r - 1L), tm$cds_start + 3L * (r - 1L) + 2L)
  hrow <- mrow <- prow_h <- prow_m <- character()
  for (j in w) {
    hc <- if (a[j] == "-") "---" else cod(tmh, hi[j]); mc <- if (s[j] == "-") "---" else cod(tmm, mi[j])
    hb <- strsplit(hc, "")[[1]]; mb <- strsplit(mc, "")[[1]]
    hidx <- if (a[j] == "-") rep(NA, 3) else tmh$cds_start + 3L * (hi[j] - 1L) + 0:2
    cls <- ifelse(hb == "-" | mb == "-", "gap", ifelse(hb == mb, "id", "mm"))
    del <- !is.na(hidx) & hidx %in% del_idx
    v <- if (j == col) " var" else ""
    hrow <- c(hrow, paste0("<span class='", cls, ifelse(del, " del", ""), v, "'>", hb, "</span>", collapse = ""), " ")
    mrow <- c(mrow, paste0("<span class='", cls, v, "'>", mb, "</span>", collapse = ""), " ")
    prow_h <- c(prow_h, sprintf(" %s  ", a[j])); prow_m <- c(prow_m, sprintf(" %s  ", s[j]))
  }
  sprintf("<div class='al'><b>%s</b><pre>human aa   %s\nhuman DNA  %s\nmouse DNA  %s\nmouse aa   %s</pre></div>",
          label, paste(prow_h, collapse = ""), paste(hrow, collapse = ""), paste(mrow, collapse = ""), paste(prow_m, collapse = ""))
}

# HTML block of an alignment window: k columns either side of the human position hpos
aln_html <- function(al, hpos, label, k = 20L) {
  a <- strsplit(al$p, "", fixed = TRUE)[[1]]; s <- strsplit(al$s, "", fixed = TRUE)[[1]]
  col <- which(cumsum(a != "-") == hpos & a != "-")[1]
  if (is.na(col)) return(sprintf("<div class='al'><b>%s</b>: position not found</div>", label))
  w <- max(1L, col - k):min(length(a), col + k)
  hstart <- sum(a[seq_len(w[1])] != "-") + (a[w[1]] == "-")
  mstart <- sum(s[seq_len(w[1])] != "-") + (s[w[1]] == "-")
  cls <- ifelse(a[w] == "-" | s[w] == "-", "gap", ifelse(a[w] == s[w], "id", "mm"))
  cell <- function(ch) paste0("<span class='", cls, ifelse(w == col, " var", ""), "'>", ch, "</span>", collapse = "")
  mid <- paste0("<span class='", cls, ifelse(w == col, " var", ""), "'>", ifelse(cls == "id", "|", ifelse(cls == "mm", ".", " ")), "</span>", collapse = "")
  sprintf("<div class='al'><b>%s</b><pre>human %5d %s\n            %s\nmouse %5d %s</pre></div>",
          label, hstart, cell(a[w]), mid, mstart, cell(s[w]))
}

main <- function() {
  a <- parse_cli(c("all-classes", "h2m-plots"))
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("check_class4", P)
  build <- if (is.null(a$build)) cfg$h2m_mouse_build else a$build   # --build GRCm38: same check on Ensembl 102 isoforms
  if (!build %chin% names(P$mouse)) stop("unknown --build ", build, call. = FALSE)
  M <- P$mouse[[build]]
  sfx <- if (build == cfg$h2m_mouse_build) "" else paste0(".", build)
  log("mouse build %s (H2M modelled in %s)", build, cfg$h2m_mouse_build)

  snv <- read_layer(P, "expected_mouse_snv.tsv")
  if (is.null(snv)) stop("expected_mouse_snv.tsv missing - run 01 first", call. = FALSE)
  x <- if (isTRUE(a$all_classes)) snv[!is.na(h2m_class) & !is.na(gene_id_m)] else snv[h2m_class %in% 4L]
  log("%d variants to check (%s)", nrow(x), if (isTRUE(a$all_classes)) "all modelled" else "H2M class 4")
  if (!nrow(x)) { write_tsv(data.table(), file.path(P$hp, paste0("class4_check", sfx, ".tsv")), log); return(invisible()) }

  pp <- parse_hgvsp(x$hgvsp_h)
  x[, `:=`(human_aa = pp$aa, protein_pos = pp$pos)]

  mtx <- load_genes(M$tx)
  iso <- mtx[gene_id %chin% x$gene_id_m & has_cds == TRUE, .(gene_id, transcript_id, canonical)]
  log("reading exons + CDS: %d human transcripts, %d mouse isoforms of %d genes", uniqueN(x$tx_h), nrow(iso), uniqueN(iso$gene_id))
  hfeat <- attach_seq(read_tx_features(P$human_gtf, tx_ids = unique(na.omit(x$tx_h))), cfg$human_fasta)
  mfeat <- attach_seq(read_tx_features(M$gtf, gene_ids = unique(iso$gene_id))[transcript_id %chin% iso$transcript_id], M$fasta)
  hprot <- proteins_from(hfeat); mprot <- proteins_from(mfeat)
  hmod <- lapply(split(hfeat, by = "transcript_id"), tx_model)
  mmod <- lapply(split(mfeat, by = "transcript_id"), tx_model)

  res <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
    r <- x[i]
    out <- data.table(variant_id = r$variant_id, gene_symbol_h = r$gene_symbol_h, hgvsp_h = r$hgvsp_h,
                      pave_coding_effect = r$pave_coding_effect, h2m_class = r$h2m_class, tx_h = r$tx_h, human_strand = NA_character_,
                      first_codon = NA_integer_, last_codon = NA_integer_, human_aa = r$human_aa, protein_pos = r$protein_pos,
                      gene_id_m = r$gene_id_m, symbol_m = r$symbol_m, tx_m = r$tx_m, n_mouse_isoforms = 0L,
                      span_aa_h = NA_character_, span_aa_m = NA_character_, span_codons_h = NA_character_, span_codons_m = NA_character_,
                      mouse_first_codon = NA_integer_, local_identity_canonical = NA_real_, best_isoform = NA_character_,
                      human_outcome = NA_character_, human_wt_aa = NA_integer_, human_mut_aa = NA_integer_, human_nmd = NA,
                      mouse_outcome = NA_character_, mouse_wt_aa = NA_integer_, mouse_mut_aa = NA_integer_, mouse_nmd = NA,
                      verdict = NA_character_, functional_equivalence = NA_character_, html = NA_character_)
    tmh <- hmod[[r$tx_h]]; hp <- hprot[r$tx_h]
    if (is.null(tmh) || is.na(hp)) return(out[, verdict := "human_transcript_not_built"])
    out[, human_strand := tmh$strand]
    mut <- apply_vcf_to_tx(tmh, as.integer(r$pos), r$ref, r$alt)
    if (is.null(mut)) return(out[, verdict := "change_outside_exons"])
    c0 <- mut$first - tmh$cds_start
    if (c0 < 0 || mut$first > tmh$cds_end + 3L) return(out[, verdict := "change_outside_cds"])
    r1 <- c0 %/% 3L + 1L; phase <- c0 %% 3L
    r2 <- if (mut$n_del > 0) (mut$first + mut$n_del - 1L - tmh$cds_start) %/% 3L + 1L
          else if (nchar(r$ref) == nchar(r$alt)) (mut$first + nchar(r$ref) - 1L - tmh$cds_start) %/% 3L + 1L else r1
    hc <- mrna_consequence(tmh, mut$mrna, mut$first, mut$delta)
    out[, `:=`(first_codon = r1, last_codon = r2, span_aa_h = substr(hp, r1, r2), human_outcome = hc$outcome,
               human_wt_aa = hc$wt_aa, human_mut_aa = hc$mut_aa, human_nmd = hc$nmd)]
    tx <- iso[gene_id == r$gene_id_m, transcript_id]
    tx <- intersect(c(r$tx_m, setdiff(tx, r$tx_m)), intersect(names(mprot), names(mmod)))   # H2M's transcript first
    out[, n_mouse_isoforms := length(tx)]
    if (!length(tx)) return(out[, verdict := "mouse_transcript_not_built"])
    # span check per isoform: aligned mouse residues of human codons r1..r2, then their codons
    span_of <- function(t) {
      al <- align_proteins(hp, mprot[[t]]); tm <- mmod[[t]]
      mp <- lapply(r1:r2, function(rr) map_aligned_pos(al$p, al$s, rr))
      maa <- vapply(mp, function(z) fcoalesce(z$mouse_aa, "-"), "")
      mpos <- vapply(mp, function(z) fcoalesce(z$mouse_pos, NA_integer_), 0L)
      cod <- function(tm, rr) if (is.na(rr)) "---" else substr(tm$mrna, tm$cds_start + 3L * (rr - 1L), tm$cds_start + 3L * (rr - 1L) + 2L)
      hcod <- vapply(r1:r2, function(rr) cod(tmh, rr), ""); mcod <- vapply(mpos, function(rr) cod(tm, rr), "")
      v <- if (any(maa == "-")) "gap_in_mouse" else if (any(maa != strsplit(substr(hp, r1, r2), "")[[1]])) "aa_differs"
           else if (any(hcod != mcod)) "nt_differs_synonymous" else "span_identical"
      list(tx = t, al = al, tm = tm, verdict = v, maa = paste(maa, collapse = ""), hcod = paste(hcod, collapse = " "),
           mcod = paste(mcod, collapse = " "), m1 = mpos[1], ident = mp[[1]]$local_identity)
    }
    sp <- lapply(setNames(tx, tx), span_of)
    can <- sp[[tx[1]]]
    other <- Filter(function(z) z$verdict == "span_identical", sp[-1])
    out[, `:=`(span_aa_m = can$maa, span_codons_h = can$hcod, span_codons_m = can$mcod, mouse_first_codon = can$m1,
               local_identity_canonical = round(can$ident, 2), best_isoform = if (length(other)) other[[1]]$tx else tx[1],
               verdict = if (can$verdict != "span_identical" && length(other)) "identical_in_other_isoform" else can$verdict)]
    # analogous mouse change at the aligned codon (same within-codon phase): same deletion length / inserted bases,
    # or a stop codon for a nonsense variant; missense has no "analogous" truncation
    trunc <- r$pave_coding_effect %chin% c("NONSENSE_OR_FRAMESHIFT")
    if (trunc && !is.na(can$m1)) {
      idx <- can$tm$cds_start + 3L * (can$m1 - 1L) + phase
      nonsense <- nchar(r$ref) == nchar(r$alt)
      am <- analog_mutation(can$tm, if (nonsense) can$tm$cds_start + 3L * (can$m1 - 1L) else idx, mut$n_del, mut$ins, stop = nonsense)
      mc <- mrna_consequence(can$tm, am$mrna, am$first, am$delta)
      out[, `:=`(mouse_outcome = mc$outcome, mouse_wt_aa = mc$wt_aa, mouse_mut_aa = mc$mut_aa, mouse_nmd = mc$nmd)]
      same_class <- hc$outcome == "premature_stop" && mc$outcome == "premature_stop" && identical(hc$nmd, mc$nmd)
      out[, functional_equivalence := if (same_class) "predicted" else "unresolved"]
    } else out[, functional_equivalence := if (trunc) "unresolved" else "not_applicable"]
    # HTML: protein (H2M-style + proper), DNA codons, other isoform
    del_idx <- if (mut$n_del > 0) mut$first:(mut$first + mut$n_del - 1L) else integer()
    blocks <- c(aln_html(align_proteins(hp, mprot[[tx[1]]], "h2m"), r1, sprintf("Protein, H2M-style alignment (match only, free gaps) vs %s", tx[1])),
                aln_html(can$al, r1, sprintf("Protein, BLOSUM62 + gap penalties vs %s", tx[1])),
                codon_html(can$al, r1, tmh, can$tm, del_idx, sprintf("Coding DNA (coding strand; human gene on %s strand), deleted bases struck through", tmh$strand)))
    if (length(other)) blocks <- c(blocks, aln_html(other[[1]]$al, r1, sprintf("Protein, BLOSUM62 vs other mouse isoform %s", other[[1]]$tx)))
    cons <- sprintf("<p>Human %s -> %s, %s aa of %s%s. Analogous mouse change -> %s%s.</p>",
                    r$hgvsp_h, hc$outcome, hc$mut_aa, hc$wt_aa, if (isTRUE(hc$nmd)) ", NMD expected" else if (isFALSE(hc$nmd)) ", no NMD" else "",
                    fcoalesce(out$mouse_outcome, "not computed"),
                    if (!is.na(out$mouse_mut_aa)) sprintf(", %s aa of %s%s", out$mouse_mut_aa, out$mouse_wt_aa,
                                                          if (isTRUE(out$mouse_nmd)) ", NMD expected" else if (isFALSE(out$mouse_nmd)) ", no NMD" else "") else "")
    out[, html := paste(c(cons, blocks), collapse = "")]
    out
  }))
  if (!"html" %chin% names(res)) res[, html := NA_character_]
  f <- write_tsv(res[, !"html"], file.path(P$hp, paste0("class4_check", sfx, ".tsv")), log)
  # ---- H2M's own plot (tutorial: model(show_sequence=True) + visualization) ------------------------------
  # for the checked variants (class 4: H2M cannot draw, the message says so) and the modelled driver variants
  res[, h2m_plot := NA_character_]
  drv_html <- character()
  if (isTRUE(a$h2m_plots) && build != cfg$h2m_mouse_build) log("--h2m-plots only for %s (H2M's build) - skipped", cfg$h2m_mouse_build, level = "WARNING")
  if (isTRUE(a$h2m_plots) && build == cfg$h2m_mouse_build) {
    lab <- function(d) gsub("[^A-Za-z0-9]+", "_", paste(d$variant_id, d$gene_id_m))
    drv <- snv[h2m_class %in% 0:2 & !is.na(gene_id_m) & (truthy(reported) | truthy(purple_driver_gene) | truthy(hotspot))]
    sel <- unique(rbind(x[, .(variant_id, gene_id_m, gene_symbol_h, hgvsp_h, pos, ref, alt, tx_h, tx_m, set = "checked")],
                        drv[, .(variant_id, gene_id_m, gene_symbol_h, hgvsp_h, pos, ref, alt, tx_h, tx_m, set = "driver")]),
                  by = c("variant_id", "gene_id_m"))[!is.na(tx_h) & !is.na(tx_m)]
    sel[, label := lab(sel)]
    vin <- cbind(sel[, .(label, tx_id_h = tx_h, tx_id_m = tx_m)], vcf_to_maf(as.integer(sel$pos), sel$ref, sel$alt))
    pdir <- file.path(P$hp, "h2m_plots"); fin <- file.path(P$tmp, "h2m_viz_input.tsv"); fst <- file.path(pdir, "status.tsv")
    write_tsv(vin, fin, log)
    py <- if (is.null(cfg$python)) "python" else cfg$python
    args <- c(file.path(H2M_SCRIPTS, "h2m_visualize.py"), "--input", fin, "--outdir", pdir, "--status", fst,
              "--human-fasta", cfg$human_fasta, "--mouse-fasta", M$fasta, "--human-gtf", P$human_gtf, "--mouse-gtf", M$gtf,
              "--workdir", file.path(P$tmp, "h2m_viz"))
    log("H2M plots for %d variants", nrow(vin))
    if (system2(py, shQuote(args)) != 0) log("h2m_visualize.py failed - no H2M plots", level = "WARNING")
    st <- if (file.exists(fst)) fread(fst) else data.table(label = character(), plotted = logical(), message = character())
    plot_html <- function(l) {
      k <- match(l, st$label)
      if (is.na(k)) return("<p class='v'>no H2M plot</p>")
      png <- file.path(pdir, paste0(l, ".png"))
      if (isTRUE(st$plotted[k]) && file.exists(png))
        sprintf("<div class='al'><b>H2M visualization</b><br><img style='max-width:100%%' src='%s'/></div>", base64enc::dataURI(file = png, mime = "image/png"))
      else sprintf("<p class='v'>H2M visualization: %s</p>", st$message[k])
    }
    res[, h2m_plot := vapply(lab(res), plot_html, "")]
    d <- sel[set == "driver"]
    drv_html <- vapply(seq_len(nrow(d)), function(i) sprintf("<h2>%s %s</h2>%s", d$gene_symbol_h[i], fcoalesce(d$hgvsp_h[i], ""), plot_html(d$label[i])), "")
  }

  fh <- file.path(P$hp, paste0("class4_alignments", sfx, ".html"))
  css <- "body{font:14px system-ui,sans-serif;max-width:1100px;margin:auto;padding:16px}pre{font:13px/1.35 ui-monospace,Menlo,monospace;overflow-x:auto}
.id{background:#d4efdf}.mm{background:#fad7a0}.gap{color:#999;background:#eee}.var{outline:2px solid #c0392b;font-weight:bold}
.v{display:inline-block;padding:1px 8px;border-radius:9px;background:#eee;font-size:12px}h2{font-size:16px;margin-top:28px;border-top:1px solid #ccc;padding-top:10px}
.al{margin:6px 0 12px}.del{text-decoration:line-through;color:#c0392b}"
  body <- vapply(seq_len(nrow(res)), function(i) with(res[i], sprintf(
    "<h2>%s %s <span class='v'>%s</span> <span class='v'>functional equivalence: %s</span></h2><p>%s | %s | H2M class %s | %s (%s strand) codons %s-%s: %s -> mouse %s (%s), codons %s -> %s, local identity %s</p>%s",
    gene_symbol_h, fcoalesce(hgvsp_h, ""), verdict, fcoalesce(functional_equivalence, ""), variant_id, fcoalesce(pave_coding_effect, ""),
    h2m_class, tx_h, fcoalesce(human_strand, "?"), fcoalesce(as.character(first_codon), "?"), fcoalesce(as.character(last_codon), "?"),
    fcoalesce(span_aa_h, "?"), fcoalesce(span_aa_m, "?"), fcoalesce(tx_m, "?"), fcoalesce(span_codons_h, "?"), fcoalesce(span_codons_m, "?"),
    fcoalesce(as.character(local_identity_canonical), "NA"), paste0(fcoalesce(html, ""), fcoalesce(h2m_plot, "")))), "")
  writeLines(c("<!doctype html><html><head><meta charset='utf-8'><title>H2M class 4 alignments</title><style>", css, "</style></head><body>",
               sprintf("<h1>H2M class 4 - alignment check (mouse %s)</h1>", build),
               "<p>Green = identical, orange = different, grey = gap, red box = first affected residue / codon, struck-through red = deleted human bases. ",
               "First block: H2M's own scoring (globalxx). Second: a proper protein alignment (BLOSUM62, gap open 10, extend 0.5). ",
               "If the residue is identical in the second block but not the first, class 4 comes from H2M's alignment; if identical only ",
               "in another isoform, from the transcript choice; if different in all, the residue really differs in mouse.</p>",
               body,
               if (length(drv_html)) c("<h1>H2M visualization of the modelled driver / reported variants</h1>",
                                       "<p>H2M's own plot (tutorial: show_sequence=True, flank 4, print size 2): human WT / mutant, mouse WT, original and alternative model.</p>", drv_html),
               "</body></html>"), fh)
  log("wrote %s", fh)
  for (v in res[, .N, by = verdict][order(-N), paste0(verdict, ": ", N)]) log("  %s", v)
  write_provenance(P, cfg, "check_class4", resources = list(human_gtf = P$human_gtf, mouse_gtf = M$gtf), outputs = f,
                   extra = list(alignment = "global, BLOSUM62, gap opening 10, extension 0.5",
                                aligner = if (identical(ALN, asNamespace("pwalign"))) "pwalign" else "Biostrings"))
  log("done")
}

main()
