# ---------------------------------------------------------------------------
# check_class4.R - are H2M class 4 calls ("mutated sequences are not identical") real, or an
# artefact of H2M's alignment / transcript choice? Run after 01.
#
#   Rscript scripts/check_class4.R [--config ...] [--release ...] [--all-classes]
#
# H2M aligns the human and mouse canonical proteins with pairwise2.globalxx (match +1, no mismatch or
# gap penalty) and calls class 4 when the mutated residue(s) are not identical in that alignment. Here,
# for every class 4 variant with a protein position (PAVE HGVSp):
#   1. the human protein of H2M's human transcript (tx_h) and EVERY protein-coding isoform of the mouse
#      ortholog are built from the release GTFs (CDS) and FASTAs (human_fasta, GRCm39);
#   2. the human residue at the HGVSp position is checked against PAVE (hgvsp_mismatch if not: tx_h is
#      not the transcript PAVE used);
#   3. human vs each mouse isoform is realigned globally with BLOSUM62, gap opening 10, extension 0.5;
#   4. verdict: artefact_identical_in_canonical (same residue in H2M's mouse transcript - class 4 is an
#      alignment artefact), identical_in_other_isoform (same residue, but only in another mouse isoform -
#      artefact of transcript choice), gap_in_mouse (no mouse residue at that position in any isoform),
#      confirmed_differs (a different residue in every isoform), not_protein_level (no HGVSp: splice /
#      non-coding), hgvsp_mismatch.
# --all-classes checks every modelled variant, not only class 4 (a sanity run of the method).
#
# Input   human_predictions/expected_mouse_snv.tsv, release GTFs (00), human_fasta, GRCm39 FASTA
# Output  human_predictions/class4_check.tsv   variant_id gene_symbol_h hgvsp_h pave_coding_effect
#           h2m_class human_aa protein_pos human_aa_check symbol_m tx_m mouse_pos_canonical
#           mouse_aa_canonical local_identity_canonical context_h context_m_canonical
#           n_mouse_isoforms best_isoform mouse_aa_best verdict
#         human_predictions/class4_alignments.html  per variant: +/-20 aa of the H2M-style alignment and the
#           BLOSUM62 alignment vs H2M's mouse transcript (+ the isoform with an identical residue), coloured
#           identical / different / gap, the mutated position boxed
#         logs/check_class4.log (counts per verdict)
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

ALN <- if (requireNamespace("pwalign", quietly = TRUE)) asNamespace("pwalign") else asNamespace("Biostrings")

# CDS rows (chrom start end strand frame transcript_id gene_id) of the given transcripts / genes
read_cds <- function(gtf, tx_ids = character(), gene_ids = character()) {
  cmd <- sprintf("gzip -dc %s | awk -F'\\t' '$3==\"CDS\"'", shQuote(gtf))
  g <- fread(cmd = cmd, header = FALSE, sep = "\t", quote = "", select = c(1, 4, 5, 7, 8, 9),
             col.names = c("chrom", "start", "end", "strand", "frame", "attr"), colClasses = list(character = c(1, 8, 9)))
  g[, `:=`(transcript_id = gtf_attr(attr, "transcript_id"), gene_id = gtf_attr(attr, "gene_id"))]
  g[transcript_id %chin% tx_ids | gene_id %chin% gene_ids, !"attr"]
}

# transcript -> protein string, reading the CDS pieces from an indexed FASTA
proteins_from <- function(cds, fasta) {
  fa <- Rsamtools::FaFile(fasta)
  refs <- GenomeInfoDb::seqnames(Rsamtools::seqinfo(fa))
  cds[, fa_chrom := fifelse(chrom %chin% refs, chrom, fifelse(paste0("chr", chrom) %chin% refs, paste0("chr", chrom), sub("^chr", "", chrom)))]
  cds <- cds[fa_chrom %chin% refs]
  if (!nrow(cds)) return(setNames(character(), character()))
  cds[, seq := as.character(Rsamtools::scanFa(fa, GenomicRanges::GRanges(fa_chrom, IRanges::IRanges(start, end))))]
  p <- cds[, .(protein = cds_protein(.SD, seq)), by = transcript_id]
  setNames(sub("\\*$", "", p$protein), p$transcript_id)
}

# proper = BLOSUM62 with gap penalties; h2m = H2M's scoring (pairwise2.globalxx: match 1, mismatch 0, free gaps).
# Ties can be broken differently than in H2M, so the h2m view reproduces its scoring, not necessarily its exact path.
align <- function(h, m, mode = c("proper", "h2m")) {
  mode <- match.arg(mode)
  if (mode == "proper") {
    a <- ALN$pairwiseAlignment(Biostrings::AAString(h), Biostrings::AAString(m), substitutionMatrix = "BLOSUM62",
                               gapOpening = 10, gapExtension = 0.5, type = "global")
  } else {
    l <- sort(unique(c(strsplit(h, "")[[1]], strsplit(m, "")[[1]])))
    idm <- diag(1, length(l)); dimnames(idm) <- list(l, l)
    a <- ALN$pairwiseAlignment(Biostrings::AAString(h), Biostrings::AAString(m), substitutionMatrix = idm,
                               gapOpening = 0, gapExtension = 0, type = "global")
  }
  list(p = as.character(ALN$alignedPattern(a)), s = as.character(ALN$alignedSubject(a)))
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

ctx <- function(p, pos, k = 5L) if (is.na(pos)) NA_character_ else
  paste0(substr(p, max(1L, pos - k), pos - 1L), "[", substr(p, pos, pos), "]", substr(p, pos + 1L, pos + k))

main <- function() {
  a <- parse_cli("all-classes")
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("check_class4", P)
  M <- P$mouse[[cfg$h2m_mouse_build]]

  snv <- read_layer(P, "expected_mouse_snv.tsv")
  if (is.null(snv)) stop("expected_mouse_snv.tsv missing - run 01 first", call. = FALSE)
  x <- if (isTRUE(a$all_classes)) snv[!is.na(h2m_class) & !is.na(gene_id_m)] else snv[h2m_class %in% 4L]
  log("%d variants to check (%s)", nrow(x), if (isTRUE(a$all_classes)) "all modelled" else "H2M class 4")
  if (!nrow(x)) { write_tsv(data.table(), file.path(P$hp, "class4_check.tsv"), log); return(invisible()) }

  pp <- parse_hgvsp(x$hgvsp_h)
  x[, `:=`(human_aa = pp$aa, protein_pos = pp$pos)]

  mtx <- load_genes(M$tx)
  iso <- mtx[gene_id %chin% x$gene_id_m & has_cds == TRUE, .(gene_id, transcript_id, canonical)]
  log("reading CDS: %d human transcripts, %d mouse isoforms of %d genes", uniqueN(x$tx_h), nrow(iso), uniqueN(iso$gene_id))
  hprot <- proteins_from(read_cds(P$human_gtf, tx_ids = unique(na.omit(x$tx_h))), cfg$human_fasta)
  mprot <- proteins_from(read_cds(M$gtf, gene_ids = unique(iso$gene_id)), M$fasta)

  res <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
    r <- x[i]
    out <- data.table(variant_id = r$variant_id, gene_symbol_h = r$gene_symbol_h, hgvsp_h = r$hgvsp_h,
                      pave_coding_effect = r$pave_coding_effect, h2m_class = r$h2m_class, human_aa = r$human_aa,
                      protein_pos = r$protein_pos, human_aa_check = NA, symbol_m = r$symbol_m, tx_m = r$tx_m,
                      mouse_pos_canonical = NA_integer_, mouse_aa_canonical = NA_character_, local_identity_canonical = NA_real_,
                      context_h = NA_character_, context_m_canonical = NA_character_, n_mouse_isoforms = 0L,
                      best_isoform = NA_character_, mouse_aa_best = NA_character_, verdict = NA_character_)
    if (is.na(r$protein_pos)) return(out[, verdict := "not_protein_level"])
    hp <- hprot[r$tx_h]
    if (is.na(hp)) return(out[, verdict := "human_protein_not_built"])
    out[, `:=`(human_aa_check = substr(hp, r$protein_pos, r$protein_pos) == r$human_aa, context_h = ctx(hp, r$protein_pos))]
    if (!isTRUE(out$human_aa_check)) return(out[, verdict := "hgvsp_mismatch"])
    tx <- iso[gene_id == r$gene_id_m, transcript_id]
    tx <- intersect(c(r$tx_m, setdiff(tx, r$tx_m)), names(mprot))      # canonical (H2M's) first
    out[, n_mouse_isoforms := length(tx)]
    if (!length(tx)) return(out[, verdict := "mouse_protein_not_built"])
    als <- lapply(setNames(tx, tx), function(t) align(hp, mprot[[t]]))
    hits <- rbindlist(lapply(tx, function(t) {
      mp <- map_aligned_pos(als[[t]]$p, als[[t]]$s, r$protein_pos)
      data.table(tx = t, mouse_pos = mp$mouse_pos, mouse_aa = mp$mouse_aa, local_identity = mp$local_identity,
                 context_m = ctx(mprot[[t]], mp$mouse_pos))
    }))
    can <- hits[tx == r$tx_m]
    if (nrow(can)) out[, `:=`(mouse_pos_canonical = can$mouse_pos, mouse_aa_canonical = can$mouse_aa,
                              local_identity_canonical = round(can$local_identity, 2), context_m_canonical = can$context_m)]
    same <- hits[mouse_aa == r$human_aa]
    best <- if (nrow(same)) same[1] else hits[order(-local_identity)][1]
    out[, `:=`(best_isoform = best$tx, mouse_aa_best = best$mouse_aa)]
    out[, verdict := fcase(nrow(can) > 0 && identical(can$mouse_aa, r$human_aa), "artefact_identical_in_canonical",
                           nrow(same) > 0, "identical_in_other_isoform",
                           all(hits$mouse_aa == "-" | is.na(hits$mouse_aa)), "gap_in_mouse",
                           default = "confirmed_differs")]
    # alignment views for the HTML: H2M-style and proper vs H2M's mouse transcript, + the isoform that matches
    blocks <- character()
    if (r$tx_m %chin% tx) blocks <- c(aln_html(align(hp, mprot[[r$tx_m]], "h2m"), r$protein_pos, sprintf("H2M-style alignment (match only, free gaps) vs %s (H2M's mouse transcript)", r$tx_m)),
                                      aln_html(als[[r$tx_m]], r$protein_pos, sprintf("BLOSUM62 + gap penalties vs %s (H2M's mouse transcript)", r$tx_m)))
    if (!is.na(best$tx) && best$tx != r$tx_m) blocks <- c(blocks, aln_html(als[[best$tx]], r$protein_pos, sprintf("BLOSUM62 vs other mouse isoform %s", best$tx)))
    out[, html := paste(blocks, collapse = "")]
    out
  }))
  if (!"html" %chin% names(res)) res[, html := NA_character_]
  f <- write_tsv(res[, !"html"], file.path(P$hp, "class4_check.tsv"), log)
  fh <- file.path(P$hp, "class4_alignments.html")
  css <- "body{font:14px system-ui,sans-serif;max-width:1100px;margin:auto;padding:16px}pre{font:13px/1.35 ui-monospace,Menlo,monospace;overflow-x:auto}
.id{background:#d4efdf}.mm{background:#fad7a0}.gap{color:#999;background:#eee}.var{outline:2px solid #c0392b;font-weight:bold}
.v{display:inline-block;padding:1px 8px;border-radius:9px;background:#eee;font-size:12px}h2{font-size:16px;margin-top:28px;border-top:1px solid #ccc;padding-top:10px}
.al{margin:6px 0 12px}"
  body <- vapply(seq_len(nrow(res)), function(i) with(res[i], sprintf(
    "<h2>%s %s <span class='v'>%s</span></h2><p>%s | %s | H2M class %s | human %s%s -> mouse (H2M transcript) %s%s, local identity %s</p>%s",
    gene_symbol_h, fcoalesce(hgvsp_h, ""), verdict, variant_id, fcoalesce(pave_coding_effect, ""), h2m_class,
    fcoalesce(human_aa, "?"), fcoalesce(as.character(protein_pos), "?"), fcoalesce(mouse_aa_canonical, "?"),
    fcoalesce(as.character(mouse_pos_canonical), ""), fcoalesce(as.character(local_identity_canonical), "NA"), fcoalesce(html, ""))), "")
  writeLines(c("<!doctype html><html><head><meta charset='utf-8'><title>H2M class 4 alignments</title><style>", css, "</style></head><body>",
               "<h1>H2M class 4 - alignment check</h1>",
               "<p>Green = identical, orange = different amino acid, grey = gap, red box = the mutated position. ",
               "First block: H2M's own scoring (globalxx). Second: a proper protein alignment (BLOSUM62, gap open 10, extend 0.5). ",
               "If the residue is identical in the second block but not the first, class 4 comes from H2M's alignment; if identical only ",
               "in another isoform, from the transcript choice; if different in all, the residue really differs in mouse.</p>",
               body, "</body></html>"), fh)
  log("wrote %s", fh)
  for (v in res[, .N, by = verdict][order(-N), paste0(verdict, ": ", N)]) log("  %s", v)
  write_provenance(P, cfg, "check_class4", resources = list(human_gtf = P$human_gtf, mouse_gtf = M$gtf), outputs = f,
                   extra = list(alignment = "global, BLOSUM62, gap opening 10, extension 0.5",
                                aligner = if (identical(ALN, asNamespace("pwalign"))) "pwalign" else "Biostrings"))
  log("done")
}

main()
