# ---------------------------------------------------------------------------
# 03_cn_to_mouse.R - PURPLE gene copy number -> expected mouse CN intervals per build.
#
#   Rscript scripts/03_cn_to_mouse.R [--config ...] [--release ...]
#
# Per gene (purple.cnv.gene.tsv, canonical rows): ploidy-relative state from config cn.*
# (AMP, GAIN, NEUTRAL, LOSS, DEEP_DEL; loss side from minCopyNumber, gain side from maxCopyNumber,
# partial = the two disagree) and LOH (minMinorAlleleCopyNumber < loh_minor_max, not with DEEP_DEL).
# Genes go to mouse through Ensembl 1:1 orthologs; 1:many / many:many go to a separate table.
# In each build, consecutive mouse genes with the same direction (gain / loss; LOH as its own
# track) merge into intervals; a run breaks at a gene of another direction, a gap >
# merge_max_gap_bp, a new chromosome or (split_on_block) a new synteny block. Synteny blocks:
# runs of 1:1 orthologs in mouse order whose human rank moves by <= synteny.max_rank_jump on one
# human chromosome (genome-wide, from the ortholog table; nothing downloaded).
#
# Mapping coverage per human segment (purple.cnv.somatic.tsv) and build: fraction of its genes with
# a 1:1 ortholog placed in the build, and fraction of its length that lifts over (segment cut into
# crosscheck_bin_bp bins, liftOver -minMatch=liftover.min_match). Confidence:
#   none      no gene mapped and < cn.conf_min_len_frac of the length lifts
#   ambiguous more of its genes are 1:many / many:many than 1:1
#   high      gene fraction >= conf_high_gene_frac (or no genes) and length fraction >= conf_high_len_frac
#   partial   anything else
# An interval takes the lowest confidence of its source segments.
# Cross-check: bins of non-neutral (and LOH) segments are lifted directly, merged per direction and
# compared (bp overlap, Jaccard) with the ortholog-based intervals.
#
# Input   <oncoanalyser>/purple/*.purple.cnv.gene.tsv, *.purple.cnv.somatic.tsv, *.purple.purity.tsv
# Output  human_predictions/
#   expected_mouse_cn.<build>.bed   #chrom start end name direction states source_human_segments genes
#                                   human_genes block_ids n_genes confidence   (0-based, half-open)
#   expected_mouse_cn_genes.tsv     every PURPLE gene: CN, state, direction, partial, loh, segment,
#                                   homology, mouse gene, per-build coords / block / in_build,
#                                   prediction (state, or no_ortholog / not_in_release / one2many /
#                                   not_in_build)
#   cn_one2many.tsv                 non-neutral or LOH genes with 1:many / many:many orthologs
#   cn_segment_mapping.tsv          per human segment x build: state, genes, mapped fractions, mouse
#                                   chromosomes, blocks, split_across_chromosomes, confidence
#   cn_liftover_crosscheck.tsv      per build x direction: bp, overlap, Jaccard
#   summary/cn.tsv
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

CONF_LEVELS <- c("high", "partial", "ambiguous", "none")

segment_confidence <- function(gene_frac, len_frac, n_one2one, n_multi, t) {
  fcase((is.na(gene_frac) | gene_frac == 0) & len_frac < t$conf_min_len_frac, "none",
        n_multi > n_one2one, "ambiguous",
        (is.na(gene_frac) | gene_frac >= t$conf_high_gene_frac) & len_frac >= t$conf_high_len_frac, "high",
        default = "partial")
}

worst_conf <- function(x) CONF_LEVELS[max(match(x, CONF_LEVELS, nomatch = 4L))]

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("03_cn_to_mouse", P)
  t <- cfg$cn
  pdir <- file.path(cfg$oncoanalyser_dir, "purple")
  f_gene <- find_one(pdir, "\\.purple\\.cnv\\.gene\\.tsv$", "PURPLE gene CN", log)
  f_seg <- find_one(pdir, "\\.purple\\.cnv\\.somatic\\.tsv$", "PURPLE segments", log)
  f_pur <- find_one(pdir, "\\.purple\\.purity\\.tsv$", "PURPLE purity", log)

  pur <- fread(f_pur)
  need_cols(pur, "ploidy", "purple.purity.tsv")
  ploidy <- as.numeric(pur$ploidy[1])
  log("ploidy %.3f; thresholds deep_del<%s, loss<%s x, gain>=%s x, amp>=%s x, LOH minor<%s", ploidy,
      t$deep_del_abs, t$loss_rel, t$gain_rel, t$amp_rel, t$loh_minor_max)

  g <- fread(f_gene, colClasses = list(character = "chromosome"))
  need_cols(g, c("chromosome", "start", "end", "gene", "minCopyNumber", "maxCopyNumber"), "purple.cnv.gene.tsv")
  if ("isCanonical" %chin% names(g)) g <- g[truthy(isCanonical)]
  g[, chromosome := norm_chrom(chromosome)]
  g <- g[chromosome %chin% STD_CHROMS][!duplicated(gene)]
  if (!"minMinorAlleleCopyNumber" %chin% names(g)) {
    log("no minMinorAlleleCopyNumber - LOH not assessable", level = "WARNING")
    g[, minMinorAlleleCopyNumber := NA_real_]
  }
  g[, c("state", "direction", "partial", "loh") := gene_cn_state(minCopyNumber, maxCopyNumber, minMinorAlleleCopyNumber, ploidy, t)]

  segs <- fread(f_seg, colClasses = list(character = "chromosome"))
  need_cols(segs, c("chromosome", "start", "end", "copyNumber"), "purple.cnv.somatic.tsv")
  segs[, chromosome := norm_chrom(chromosome)]
  segs <- segs[chromosome %chin% STD_CHROMS]
  segs[, seg_id := sprintf("%s:%.0f-%.0f", chromosome, start, end)]
  segs[, state := cn_state(copyNumber, ploidy, t)]
  segs[, direction := fcase(state %chin% GAIN_SIDE, "gain", state %chin% LOSS_SIDE, "loss", default = "neutral")]
  if (!"minorAlleleCopyNumber" %chin% names(segs)) segs[, minorAlleleCopyNumber := NA_real_]
  segs[, loh := !is.na(minorAlleleCopyNumber) & minorAlleleCopyNumber < t$loh_minor_max & state != "DEEP_DEL"]

  # gene midpoint -> segment
  g[, `:=`(mid = (start + end) %/% 2, mid2 = (start + end) %/% 2)]
  sk <- segs[, .(chromosome, s = start, e = end, seg_id)]
  setkey(sk, chromosome, s, e)
  ov <- foverlaps(g[, .(gene, chromosome, mid, mid2)], sk, by.x = c("chromosome", "mid", "mid2"), nomatch = NULL)
  g[ov[!duplicated(gene)], on = "gene", segment := i.seg_id]
  g[, c("mid", "mid2") := NULL]

  hg <- load_genes(P$human_genes); htx <- load_genes(P$human_tx)
  g[, gene_id_h := gene_ids_for(gene, hg, htx, if ("transcriptId" %chin% names(g)) transcriptId else NULL)]
  orth <- read_orthologs(P$orthologs)
  cls <- classify_orthologs(orth, g$gene_id_h)
  g[cls, on = "gene_id_h", `:=`(homology = i.homology, gene_id_m = i.gene_id_m, mouse_targets = i.mouse_targets)]
  g[is.na(homology), homology := "not_in_release"]
  log("genes %d; states %s; LOH %d; homology %s", nrow(g), paste(g[, .N, by = state][, paste0(state, "=", N)], collapse = ","),
      g[loh == TRUE, .N], paste(g[, .N, by = homology][, paste0(homology, "=", N)], collapse = ","))

  # genome-wide 1:1 orthologs with human coordinates (synteny blocks)
  o1 <- cls[homology == "ortholog_one2one", .(gene_id_h, gene_id_m)]
  o1[hg, on = c(gene_id_h = "gene_id"), `:=`(h_chrom = i.chrom, h_start = i.start)]

  # human segments cut into bins (coverage + cross-check)
  bs <- t$crosscheck_bin_bp
  bins <- segs[, {
    b0 <- seq(start - 1, end - 1, by = bs)
    .(start = b0, end = pmin(b0 + bs, end), k = seq_along(b0) - 1L)
  }, by = .(chrom = chromosome, seg_id, direction, loh)]
  bins[, name := paste0(seg_id, "|", k)]

  seg_rows <- list(); cross_rows <- list(); out_files <- character()
  for (build in names(P$mouse)) {
    M <- P$mouse[[build]]
    mg <- load_genes(M$genes)
    g <- attach_mouse_coords(g, mg, build)
    inb <- paste0(build, "_in_build"); bcol <- paste0(build, "_block")
    g[, (inb) := !is.na(gene_id_m) & gene_id_m %chin% mg$gene_id]

    oo <- o1[gene_id_m %chin% mg$gene_id & h_chrom %chin% STD_CHROMS]
    oo[mg, on = c(gene_id_m = "gene_id"), `:=`(m_chrom = i.chrom, m_start = i.start)]
    oo[, block := synteny_blocks(oo, cfg$synteny$max_rank_jump, prefix = paste0(build, "_B"))]
    g[oo, on = "gene_id_m", (bcol) := i.block]
    log("%s: %d synteny blocks from %d 1:1 orthologs", build, uniqueN(oo$block), nrow(oo))

    # ---- mapping coverage per human segment ------------------------------------------------
    lo <- run_liftover(bins, M$chain, cfg$liftover$min_match, P$tmp, paste0("cnbins_", build))
    bb <- lo$status[bins, on = "name"]
    lenf <- bb[, .(length_frac_mapped = sum((end - start) * (status == "mapped")) / sum(end - start)), by = seg_id]
    sg <- g[!is.na(segment), .(n_genes = .N, n_one2one = sum(homology == "ortholog_one2one"),
                               n_one2many = sum(homology %chin% c("ortholog_one2many", "ortholog_many2many")),
                               n_mapped_genes = sum(get(inb)),
                               mouse_chroms = {
                                 mc <- unique(get(paste0(build, "_chrom"))[get(inb)])
                                 paste(mc[order(chrom_order(mc))], collapse = ",")
                               },
                               n_blocks = uniqueN(na.omit(get(bcol)))), by = .(seg_id = segment)]
    sm <- sg[lenf[segs[, .(seg_id, copyNumber, state, direction, loh, length_bp = end - start + 1)], on = "seg_id"], on = "seg_id"]
    for (c in c("n_genes", "n_one2one", "n_one2many", "n_mapped_genes", "n_blocks")) set(sm, which(is.na(sm[[c]])), c, 0L)
    sm[is.na(mouse_chroms), mouse_chroms := ""]
    sm[is.na(length_frac_mapped), length_frac_mapped := 0]
    sm[, gene_frac_mapped := fifelse(n_genes > 0, n_mapped_genes / n_genes, NA_real_)]
    sm[, n_mouse_chroms := lengths(strsplit(mouse_chroms, ",", fixed = TRUE))]
    sm[, split_across_chromosomes := n_mouse_chroms > 1]
    sm[, confidence := segment_confidence(gene_frac_mapped, length_frac_mapped, n_one2one, n_one2many, t)]
    sm[, build := build]
    seg_rows[[build]] <- sm

    # ---- ortholog-based intervals ------------------------------------------------------------
    mm <- g[get(inb) == TRUE]
    track <- mm[, .(chrom = get(paste0(build, "_chrom")), start = get(paste0(build, "_start")), end = get(paste0(build, "_end")),
                    gene = gene_id_m, human_gene = gene_id_h, segment = fcoalesce(segment, "none"),
                    direction, lohlab = fifelse(loh, "LOH", "no"), block = get(bcol))]
    blk <- if (isTRUE(t$split_on_block)) "block" else NULL
    iv <- rbind(merge_runs(track, "direction", c("gain", "loss"), t$merge_max_gap_bp, blk),
                merge_runs(track, "lohlab", "LOH", t$merge_max_gap_bp, blk))
    state_of <- setNames(mm$state, mm$gene_id_m)
    sym_m <- setNames(mg$gene_name, mg$gene_id); sym_h <- setNames(hg$gene_name, hg$gene_id)
    conf_of <- setNames(sm$confidence, sm$seg_id)
    split_map <- function(s, f) vapply(strsplit(s, ",", fixed = TRUE), f, "")
    iv[, states := fifelse(direction == "LOH", "LOH", split_map(genes, function(x) paste(sort(unique(state_of[x])), collapse = ",")))]
    iv[, confidence := split_map(segments, function(x) worst_conf(fcoalesce(unname(conf_of[x]), "none")))]
    iv[, genes := split_map(genes, function(x) paste(fcoalesce(unname(sym_m[x]), x), collapse = ","))]
    iv[, human_genes := split_map(human_genes, function(x) paste(fcoalesce(unname(sym_h[x]), x), collapse = ","))]
    setorder(iv, chrom, start)
    iv <- iv[order(chrom_order(chrom), start)]
    iv[, name := sprintf("%s_CN%04d", build, seq_len(.N))]
    bed <- iv[, .(`#chrom` = chrom, start = start - 1L, end, name, direction, states, source_human_segments = segments,
                  genes, human_genes, block_ids = blocks, n_genes, confidence)]
    out_files <- c(out_files, write_tsv(bed, file.path(P$hp, sprintf("expected_mouse_cn.%s.bed", build)), log))
    log("%s intervals: %s", build, paste(iv[, .N, by = direction][, paste0(direction, "=", N)], collapse = ", "))

    # ---- cross-check: direct liftOver of the segments ---------------------------------------
    lifted <- lo$mapped[bins[, .(name, direction, loh)], on = "name", nomatch = NULL]
    for (d in c("gain", "loss", "LOH")) {
      lb <- if (d == "LOH") lifted[loh == TRUE] else lifted[direction == d]
      lb <- merge_overlapping(lb[, .(chrom, start, end)], gap = bs)
      ob <- iv[direction == d, .(chrom, start = start - 1L, end)]
      A <- total_bp(ob); B <- total_bp(lb); O <- overlap_bp(ob, lb)
      cross_rows[[paste(build, d)]] <- data.table(build = build, direction = d, ortholog_bp = A, liftover_bp = B, overlap_bp = O,
                                                  jaccard = if (A + B - O > 0) O / (A + B - O) else NA_real_,
                                                  frac_ortholog_in_liftover = if (A > 0) O / A else NA_real_,
                                                  frac_liftover_in_ortholog = if (B > 0) O / B else NA_real_)
    }
  }

  # ---- per-gene prediction -----------------------------------------------------------------
  hb <- cfg$h2m_mouse_build
  g[, prediction := fcase(homology %chin% c("no_ortholog", "not_in_release"), homology,
                          homology %chin% c("ortholog_one2many", "ortholog_many2many"), "one2many",
                          !get(paste0(hb, "_in_build")), "not_in_build", default = state)]
  cols <- c("gene", "gene_id_h", "chromosome", "start", "end", "minCopyNumber", "maxCopyNumber", "minMinorAlleleCopyNumber",
            "state", "direction", "partial", "loh", "segment", "homology", "gene_id_m", "mouse_targets",
            unlist(lapply(names(P$mouse), function(b) grep(paste0("^", b, "_"), names(g), value = TRUE))), "prediction")
  g <- g[order(chrom_order(chromosome), start)]
  f1 <- write_tsv(g[, cols, with = FALSE], file.path(P$hp, "expected_mouse_cn_genes.tsv"), log)
  f2 <- write_tsv(g[(direction %chin% c("gain", "loss") | loh) & homology %chin% c("ortholog_one2many", "ortholog_many2many"),
                    .(gene, gene_id_h, state, loh, segment, homology, mouse_targets)], file.path(P$hp, "cn_one2many.tsv"), log)
  f3 <- write_tsv(rbindlist(seg_rows, use.names = TRUE)[, .(build, segment = seg_id, copyNumber, state, direction, loh, length_bp,
                                                            n_genes, n_one2one, n_one2many, n_mapped_genes, gene_frac_mapped,
                                                            length_frac_mapped, mouse_chroms, n_mouse_chroms, n_blocks,
                                                            split_across_chromosomes, confidence)],
                  file.path(P$hp, "cn_segment_mapping.tsv"), log)
  f4 <- write_tsv(rbindlist(cross_rows), file.path(P$hp, "cn_liftover_crosscheck.tsv"), log)

  ev <- g[direction %chin% c("gain", "loss")]
  inb <- paste0(hb, "_in_build")
  rows <- rbind(srow("genes_assessed", g[state != "not_assessable", .N], "info"),
                srow("genes_not_assessable", g[state == "not_assessable", .N], "info"),
                g[, .N, by = state][order(state), .(category = paste0("state_", state), n = N, counterpart = "info")],
                srow("cn_event_genes", nrow(ev), "info"),
                srow("cn_event_one2one_in_build", ev[homology == "ortholog_one2one" & get(inb) == TRUE, .N], "yes"),
                srow("cn_event_one2many", ev[homology %chin% c("ortholog_one2many", "ortholog_many2many"), .N], "uncertain"),
                srow("cn_event_not_in_build", ev[homology == "ortholog_one2one" & get(inb) == FALSE, .N], "no"),
                srow("cn_event_no_ortholog", ev[homology %chin% c("no_ortholog", "not_in_release"), .N], "no"))
  write_summary(P, "cn", rows, log)
  # LOH: an ortholog exists, but inbred mice have no heterozygous SNPs -> uncertain (not assessable)
  lh <- g[loh == TRUE]
  write_summary(P, "loh", rbind(srow("loh_genes", nrow(lh), "info"),
                                srow("loh_one2one_in_build", lh[homology == "ortholog_one2one" & get(inb) == TRUE, .N], "uncertain"),
                                srow("loh_one2many", lh[homology %chin% c("ortholog_one2many", "ortholog_many2many"), .N], "uncertain"),
                                srow("loh_not_in_build", lh[homology == "ortholog_one2one" & get(inb) == FALSE, .N], "no"),
                                srow("loh_no_ortholog", lh[homology %chin% c("no_ortholog", "not_in_release"), .N], "no")), log)
  write_provenance(P, cfg, "03_cn_to_mouse", inputs = list(purple_cnv_gene = f_gene, purple_cnv_somatic = f_seg, purple_purity = f_pur),
                   resources = c(list(orthologs = P$orthologs), setNames(lapply(P$mouse, `[[`, "chain"), paste0(names(P$mouse), "_chain"))),
                   outputs = c(out_files, f1, f2, f3, f4), extra = list(ploidy = ploidy, thresholds = t))
  log("done")
}

main()
