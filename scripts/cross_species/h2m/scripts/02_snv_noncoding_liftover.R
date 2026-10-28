# ---------------------------------------------------------------------------
# 02_snv_noncoding_liftover.R - non-coding PASS SNVs/indels lifted to each mouse build.
#
#   Rscript scripts/02_snv_noncoding_liftover.R [--config ...] [--release ...]
#
# Everything 01 does not model (PAVE canonical coding effect not in snv.coding_effects) is classed,
# first match wins: TERT_promoter (config window), promoter (TSS - promoter_upstream_bp .. TSS +
# promoter_downstream_bp of an Ensembl canonical transcript, or PAVE upstream_gene_variant),
# splice_region (PAVE splice-region flag / effect), UTR5, UTR3, intronic, other.
# The REF span is lifted with liftOver -minMatch=<liftover.min_match> -multiple, so unmapped and
# multi-mapped are counted separately. A single hit gets the mouse REF from the build's FASTA
# (reverse-complemented when the chain flips strand): ref_conserved. in_ortholog_promoter /
# in_ortholog_gene: the hit lies in the promoter window / body of a mouse ortholog of the human
# gene (consistency check of the liftOver).
#
# Input   <oncoanalyser>/purple/*.purple.somatic.vcf.gz
# Output  human_predictions/
#   expected_mouse_noncoding.tsv  variant_id chrom pos ref alt vtype class gene_symbol_h gene_id_h
#     pave_effect reported homology_type gene_id_m, then per build <b>_status (mapped / multi_mapped /
#     unmapped) <b>_n_hits <b>_reason <b>_chrom <b>_pos <b>_strand <b>_ref <b>_ref_conserved
#     <b>_in_ortholog_promoter <b>_in_ortholog_gene; counterpart
#   expected_mouse_noncoding_multihits.tsv  every hit of multi-mapped variants
#   summary/snv_noncoding.tsv
# counterpart (H2M build): yes = mapped once with the same REF; uncertain = other REF or
# multi-mapped; no = unmapped.
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

promoter_windows <- function(tx, up, down) {
  t <- tx[canonical == TRUE]
  t[, tss := fifelse(strand == "+", start, end)]
  t[, `:=`(w_start = fifelse(strand == "+", tss - up, tss - down), w_end = fifelse(strand == "+", tss + down, tss + up))]
  t[, .(chrom, w_start, w_end, gene_id, transcript_id)]
}

in_windows <- function(chrom, pos, win) {
  q <- data.table(chrom = chrom, s = pos, e = pos, i = seq_along(pos))
  w <- copy(win)[, `:=`(s = w_start, e = w_end)]
  setkey(w, chrom, s, e)
  o <- foverlaps(q, w, by.x = c("chrom", "s", "e"), nomatch = NULL)
  seq_along(pos) %in% o$i
}

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("02_snv_noncoding_liftover", P)
  nc <- cfg$noncoding
  vcf <- find_one(file.path(cfg$oncoanalyser_dir, "purple"), "\\.purple\\.somatic\\.vcf\\.gz$", "PURPLE somatic VCF", log)
  v <- read_purple_somatic(vcf, log)[!pave_coding_effect %chin% cfg$snv$coding_effects]
  log("non-coding: %d variants", nrow(v))

  hg <- load_genes(P$human_genes); htx <- load_genes(P$human_tx)
  tert <- parse_region(nc$tert_promoter)
  prom <- in_windows(v$chrom, v$pos, promoter_windows(htx, nc$promoter_upstream_bp, nc$promoter_downstream_bp))
  eff <- tolower(fcoalesce(v$pave_effect, ""))
  v[, class := fcase(chrom == tert$chrom & pos >= tert$start & pos <= tert$end, "TERT_promoter",
                     prom | grepl("upstream", eff), "promoter",
                     splice_region | grepl("splice", eff), "splice_region",
                     grepl("5_prime_utr", eff), "UTR5", grepl("3_prime_utr", eff), "UTR3",
                     grepl("intron", eff), "intronic", default = "other")]
  v[, gene_id_h := gene_ids_for(gene_symbol_h, hg, htx, transcript_h)]
  log("classes: %s", paste(v[, .N, by = class][, paste0(class, "=", N)], collapse = ", "))

  orth <- read_orthologs(P$orthologs)
  cls <- classify_orthologs(orth, v$gene_id_h)
  v[cls, on = "gene_id_h", `:=`(homology_type = i.homology, mouse_targets = i.mouse_targets, gene_id_m = i.gene_id_m)]

  bed <- v[, .(chrom, start = pos - 1L, end = pos - 1L + nchar(ref), name = variant_id, strand = "+")]
  multi_all <- list()
  for (build in names(P$mouse)) {
    M <- P$mouse[[build]]
    lo <- run_liftover(bed, M$chain, cfg$liftover$min_match, P$tmp, paste0("noncoding_", build), multiple = TRUE)
    st <- lo$status
    v[st, on = c(variant_id = "name"), c("status", "n_hits", "reason") := .(i.status, i.n_hits, i.reason)]
    one <- lo$mapped[name %chin% st[status == "mapped", name]]
    v[one, on = c(variant_id = "name"), c("m_chrom", "m_pos", "m_strand") := .(i.chrom, i.start + 1L, i.strand)]
    if (!"m_chrom" %chin% names(v)) v[, c("m_chrom", "m_pos", "m_strand") := .(NA_character_, NA_integer_, NA_character_)]
    multi_all[[build]] <- lo$mapped[name %chin% st[status == "multi_mapped", name]][, build := build]

    # mouse REF at a single hit
    if (file.exists(M$fasta) && file.exists(paste0(M$fasta, ".fai"))) {
      fa <- Rsamtools::FaFile(M$fasta)
      seqs <- GenomeInfoDb::seqnames(Rsamtools::seqinfo(fa))
      ok <- which(!is.na(v$m_chrom))
      nm <- fifelse(v$m_chrom[ok] %chin% seqs, v$m_chrom[ok], sub("^chr", "", v$m_chrom[ok]))
      keep <- nm %chin% seqs
      v[, m_ref := NA_character_]
      if (any(keep)) {
        gr <- GenomicRanges::GRanges(nm[keep], IRanges::IRanges(v$m_pos[ok][keep], v$m_pos[ok][keep] + nchar(v$ref[ok][keep]) - 1L))
        s <- toupper(as.character(Rsamtools::scanFa(fa, gr)))
        flip <- v$m_strand[ok][keep] == "-"
        s[flip] <- revcomp(s[flip])
        set(v, i = ok[keep], j = "m_ref", value = s)
      }
      v[, m_ref_conserved := fifelse(is.na(m_ref), NA, m_ref == ref)]
    } else {
      log("%s FASTA missing - no mouse REF check", build, level = "WARNING")
      v[, c("m_ref", "m_ref_conserved") := .(NA_character_, NA)]
    }

    # does the hit fall in the promoter / body of a mouse ortholog of the human gene?
    mg <- load_genes(M$genes); mtx <- load_genes(M$tx)
    pw <- promoter_windows(mtx, nc$promoter_upstream_bp, nc$promoter_downstream_bp)
    tg <- v[!is.na(m_chrom) & !is.na(mouse_targets), .(gene_id = unlist(strsplit(mouse_targets, ";", fixed = TRUE))),
            by = .(variant_id, m_chrom, m_pos)]
    hit_p <- pw[tg, on = "gene_id", nomatch = NULL][chrom == m_chrom & m_pos >= w_start & m_pos <= w_end, unique(variant_id)]
    hit_g <- mg[tg, on = "gene_id", nomatch = NULL][chrom == m_chrom & m_pos >= start & m_pos <= end, unique(variant_id)]
    v[, `:=`(m_in_ortholog_promoter = fifelse(is.na(m_chrom) | is.na(mouse_targets), NA, variant_id %chin% hit_p),
             m_in_ortholog_gene = fifelse(is.na(m_chrom) | is.na(mouse_targets), NA, variant_id %chin% hit_g))]

    tmp <- c("status", "n_hits", "reason", "m_chrom", "m_pos", "m_strand", "m_ref", "m_ref_conserved",
             "m_in_ortholog_promoter", "m_in_ortholog_gene")
    setnames(v, tmp, paste0(build, "_", c("status", "n_hits", "reason", "chrom", "pos", "strand", "ref", "ref_conserved",
                                         "in_ortholog_promoter", "in_ortholog_gene")))
    log("%s: %s", build, paste(v[, .N, by = c(paste0(build, "_status"))][, paste0(get(paste0(build, "_status")), "=", N)], collapse = ", "))
  }

  hb <- cfg$h2m_mouse_build
  v[, counterpart := fcase(get(paste0(hb, "_status")) == "unmapped", "no",
                           get(paste0(hb, "_status")) == "mapped" & get(paste0(hb, "_ref_conserved")) %in% TRUE, "yes",
                           default = "uncertain")]
  keep <- c("variant_id", "chrom", "pos", "ref", "alt", "vtype", "class", "gene_symbol_h", "gene_id_h", "pave_effect",
            "reported", "homology_type", "gene_id_m")
  per_build <- unlist(lapply(names(P$mouse), function(b) grep(paste0("^", b, "_"), names(v), value = TRUE)))
  f <- write_tsv(v[, c(keep, per_build, "counterpart"), with = FALSE], file.path(P$hp, "expected_mouse_noncoding.tsv"), log)
  fm <- write_tsv(rbindlist(multi_all, use.names = TRUE, fill = TRUE), file.path(P$hp, "expected_mouse_noncoding_multihits.tsv"), log)

  rows <- rbind(srow("variants_total", nrow(v), "info"),
                v[, .N, by = class][order(class), .(category = paste0("class_", class), n = N, counterpart = "info")],
                rbindlist(lapply(names(P$mouse), function(b) v[, .N, by = c(paste0(b, "_status"))][
                  , .(category = paste0(b, "_", get(paste0(b, "_status"))), n = N, counterpart = "info")])),
                v[, .N, by = counterpart][order(counterpart), .(category = paste0("counterpart_", counterpart), n = N, counterpart)])
  write_summary(P, "snv_noncoding", rows, log)
  write_provenance(P, cfg, "02_snv_noncoding_liftover", inputs = list(purple_somatic_vcf = vcf),
                   resources = setNames(lapply(P$mouse, `[[`, "chain"), paste0(names(P$mouse), "_chain")),
                   outputs = c(f, fm), extra = list(liftover_min_match = cfg$liftover$min_match))
  log("done")
}

main()
