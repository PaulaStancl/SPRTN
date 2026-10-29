# ---------------------------------------------------------------------------
# 04_sv_to_mouse.R - LINX SVs, disrupted genes and fusions -> mouse.
#
#   Rscript scripts/04_sv_to_mouse.R [--config ...] [--release ...]
#
# Breakends from linx.vis_sv_data.tsv (SvId, ClusterId, Type, ResolvedType, Chr/Pos/Orient Start/End
# - the columns analysis/run_chromothripsis.R already uses); orientation 1 -> '+', -1 -> '-'.
# Each breakend +/- sv.breakend_flank_bp is lifted (-minMatch=liftover.min_match); the chain can
# flip the strand, which flips the orientation in mouse. SVs in sv.chromothripsis_regions (config)
# are flagged. Genes hit by a breakend (linx.breakend.tsv) go to mouse via Ensembl orthologs - the
# prediction is "disruption of the ortholog", never the breakpoint coordinates.
#
# Fusions (linx.fusion.tsv): feasibility is a FLAG, not a prediction - no mouse fusion is ever
# predicted because both partners have orthologs (mouse_fusion_predicted is always FALSE).
# sv_class = the SV type that would join the 5' and 3' genes given their positions and strands
# (TRA / INV / DEL / DUP), in human and in each mouse build; fusion_feasibility =
# not_possible_no_ortholog / feasible_same_sv_class / feasible_needs_<class> / not_assessable.
#
# Input   <oncoanalyser>/linx/*.linx.{vis_sv_data,clusters,breakend,fusion,driver.catalog}.tsv
# Output  human_predictions/
#   expected_mouse_sv.<build>.bedpe  chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2
#       (lifted +/- flank windows; '.' / -1 when unmapped) sv_type resolved_type cluster_id
#       cluster_resolved_type chromothripsis human_locus1 human_locus2 status1 status2
#       same_mouse_chrom genes1_h genes2_h
#   expected_mouse_sv_genes.tsv      per human gene hit: sv_ids, event (HOM_DISRUPTION / DISRUPTION /
#       FUSION_PARTNER / NON_DISRUPTIVE), reported, region/coding type, homology, mouse gene(s) +
#       coords per build, prediction, counterpart
#   expected_mouse_fusions.tsv       name reported reportedType gene5 gene3 (+IDs) human_sv_class
#       n_mouse5 n_mouse3, per build <b>_sv_class <b>_fusion_feasibility; mouse_fusion_predicted, note
#   summary/sv.tsv
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("04_sv_to_mouse", P)
  ld <- file.path(cfg$oncoanalyser_dir, "linx")
  F <- list(vis = find_one(ld, "\\.linx\\.vis_sv_data\\.tsv$", "LINX vis_sv_data", log),
            clusters = find_one(ld, "\\.linx\\.clusters\\.tsv$", "LINX clusters", log, required = FALSE),
            breakend = find_one(ld, "\\.linx\\.breakend\\.tsv$", "LINX breakends", log),
            fusion = find_one(ld, "\\.linx\\.fusion\\.tsv$", "LINX fusions", log, required = FALSE),
            driver = find_one(ld, "\\.linx\\.driver\\.catalog\\.tsv$", "LINX driver catalog", log, required = FALSE))
  flank <- cfg$sv$breakend_flank_bp
  regions <- lapply(cfg$sv$chromothripsis_regions, parse_region)

  vis <- fread(F$vis, colClasses = list(character = c("ChrStart", "ChrEnd")))
  need_cols(vis, c("ClusterId", "SvId", "Type", "ResolvedType", "ChrStart", "ChrEnd", "PosStart", "PosEnd",
                   "OrientStart", "OrientEnd"), "linx.vis_sv_data.tsv")
  vis <- vis[!duplicated(SvId)]
  vis[, c1 := norm_chrom(ChrStart)]
  vis[, has2 := !ChrEnd %chin% c("-1", "0") & !is.na(ChrEnd) & !Type %chin% c("SGL", "INF")]
  vis[, c2 := fifelse(has2, norm_chrom(ChrEnd), NA_character_)]
  in_ct <- function(ch, p) {
    hit <- rep(FALSE, length(ch))
    for (r in regions) hit <- hit | (!is.na(ch) & ch == r$chrom & p >= r$start & p <= r$end)
    hit
  }
  vis[, chromothripsis := in_ct(c1, PosStart) | in_ct(c2, PosEnd)]
  vis[, cluster_resolved_type := NA_character_]
  if (!is.null(F$clusters)) {
    clu <- fread(F$clusters)
    if (all(c("clusterId", "resolvedType") %chin% names(clu)))
      vis[clu[!duplicated(clusterId)], on = c(ClusterId = "clusterId"), cluster_resolved_type := i.resolvedType]
  }
  log("SVs: %d (%s); in chromothripsis regions: %d", nrow(vis), paste(vis[, .N, by = Type][, paste0(Type, "=", N)], collapse = ","),
      vis[chromothripsis == TRUE, .N])

  # ---- genes at breakends --------------------------------------------------------------------
  be <- fread(F$breakend)
  need_cols(be, c("svId", "isStart", "gene"), "linx.breakend.tsv")
  for (c in c("disruptive", "reportedDisruption")) set(be, j = c, value = if (c %chin% names(be)) truthy(be[[c]]) else FALSE)
  for (c in c("regionType", "codingType")) if (!c %chin% names(be)) set(be, j = c, value = NA_character_)
  be <- be[!is.na(gene) & gene != ""]
  be[, end_tag := fifelse(truthy(isStart), 1L, 2L)]
  genes_end <- be[, .(genes = paste(sort(unique(gene)), collapse = ",")), by = .(svId, end_tag)]

  hg <- load_genes(P$human_genes); htx <- load_genes(P$human_tx)
  orth <- read_orthologs(P$orthologs)

  drv <- if (!is.null(F$driver)) fread(F$driver) else data.table(gene = character(), driver = character())
  drv_of <- drv[, .(linx_driver = paste(sort(unique(driver)), collapse = ",")), by = gene]

  fus <- if (!is.null(F$fusion)) fread(F$fusion) else data.table()
  if (nrow(fus)) {
    need_cols(fus, "name", "linx.fusion.tsv")
    parts <- tstrsplit(fus$name, "_", fixed = TRUE, keep = 1:2)
    fus[, gene5 := if ("geneStart" %chin% names(fus)) as.character(geneStart) else parts[[1]]]
    fus[, gene3 := if ("geneEnd" %chin% names(fus)) as.character(geneEnd) else parts[[2]]]
    fus[, reported := if ("reported" %chin% names(fus)) truthy(reported) else FALSE]
  }
  fusion_genes <- if (nrow(fus)) unique(c(fus$gene5, fus$gene3)) else character()

  be[, gene_id_h := gene_ids_for(gene, hg, htx, if ("transcriptId" %chin% names(be)) transcriptId else NULL)]
  ge <- be[, .(sv_ids = paste(sort(unique(svId)), collapse = ","), disruptive = any(disruptive),
               reported_disruption = any(reportedDisruption),
               region_type = paste(sort(unique(na.omit(as.character(regionType)))), collapse = ","),
               coding_type = paste(sort(unique(na.omit(as.character(codingType)))), collapse = ",")),
           by = .(gene, gene_id_h)]
  extra <- setdiff(fusion_genes, ge$gene)
  if (length(extra))
    ge <- rbind(ge, data.table(gene = extra, gene_id_h = gene_ids_for(extra, hg, htx), sv_ids = "", disruptive = FALSE,
                               reported_disruption = FALSE, region_type = "", coding_type = ""))
  ge[drv_of, on = "gene", linx_driver := i.linx_driver]
  if (!"linx_driver" %chin% names(ge)) ge[, linx_driver := NA_character_]
  ge[, fusion_partner := gene %chin% fusion_genes]
  ge[, event := fcase(grepl("HOM_", fcoalesce(linx_driver, "")), "HOM_DISRUPTION", disruptive, "DISRUPTION",
                      fusion_partner, "FUSION_PARTNER", default = "NON_DISRUPTIVE")]
  ge[, reported := reported_disruption | !is.na(linx_driver)]

  cls <- classify_orthologs(orth, ge$gene_id_h)
  ge[cls, on = "gene_id_h", `:=`(homology = i.homology, n_mouse = i.n_mouse, mouse_targets = i.mouse_targets, gene_id_m = i.gene_id_m)]
  ge[is.na(homology), `:=`(homology = "not_in_release", n_mouse = 0L)]

  if (nrow(fus)) {
    for (side in c("5", "3")) {
      gid <- gene_ids_for(fus[[paste0("gene", side)]], hg, htx)
      set(fus, j = paste0("gene_id", side, "_h"), value = gid)
      set(fus, j = paste0("n_mouse", side), value = fcoalesce(cls$n_mouse[match(gid, cls$gene_id_h)], 0L))
      set(fus, j = paste0("gene_id", side, "_m"), value = cls$gene_id_m[match(gid, cls$gene_id_h)])
      for (c in c("chrom", "start", "strand")) set(fus, j = paste0("h", side, "_", c), value = hg[[c]][match(gid, hg$gene_id)])
    }
    fus[, human_sv_class := sv_class(h5_chrom, h5_start, h5_strand, h3_chrom, h3_start, h3_strand)]
  }

  hb <- cfg$primary_mouse_build      # main build for counterpart calls, event table, report
  out_files <- character(); sv_rows <- list()
  for (build in names(P$mouse)) {
    M <- P$mouse[[build]]
    mg <- load_genes(M$genes)
    ge <- attach_mouse_coords(ge, mg, build)

    ends <- rbind(vis[, .(chrom = c1, pos = PosStart, name = paste0(SvId, "_1"), strand = fifelse(OrientStart == 1, "+", "-"))],
                  vis[has2 == TRUE, .(chrom = c2, pos = PosEnd, name = paste0(SvId, "_2"), strand = fifelse(OrientEnd == 1, "+", "-"))])
    ends[, `:=`(start = pmax(0, pos - 1 - flank), end = pos + flank)]
    lo <- run_liftover(ends, M$chain, cfg$liftover$min_match, P$tmp, paste0("sv_", build))
    mp <- lo$mapped[!duplicated(name)]
    end_cols <- function(k) {
      n <- paste0(vis$SvId, "_", k)
      i <- match(n, mp$name)
      list(chrom = fcoalesce(mp$chrom[i], "."), start = fcoalesce(as.integer(mp$start[i]), -1L), end = fcoalesce(as.integer(mp$end[i]), -1L),
           strand = fcoalesce(mp$strand[i], "."), status = fcoalesce(lo$status$status[match(n, lo$status$name)], "no_breakend"))
    }
    e1 <- end_cols(1); e2 <- end_cols(2)
    g1 <- genes_end[end_tag == 1L][match(vis$SvId, svId), fcoalesce(genes, "")]
    g2 <- genes_end[end_tag == 2L][match(vis$SvId, svId), fcoalesce(genes, "")]
    bedpe <- data.table(`#chrom1` = e1$chrom, start1 = e1$start, end1 = e1$end, chrom2 = e2$chrom, start2 = e2$start, end2 = e2$end,
                        name = paste0("SV", vis$SvId), score = ".", strand1 = e1$strand, strand2 = e2$strand,
                        sv_type = vis$Type, resolved_type = vis$ResolvedType, cluster_id = vis$ClusterId,
                        cluster_resolved_type = vis$cluster_resolved_type, chromothripsis = vis$chromothripsis,
                        human_locus1 = paste(vis$c1, vis$PosStart, vis$OrientStart, sep = ":"),
                        human_locus2 = fifelse(vis$has2, paste(vis$c2, vis$PosEnd, vis$OrientEnd, sep = ":"), "."),
                        status1 = e1$status, status2 = e2$status,
                        same_mouse_chrom = fifelse(e1$chrom != "." & e2$chrom != ".", e1$chrom == e2$chrom, NA),
                        genes1_h = g1, genes2_h = g2)
    out_files <- c(out_files, write_tsv(bedpe, file.path(P$hp, sprintf("expected_mouse_sv.%s.bedpe", build)), log))
    st <- c(e1$status, e2$status[vis$has2])
    sv_rows[[build]] <- data.table(status = st)[, .N, by = status][order(status), .(category = paste0(build, "_breakends_", status), n = N, counterpart = "info")]

    if (nrow(fus)) {
      for (side in c("5", "3")) {
        gm <- fus[[paste0("gene_id", side, "_m")]]
        for (c in c("chrom", "start", "strand")) set(fus, j = paste0(build, "_m", side, "_", c), value = mg[[c]][match(gm, mg$gene_id)])
      }
      fus[, (paste0(build, "_sv_class")) := sv_class(get(paste0(build, "_m5_chrom")), get(paste0(build, "_m5_start")), get(paste0(build, "_m5_strand")),
                                                     get(paste0(build, "_m3_chrom")), get(paste0(build, "_m3_start")), get(paste0(build, "_m3_strand")))]
      fus[, (paste0(build, "_fusion_feasibility")) := fusion_feasibility(n_mouse5, n_mouse3, human_sv_class, get(paste0(build, "_sv_class")))]
    }
  }

  ge[, prediction := fcase(homology %chin% c("no_ortholog", "not_in_release"), homology,
                           event == "NON_DISRUPTIVE", "no_prediction_non_disruptive",
                           homology != "ortholog_one2one", "one2many",
                           is.na(get(paste0(hb, "_chrom"))), "not_in_build", default = "disruption_of_ortholog")]
  cp <- c(disruption_of_ortholog = "yes", one2many = "uncertain", no_prediction_non_disruptive = "info")
  ge[, counterpart := fcoalesce(unname(cp[prediction]), "no")]
  cols <- c("gene", "gene_id_h", "sv_ids", "event", "disruptive", "reported", "linx_driver", "fusion_partner", "region_type",
            "coding_type", "homology", "n_mouse", "gene_id_m", "mouse_targets",
            unlist(lapply(names(P$mouse), function(b) grep(paste0("^", b, "_"), names(ge), value = TRUE))), "prediction", "counterpart")
  f_g <- write_tsv(ge[order(event, gene), cols, with = FALSE], file.path(P$hp, "expected_mouse_sv_genes.tsv"), log)

  if (nrow(fus)) {
    fus[, `:=`(mouse_fusion_predicted = FALSE, note = "feasibility flag only - no mouse fusion is predicted")]
    keep <- c("name", "reported", "reportedType", "gene5", "gene3", "gene_id5_h", "gene_id3_h", "human_sv_class", "n_mouse5",
              "n_mouse3", "gene_id5_m", "gene_id3_m",
              unlist(lapply(names(P$mouse), function(b) grep(paste0("^", b, "_"), names(fus), value = TRUE))),
              "mouse_fusion_predicted", "note")
    fus <- fus[, intersect(keep, names(fus)), with = FALSE]
  }
  f_f <- write_tsv(fus, file.path(P$hp, "expected_mouse_fusions.tsv"), log)

  ev <- ge[event != "NON_DISRUPTIVE"]
  rows <- rbind(srow("svs_total", nrow(vis), "info"), srow("svs_in_chromothripsis_regions", vis[chromothripsis == TRUE, .N], "info"),
                rbindlist(sv_rows), srow("genes_hit_total", nrow(ge), "info"),
                ev[, .N, by = prediction][order(prediction), .(category = paste0("event_gene_", prediction), n = N,
                                                               counterpart = fcoalesce(unname(cp[prediction]), "no"))])
  if (nrow(fus))
    rows <- rbind(rows, srow("fusions_total", nrow(fus), "info"),
                  fus[, .N, by = c(paste0(hb, "_fusion_feasibility"))][, .(category = paste0("fusion_", hb, "_", get(paste0(hb, "_fusion_feasibility"))),
                                                                          n = N, counterpart = "info")])
  write_summary(P, "sv", rows, log)
  write_provenance(P, cfg, "04_sv_to_mouse", inputs = F,
                   resources = c(list(orthologs = P$orthologs), setNames(lapply(P$mouse, `[[`, "chain"), paste0(names(P$mouse), "_chain"))),
                   outputs = c(out_files, f_g, f_f), extra = list(breakend_flank_bp = flank, chromothripsis_regions = cfg$sv$chromothripsis_regions))
  log("done")
}

main()
