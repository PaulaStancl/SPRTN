# ---------------------------------------------------------------------------
# 05b_pathways.R - every affected human gene -> COSMIC Cancer Gene Census role and Reactome
# pathways, with the expected functional direction of its events. Run after 01-04.
#
#   Rscript scripts/05b_pathways.R [--config ...] [--release ...]
#
# Affected genes: coding SNVs/indels (01), promoter / UTR / splice-region / TERT-promoter variants
# (02), AMP / GAIN / LOSS / DEEP_DEL and LOH genes (03), disrupted genes and fusion partners (04).
# CGC: config cosmic_cgc (downloaded by you - COSMIC licence), matched on Ensembl gene ID in
# 'Synonyms', then on 'Gene Symbol'. Without it roles stay empty and direction comes from the event.
# Reactome: lowest-level pathways (Ensembl2Reactome.txt) of the human gene (R-HSA) and of its 1:1
# mouse ortholog (R-MMU, for the later "related gene or pathway" comparison), plus top-level
# categories from ReactomePathwaysRelation.txt. Release = config reactome.version, else the newest
# folder 00_resources downloaded; recorded in provenance.
# Direction per event: inactivating (truncating, splice, LOSS, DEEP_DEL, LOH, disruption), activating
# (AMP, GAIN); MISSENSE / promoter -> activating for a CGC oncogene, inactivating for a CGC TSG,
# else unknown. Per gene: combined over its events (mixed if both); role_consistency = consistent /
# discordant with the CGC role.
#
# Output  human_predictions/
#   gene_pathways.tsv      gene_id_h symbol_h layers events expected_direction event_directions
#                          cgc_role cgc_tier role_consistency gene_id_m reactome_hsa_ids
#                          reactome_hsa_names reactome_top_level reactome_mmu_ids n_pathways
#   pathways_affected.tsv  pathway_id name top_level n_genes genes directions
#   summary/pathways.tsv
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

load_cgc <- function(path, log) {
  if (is.null(path)) {
    log("cosmic_cgc not set - CGC roles empty (download the Cancer Gene Census from COSMIC, set the path in config)", level = "WARNING")
    return(data.table(symbol = character(), ensg = character(), cgc_role = character(), cgc_tier = character()))
  }
  if (!file.exists(path)) stop("cosmic_cgc not found: ", path, call. = FALSE)
  d <- fread(path, colClasses = "character")
  # website export ('Gene Symbol', 'Role in Cancer', ...) or COSMIC download file (GENE_SYMBOL, ROLE_IN_CANCER, ...)
  setnames(d, gsub("[^a-z]", "", tolower(names(d))))
  need_cols(d, c("genesymbol", "roleincancer"), "COSMIC CGC (Gene Symbol / Role in Cancer)")
  syn <- if ("synonyms" %chin% names(d)) fcoalesce(d$synonyms, "") else rep("", nrow(d))
  m <- regexpr("ENSG[0-9]+", syn)
  ensg <- rep(NA_character_, length(syn))
  ensg[m > 0] <- regmatches(syn, m)
  data.table(symbol = d$genesymbol, ensg = ensg,
             cgc_role = d$roleincancer, cgc_tier = if ("tier" %chin% names(d)) d$tier else NA_character_)
}

# pathway -> all top-level ancestors (pathways without a parent)
make_top_levels <- function(rel) {
  parents <- split(rel$parent, rel$child)
  memo <- new.env()
  up <- function(pid, seen = character()) {
    if (!is.null(memo[[pid]])) return(memo[[pid]])
    ps <- setdiff(parents[[pid]], seen)
    res <- if (!length(ps)) pid else unique(unlist(lapply(ps, up, seen = c(seen, pid))))
    memo[[pid]] <- res
    res
  }
  up
}

main <- function() {
  a <- parse_cli()
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("05b_pathways", P)

  ev <- collect_gene_events(P)[!is.na(gene_id_h)]
  if (!nrow(ev)) stop("no events - run 01-04 first", call. = FALSE)
  log("affected genes: %d (%d events)", uniqueN(ev$gene_id_h), nrow(ev))

  cgc <- load_cgc(cfg$cosmic_cgc, log)
  want <- cfg$reactome$version
  rdir <- if (!is.null(want)) reactome_dir(P, want) else latest_reactome_dir(P)
  if (is.null(rdir) || !file.exists(file.path(rdir, "Ensembl2Reactome.txt"))) stop("Reactome files missing (", rdir, ") - run 00_resources.R", call. = FALSE)
  reactome_release <- as.integer(sub("^v", "", basename(rdir)))
  log("Reactome release %d (%s)", reactome_release, rdir)
  e2r <- fread(file.path(rdir, "Ensembl2Reactome.txt"), header = FALSE, sep = "\t", quote = "",
               col.names = c("ensembl_id", "pathway_id", "url", "name", "evidence", "species"))[species %chin% c("Homo sapiens", "Mus musculus")]
  rel <- fread(file.path(rdir, "ReactomePathwaysRelation.txt"), header = FALSE, col.names = c("parent", "child"))
  pnames <- fread(file.path(rdir, "ReactomePathways.txt"), header = FALSE, sep = "\t", quote = "", col.names = c("pathway_id", "name", "species"))
  pname <- setNames(pnames$name, pnames$pathway_id)
  up <- make_top_levels(rel)
  hsa <- e2r[species == "Homo sapiens", .(ids = list(sort(unique(pathway_id)))), by = ensembl_id]
  mmu <- e2r[species == "Mus musculus", .(ids = list(sort(unique(pathway_id)))), by = ensembl_id]

  cls <- classify_orthologs(read_orthologs(P$orthologs), ev$gene_id_h)

  gp <- ev[, {
    sym <- fcoalesce(na.omit(symbol_h)[1], gene_id_h)
    k <- match(gene_id_h, cgc$ensg)
    if (is.na(k)) k <- match(sym, cgc$symbol)
    role <- cgc$cgc_role[k]
    dirs <- vapply(event, event_direction, "", role = role)
    direction <- combine_directions(dirs)
    pw <- hsa$ids[match(gene_id_h, hsa$ensembl_id)][[1]]
    if (is.null(pw)) pw <- character()
    gm <- cls$gene_id_m[match(gene_id_h, cls$gene_id_h)]
    pm <- if (is.na(gm)) character() else { x <- mmu$ids[match(gm, mmu$ensembl_id)][[1]]; if (is.null(x)) character() else x }
    tops <- as.character(sort(unique(unlist(lapply(pw, up)))))
    ed <- unique(paste0(event, ":", dirs))
    .(symbol_h = sym, layers = paste(sort(unique(layer)), collapse = ","), events = paste(sort(unique(event)), collapse = ","),
      expected_direction = direction, event_directions = paste(sort(ed), collapse = ","),
      cgc_role = role, cgc_tier = cgc$cgc_tier[k], role_consistency = role_consistency(direction, role), gene_id_m = gm,
      reactome_hsa_ids = paste(pw, collapse = ";"), reactome_hsa_names = paste(fcoalesce(unname(pname[pw]), pw), collapse = ";"),
      reactome_top_level = paste(fcoalesce(unname(pname[tops]), tops), collapse = ";"),
      reactome_mmu_ids = paste(pm, collapse = ";"), n_pathways = length(pw))
  }, by = gene_id_h]
  setorder(gp, symbol_h)
  f1 <- write_tsv(gp, file.path(P$hp, "gene_pathways.tsv"), log)

  long <- gp[n_pathways > 0, .(pathway_id = unlist(strsplit(reactome_hsa_ids, ";", fixed = TRUE))), by = .(symbol_h, expected_direction)]
  pa <- long[, .(n_genes = uniqueN(symbol_h), genes = paste(sort(unique(symbol_h)), collapse = ","),
                 directions = paste(sort(unique(expected_direction)), collapse = ",")), by = pathway_id]
  pa[, `:=`(name = unname(pname[pathway_id]),
            top_level = vapply(pathway_id, function(x) { t <- up(x); paste(sort(fcoalesce(unname(pname[t]), t)), collapse = ";") }, ""))]
  setcolorder(pa, c("pathway_id", "name", "top_level"))
  setorder(pa, -n_genes, pathway_id)
  f2 <- write_tsv(pa, file.path(P$hp, "pathways_affected.tsv"), log)

  write_summary(P, "pathways", rbind(srow("genes_affected", nrow(gp), "info"), srow("genes_with_cgc_role", gp[!is.na(cgc_role), .N], "info"),
                                     srow("genes_in_reactome", gp[n_pathways > 0, .N], "info"), srow("pathways_hit", nrow(pa), "info"),
                                     gp[, .N, by = expected_direction][order(expected_direction),
                                        .(category = paste0("direction_", expected_direction), n = N, counterpart = "info")]), log)
  write_provenance(P, cfg, "05b_pathways", inputs = list(cosmic_cgc = cfg$cosmic_cgc),
                   resources = list(reactome_ensembl2reactome = file.path(rdir, "Ensembl2Reactome.txt")), outputs = c(f1, f2),
                   extra = list(reactome_release = reactome_release, cosmic_cgc_version = cfg$cosmic_cgc_version))
  log("done")
}

main()
