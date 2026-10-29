# ---------------------------------------------------------------------------
# common.R - shared helpers for cross_species/h2m, sourced by every script:
# config, paths, CLI, logging, provenance, GTF / ortholog tables, PURPLE VCF (bcftools),
# liftOver (UCSC binary), interval arithmetic, CN states, synteny blocks, SV classes,
# trinucleotide contexts, gene-level events.
#
# The pure functions (no I/O) are covered by tests/testthat: classify_orthologs, map_orthologs,
# cn_state, gene_cn_state, synteny_blocks, merge_runs, merge_overlapping, overlap_bp,
# sv_class, fusion_feasibility, vcf_to_maf, collapse_trinuc, sbs96_label, rescale_spectrum,
# event_direction, role_consistency, combine_directions, short_protein.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(jsonlite)
})

STD_CHROMS <- c(paste0("chr", 1:22), "chrX", "chrY")
CLOUD_RE <- "(^|/)(OneDrive[^/]*|CloudStorage|Dropbox|iCloud Drive|Mobile Documents|Google Drive)(/|$)"   # whole path components
HOMOLOGY_ORDER <- c("ortholog_one2one", "ortholog_one2many", "ortholog_many2many")

script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
H2M_SCRIPTS <- script_dir()
DEFAULT_CONFIG <- file.path(dirname(H2M_SCRIPTS), "config", "config.yaml")

# ============================================================================
# CLI, config, paths, logging, provenance
# ============================================================================

# --key value and --flag (flags listed in `flags`); returns a named list
parse_cli <- function(flags = character(), args = commandArgs(TRUE)) {
  out <- list(config = DEFAULT_CONFIG, release = NULL)
  i <- 1L
  while (i <= length(args)) {
    k <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--")) stop("unexpected argument: ", args[i], call. = FALSE)
    if (k %chin% flags) { out[[gsub("-", "_", k)]] <- TRUE; i <- i + 1L; next }
    if (i == length(args)) stop("missing value for ", args[i], call. = FALSE)
    out[[gsub("-", "_", k)]] <- args[i + 1L]; i <- i + 2L
  }
  for (f in flags) if (is.null(out[[gsub("-", "_", f)]])) out[[gsub("-", "_", f)]] <- FALSE
  out
}

subst <- function(x, scalars) {
  if (is.character(x)) {
    for (k in seq_len(5)) {
      y <- x
      for (n in names(scalars)) y <- gsub(paste0("{", n, "}"), scalars[[n]], y, fixed = TRUE)
      if (identical(y, x)) break
      x <- y
    }
    return(x)
  }
  if (is.list(x)) return(lapply(x, subst, scalars = scalars))
  x
}

latest_release_file <- function(resources_dir) file.path(resources_dir, "ensembl_latest_release.txt")

query_latest_release <- function() {
  r <- fromJSON("https://rest.ensembl.org/info/data/?content-type=application/json")
  as.integer(max(r$releases))
}

resolve_release <- function(tag, resources_dir, query = FALSE) {
  if (as.character(tag) != "latest") return(as.integer(tag))
  f <- latest_release_file(resources_dir)
  if (query) {
    n <- query_latest_release()
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    writeLines(as.character(n), f)
    return(n)
  }
  if (!file.exists(f)) stop("Ensembl release 'latest' not resolved yet - run 00_resources.R --release latest (", f, ")", call. = FALSE)
  as.integer(readLines(f, n = 1))
}

load_config <- function(path, release = NULL, query_latest = FALSE) {
  cfg <- read_yaml(path)
  if (!is.null(release)) cfg$ensembl_release <- release
  cfg$ensembl_release_tag <- as.character(cfg$ensembl_release)
  cfg$ensembl_release <- resolve_release(cfg$ensembl_release, cfg$resources_dir, query = query_latest)
  scalars <- Filter(function(v) length(v) == 1 && (is.character(v) || is.numeric(v)), cfg)
  scalars <- lapply(scalars, as.character)
  cfg <- subst(cfg, scalars)
  for (k in c("out_dir", "resources_dir"))
    if (grepl(CLOUD_RE, normalizePath(cfg[[k]], mustWork = FALSE), ignore.case = TRUE))
      stop("refusing to write patient-derived data to a cloud-synced folder: ", cfg[[k]], call. = FALSE)
  cfg$config_path <- normalizePath(path)
  cfg
}

ucsc_chain_name <- function(src, dst) paste0(src, "To", toupper(substr(dst, 1, 1)), substring(dst, 2), ".over.chain.gz")

paths <- function(cfg) {
  R <- cfg$ensembl_release
  rd <- cfg$resources_dir
  out <- cfg$out_dir
  rel <- file.path(rd, paste0("ensembl", R))
  P <- list(
    resources = rd, release_dir = rel,
    orthologs = file.path(rel, "orthologs_hsapiens_mmusculus.tsv"),
    human_gtf = file.path(rel, sprintf("Homo_sapiens.GRCh38.%s.gtf.gz", R)),
    human_genes = file.path(rel, "human_genes.tsv.gz"),
    human_tx = file.path(rel, "human_transcripts.tsv.gz"),
    chains_dir = file.path(rd, "chains"),
    mm39_to_mm10 = file.path(rd, "chains", ucsc_chain_name("mm39", "mm10")),
    trinuc_dir = file.path(rd, "trinuc"),
    out = out, hp = file.path(out, "human_predictions"),
    summary = file.path(out, "human_predictions", "summary"),
    frozen = file.path(out, "human_predictions", "FROZEN"),
    comparison = file.path(out, "comparison"),
    tmp = file.path(out, "tmp"), logs = file.path(out, "logs"),
    provenance = file.path(out, "provenance.json"),
    mouse = list())
  for (build in names(cfg$mouse_builds)) {
    b <- cfg$mouse_builds[[build]]
    gr <- if (as.character(b$gtf_release) == "current") R else as.integer(b$gtf_release)
    gdir <- file.path(rd, paste0("ensembl", gr))
    P$mouse[[build]] <- list(
      ucsc = b$ucsc, gtf_release = gr,
      gtf = file.path(gdir, sprintf("Mus_musculus.%s.%s.gtf.gz", build, gr)),
      genes = file.path(gdir, sprintf("mouse_%s_genes.tsv.gz", build)),
      tx = file.path(gdir, sprintf("mouse_%s_transcripts.tsv.gz", build)),
      fasta = file.path(rd, "fasta", sprintf("Mus_musculus.%s.dna.primary_assembly.fa", build)),
      chain = file.path(rd, "chains", ucsc_chain_name("hg38", b$ucsc)))
  }
  P
}

reactome_dir <- function(P, version) file.path(P$resources, "reactome", paste0("v", as.integer(version)))

latest_reactome_dir <- function(P) {
  d <- list.dirs(file.path(P$resources, "reactome"), recursive = FALSE)
  d <- d[grepl("^v[0-9]+$", basename(d))]
  if (!length(d)) return(NULL)
  d[which.max(as.integer(sub("^v", "", basename(d))))]
}

make_dirs <- function(P) for (d in c(P$hp, P$summary, P$tmp, P$logs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# log("fmt", ...) -> stderr and logs/<name>.log
make_logger <- function(name, P) {
  f <- file.path(P$logs, paste0(name, ".log"))
  cat("", file = f)
  function(fmt, ..., level = "INFO") {
    line <- sprintf("%s %s %s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), name, level, sprintf(fmt, ...))
    message(line)
    cat(line, "\n", file = f, append = TRUE, sep = "")
  }
}

sha256 <- function(path, limit_gb = 2) {
  if (is.null(path) || !file.exists(path) || file.size(path) > limit_gb * 1e9) return(NULL)
  digest::digest(file = path, algo = "sha256")
}

software_versions <- function() {
  pk <- c("data.table", "yaml", "jsonlite", "Biostrings", "Rsamtools", "ggplot2", "httr", "digest")
  v <- lapply(pk, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NULL)
  names(v) <- pk
  tool <- function(cmd, args) {
    r <- tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE)), error = function(e) NULL)
    if (is.null(r)) NULL else r[1]
  }
  c(list(R = R.version.string), v,
    list(bcftools = tool("bcftools", "--version"), liftOver = if (nzchar(Sys.which("liftOver"))) "present" else NULL,
         bioh2m = tool("python", c("-c", shQuote("from importlib.metadata import version; print(version('bioh2m'))")))))
}

file_meta <- function(path, checksum = FALSE) {
  d <- list(path = path, bytes = if (file.exists(path)) file.size(path) else NULL)
  if (checksum) d$sha256 <- sha256(path)
  d
}

# merge this script's record into <out_dir>/provenance.json (patient inputs: path + size only)
write_provenance <- function(P, cfg, script, inputs = list(), resources = list(), outputs = character(), extra = list()) {
  f <- P$provenance
  doc <- if (file.exists(f)) fromJSON(f, simplifyVector = FALSE) else list()
  doc$ensembl_release <- cfg$ensembl_release
  doc$ensembl_release_tag <- cfg$ensembl_release_tag
  doc$seed <- cfg$seed
  doc$config <- list(path = cfg$config_path, sha256 = sha256(cfg$config_path))
  if (is.null(doc$scripts)) doc$scripts <- list()
  inputs <- Filter(Negate(is.null), inputs)
  resources <- Filter(Negate(is.null), resources)
  doc$scripts[[script]] <- c(list(
    run_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    software = software_versions(),
    inputs = lapply(inputs, file_meta),
    resources = lapply(resources, file_meta, checksum = TRUE),
    outputs = as.list(as.character(outputs))), extra)
  writeLines(toJSON(doc, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"), f)
}

write_tsv <- function(dt, path, log = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path, sep = "\t", na = "NA", quote = FALSE)
  if (!is.null(log)) log("wrote %s (%d rows)", path, nrow(dt))
  invisible(path)
}

# rows: data.table(category, n, counterpart) - counterpart yes / no / uncertain / info
write_summary <- function(P, layer, rows, log = NULL) {
  rows <- copy(rows)
  rows[, layer := layer]
  setcolorder(rows, "layer")
  write_tsv(rows, file.path(P$summary, paste0(layer, ".tsv")), log)
}

srow <- function(category, n, counterpart) data.table(category = category, n = as.integer(n), counterpart = counterpart)

find_one <- function(dir, pattern, label, log, required = TRUE) {
  hits <- if (dir.exists(dir)) sort(list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)) else character()
  if (!length(hits)) {
    msg <- sprintf("%s: no file matching %s under %s", label, pattern, dir)
    if (required) { log(msg, level = "ERROR"); stop(msg, call. = FALSE) }
    log("MISSING (optional) %s", msg, level = "WARNING")
    return(NULL)
  }
  if (length(hits) > 1) log("%s: %d matches, using %s", label, length(hits), hits[1], level = "WARNING")
  log("%-28s %s", label, hits[1])
  hits[1]
}

need_cols <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(label, ": missing columns ", paste(miss, collapse = ", "),
                         "; found ", paste(names(dt), collapse = ", "), call. = FALSE)
}

truthy <- function(x) tolower(as.character(x)) %chin% c("true", "1", "yes")

# ============================================================================
# chromosomes, GTF, orthologs
# ============================================================================

norm_chrom <- function(x) {
  x <- sub("^chr", "", as.character(x))
  fifelse(x %chin% c("MT", "M"), "chrM", paste0("chr", x))
}

chrom_order <- function(x) {
  k <- sub("^chr", "", x)
  n <- suppressWarnings(as.integer(k))
  fifelse(!is.na(n), n, fcase(k == "X", 100L, k == "Y", 101L, k == "M", 102L, default = 200L))
}

gtf_attr <- function(a, key) {
  r <- regexpr(paste0(key, ' "[^"]*"'), a)
  out <- rep(NA_character_, length(a))
  out[r > 0] <- sub(paste0(key, ' "([^"]*)"'), "\\1", regmatches(a, r))
  out
}

# gene and transcript tables from an Ensembl GTF (chr-prefixed, 1-based closed)
read_gtf_tables <- function(gtf) {
  cmd <- sprintf("gzip -dc %s | awk -F'\\t' '$3==\"gene\"||$3==\"transcript\"||$3==\"CDS\"'", shQuote(gtf))
  g <- fread(cmd = cmd, sep = "\t", header = FALSE, quote = "",
             col.names = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attr"),
             colClasses = list(character = c(1, 9)))
  g[, chrom := norm_chrom(chrom)]
  genes <- g[feature == "gene", .(gene_id = gtf_attr(attr, "gene_id"), gene_name = gtf_attr(attr, "gene_name"),
                                  chrom, start, end, strand, biotype = gtf_attr(attr, "gene_biotype"))]
  genes[is.na(gene_name), gene_name := gene_id]
  cds_tx <- unique(gtf_attr(g[feature == "CDS", attr], "transcript_id"))
  tx <- g[feature == "transcript", .(transcript_id = gtf_attr(attr, "transcript_id"), gene_id = gtf_attr(attr, "gene_id"),
                                     gene_name = gtf_attr(attr, "gene_name"), chrom, start, end, strand,
                                     biotype = gtf_attr(attr, "transcript_biotype"),
                                     canonical = grepl('tag "Ensembl_canonical"', attr, fixed = TRUE))]
  tx[, has_cds := transcript_id %chin% cds_tx]
  list(genes = genes, tx = tx)
}

load_genes <- function(f) fread(f, colClasses = list(character = c("chrom")))

strip_version <- function(x) sub("\\.[0-9]+$", "", x)

# symbol -> gene ID; ties: standard chromosome first, then lowest ID
symbol_to_gene_id <- function(symbols, genes) {
  g <- genes[gene_name %chin% unique(symbols)]
  g[, std := chrom %chin% STD_CHROMS]
  setorder(g, gene_name, -std, gene_id)
  g <- g[!duplicated(gene_name)]
  setNames(g$gene_id, g$gene_name)
}

# human gene ID per row: transcript ID (version stripped) first, then symbol
gene_ids_for <- function(symbol, genes, tx, transcript = NULL) {
  out <- rep(NA_character_, length(symbol))
  if (!is.null(transcript)) {
    t2g <- setNames(tx$gene_id, strip_version(tx$transcript_id))
    out <- unname(t2g[strip_version(as.character(transcript))])
  }
  by_sym <- symbol_to_gene_id(symbol[!is.na(symbol)], genes)
  fcoalesce(out, unname(by_sym[symbol]))
}

COMPARA_URL <- "https://ftp.ensembl.org/pub/release-%s/tsv/ensembl-compara/homologies/%s/Compara.%s.%s_default.homologies.tsv.gz"
COMPARA_SETS <- data.table(dir = rep(c("homo_sapiens", "mus_musculus"), each = 2), set = rep(c("protein", "ncrna"), 2))
compara_url <- function(release, dir, set) sprintf(COMPARA_URL, release, dir, release, set)
compara_files <- function(P, release) file.path(P$resources, paste0("ensembl", release),
                                                sprintf("Compara.%s.%s.%s_default.homologies.tsv.gz", release, COMPARA_SETS$dir, COMPARA_SETS$set))

# human -> mouse orthologs from the Ensembl Compara homology dumps (protein + ncRNA) of one release.
# Compara stores each human-mouse pair once, in either the homo_sapiens or the mus_musculus folder
# (release 110: 16,242 + 6,877 protein pairs, no overlap), so both are read; mouse-side rows are swapped.
# Every human gene of that release's GTF gets a row (gene_id_m NA when it has no mouse ortholog).
# identity = % of the human protein identical in the alignment (BioMart perc_id_r1), homology_identity
# = the same for the mouse protein (BioMart perc_id); is_high_confidence -> orthology_confidence.
# files: named by the species folder they come from ("homo_sapiens" / "mus_musculus").
orthologs_from_compara <- function(files, human_genes, mouse_genes, release) {
  read_side <- function(f, other) {
    cmd <- sprintf("gzip -dc %s | awk -F'\\t' 'NR==1 || ($8==\"%s\" && $5 ~ /^ortholog/)'", shQuote(f), other)
    fread(cmd = cmd, sep = "\t", quote = "", na.strings = "NULL", colClasses = "character")
  }
  hs <- rbindlist(lapply(files[names(files) == "homo_sapiens"], read_side, other = "mus_musculus"))
  ms <- rbindlist(lapply(files[names(files) == "mus_musculus"], read_side, other = "homo_sapiens"))
  h <- rbind(hs[, .(gene_id_h = gene_stable_id, gene_id_m = homology_gene_stable_id, homology_type,
                    perc_id_h_to_m = as.numeric(identity), perc_id_m_to_h = as.numeric(homology_identity),
                    orthology_confidence = is_high_confidence)],
             ms[, .(gene_id_h = homology_gene_stable_id, gene_id_m = gene_stable_id, homology_type,
                    perc_id_h_to_m = as.numeric(homology_identity), perc_id_m_to_h = as.numeric(identity),
                    orthology_confidence = is_high_confidence)])
  h <- unique(h, by = c("gene_id_h", "gene_id_m"))
  o <- h[human_genes[, .(gene_id_h = gene_id, symbol_h = gene_name, chrom_h = chrom, start_h = start, end_h = end)], on = "gene_id_h"]
  o[mouse_genes, on = c(gene_id_m = "gene_id"), `:=`(symbol_m = i.gene_name, chrom_m = i.chrom, start_m = i.start, end_m = i.end)]
  o[, ensembl_release := as.character(release)]
  setcolorder(o, c("ensembl_release", "gene_id_h", "symbol_h", "gene_id_m", "symbol_m", "homology_type", "perc_id_h_to_m",
                   "perc_id_m_to_h", "orthology_confidence", "chrom_h", "start_h", "end_h", "chrom_m", "start_m", "end_m"))
  o[order(gene_id_h, gene_id_m)]
}

read_orthologs <- function(f) {
  o <- fread(f, colClasses = "character", na.strings = c("", "NA"))
  o[]
}

# one row per human gene: homology (one2one / one2many / many2many / no_ortholog / not_in_release),
# n_mouse, mouse_targets (';'), gene_id_m (1:1 only)
classify_orthologs <- function(orth, human_ids = NULL) {
  o <- orth[!is.na(gene_id_h)]
  hit <- o[!is.na(gene_id_m)]
  agg <- hit[, .(n_mouse = uniqueN(gene_id_m),
                 homology = homology_type[which.min(match(homology_type, HOMOLOGY_ORDER, nomatch = 9L))],
                 mouse_targets = paste(sort(unique(gene_id_m)), collapse = ";")), by = gene_id_h]
  all_ids <- data.table(gene_id_h = sort(unique(c(o$gene_id_h, human_ids[!is.na(human_ids)]))))
  res <- agg[all_ids, on = "gene_id_h"]
  res[is.na(n_mouse), n_mouse := 0L]
  res[is.na(homology), homology := fifelse(gene_id_h %chin% o$gene_id_h, "no_ortholog", "not_in_release")]
  res[, gene_id_m := fifelse(homology == "ortholog_one2one" & n_mouse == 1L, mouse_targets, NA_character_)]
  setkey(res, gene_id_h)
  res[]
}

# long table: every human gene x each mouse target (or one row without target)
map_orthologs <- function(gene_ids_h, orth) {
  ids <- unique(gene_ids_h[!is.na(gene_ids_h)])
  cls <- classify_orthologs(orth, ids)
  q <- data.table(gene_id_h = ids)
  hits <- orth[!is.na(gene_id_m), !"symbol_h"]
  m <- hits[q, on = "gene_id_h"]
  m[cls, on = "gene_id_h", `:=`(homology_type = fcoalesce(homology_type, i.homology), n_mouse = i.n_mouse)]
  setorder(m, gene_id_h, gene_id_m, na.last = TRUE)
  m[]
}

ortholog_release_diff <- function(o1, o2, r1, r2) {
  c1 <- classify_orthologs(o1)[, .(gene_id_h, h1 = homology, t1 = mouse_targets)]
  c2 <- classify_orthologs(o2)[, .(gene_id_h, h2 = homology, t2 = mouse_targets)]
  d <- merge(c1, c2, by = "gene_id_h", all = TRUE)
  d <- d[fcoalesce(h1, "") != fcoalesce(h2, "") | fcoalesce(t1, "") != fcoalesce(t2, "")]
  d[, change := fcase(is.na(h1), "new_gene", is.na(h2), "gene_removed", h1 != h2, "homology_type_changed",
                      default = "mouse_target_changed")]
  setnames(d, c("h1", "t1", "h2", "t2"), c(paste0(c("homology_", "mouse_targets_"), r1), paste0(c("homology_", "mouse_targets_"), r2)))
  d[order(gene_id_h)]
}

# add <build>_chrom/_start/_end/_strand/_symbol for the mouse gene IDs in `id_col`
attach_mouse_coords <- function(dt, mouse_genes, build, id_col = "gene_id_m") {
  mg <- mouse_genes[, .(gene_id, chrom, start, end, strand, gene_name)]
  setnames(mg, c("gene_id", paste0(build, "_", c("chrom", "start", "end", "strand", "symbol"))))
  setnames(mg, "gene_id", id_col)
  drop <- intersect(names(dt), setdiff(names(mg), id_col))
  if (length(drop)) dt[, (drop) := NULL]
  mg[dt, on = id_col]
}

# ============================================================================
# copy number
# ============================================================================

GAIN_SIDE <- c("AMP", "GAIN")
LOSS_SIDE <- c("LOSS", "DEEP_DEL")

cn_state <- function(cn, ploidy, t) {
  r <- cn / ploidy
  fcase(is.na(cn), "not_assessable", cn < t$deep_del_abs, "DEEP_DEL", r >= t$amp_rel, "AMP",
        r >= t$gain_rel, "GAIN", r < t$loss_rel, "LOSS", default = "NEUTRAL")
}

# loss side from the gene's minimum CN, gain side from its maximum; partial = they disagree
gene_cn_state <- function(min_cn, max_cn, min_minor, ploidy, t) {
  s_min <- cn_state(min_cn, ploidy, t)
  s_max <- cn_state(max_cn, ploidy, t)
  state <- fcase(s_min == "not_assessable" & s_max == "not_assessable", "not_assessable",
                 s_min %chin% LOSS_SIDE, s_min, s_max %chin% GAIN_SIDE, s_max, default = "NEUTRAL")
  data.table(state = state,
             direction = fcase(state %chin% GAIN_SIDE, "gain", state %chin% LOSS_SIDE, "loss",
                               state == "not_assessable", "not_assessable", default = "neutral"),
             partial = state != "not_assessable" & s_min != s_max,
             loh = !is.na(min_minor) & min_minor < t$loh_minor_max & state != "DEEP_DEL" & state != "not_assessable")
}

# ============================================================================
# synteny blocks and intervals
# ============================================================================

# dt: 1:1 orthologs with gene_id_h, h_chrom, h_start, m_chrom, m_start. In mouse order a block
# continues while the human chromosome is the same and the human gene rank moves by <= max_rank_jump.
# Returns block IDs in the row order of dt.
synteny_blocks <- function(dt, max_rank_jump, prefix = "B") {
  d <- copy(dt)[, row := .I]
  setorder(d, h_chrom, h_start, gene_id_h)
  d[, h_rank := seq_len(.N), by = h_chrom]
  setorder(d, m_chrom, m_start, gene_id_h)
  d[, new := is.na(shift(m_chrom)) | m_chrom != shift(m_chrom) | h_chrom != shift(h_chrom) |
        abs(h_rank - shift(h_rank)) > max_rank_jump]
  d[, block := sprintf("%s%05d", prefix, cumsum(new))]
  d[order(row), block]
}

# Consecutive genes (chrom, start order) sharing an event label merge into intervals. A run breaks
# at another label (incl. neutral), a new chromosome, a gap > max_gap or (block_col) a new block.
# dt: chrom start end gene human_gene segment <label_col> [block_col]
merge_runs <- function(dt, label_col, events, max_gap, block_col = NULL) {
  empty <- data.table(chrom = character(), start = integer(), end = integer(), direction = character(),
                      genes = character(), human_genes = character(), segments = character(),
                      blocks = character(), n_genes = integer())
  if (!nrow(dt)) return(empty)
  d <- copy(dt)
  setorderv(d, c("chrom", "start", "end", "gene"))
  lab <- d[[label_col]]
  blk <- if (is.null(block_col)) rep(NA_character_, nrow(d)) else d[[block_col]]
  run <- integer(nrow(d))
  cur <- 0L; cend <- -Inf; clab <- NA_character_; cchr <- NA_character_; cblk <- NA_character_; open <- FALSE
  for (i in seq_len(nrow(d))) {
    if (!(lab[i] %chin% events)) { open <- FALSE; run[i] <- NA_integer_; next }
    same <- open && d$chrom[i] == cchr && lab[i] == clab && d$start[i] - cend <= max_gap &&
      (is.null(block_col) || identical(blk[i], cblk))
    if (!same) { cur <- cur + 1L; cchr <- d$chrom[i]; clab <- lab[i]; cblk <- blk[i]; cend <- d$end[i]; open <- TRUE }
    cend <- max(cend, d$end[i])
    run[i] <- cur
  }
  d[, run := run]
  d[, lab := lab][, blk := blk]
  d <- d[!is.na(run)]
  if (!nrow(d)) return(empty)
  d[, .(chrom = chrom[1], start = min(start), end = max(end), direction = lab[1],
        genes = paste(gene, collapse = ","), human_genes = paste(human_gene, collapse = ","),
        segments = paste(sort(unique(as.character(segment))), collapse = ","),
        blocks = paste(sort(unique(na.omit(blk))), collapse = ","), n_genes = .N), by = run][, !"run"]
}

# union of half-open intervals per chrom; intervals closer than `gap` merge
merge_overlapping <- function(dt, gap = 0) {
  if (!nrow(dt)) return(data.table(chrom = character(), start = numeric(), end = numeric()))
  d <- dt[, .(chrom, start, end)]
  setorder(d, chrom, start)
  d[, prev := shift(cummax(end)), by = chrom]
  d[, grp := cumsum(is.na(prev) | start > prev + gap), by = chrom]
  d[, .(start = min(start), end = max(end)), by = .(chrom, grp)][, !"grp"]
}

overlap_bp <- function(a, b) {
  if (!nrow(a) || !nrow(b)) return(0)
  A <- merge_overlapping(a)[, `:=`(s = start, e = end - 1)]
  B <- merge_overlapping(b)[, `:=`(s = start, e = end - 1)]
  setkey(B, chrom, s, e)
  o <- foverlaps(A, B, by.x = c("chrom", "s", "e"), nomatch = NULL)
  if (!nrow(o)) return(0)
  sum(pmax(0, pmin(o$end, o$i.end) - pmax(o$start, o$i.start)))
}

total_bp <- function(a) if (!nrow(a)) 0 else merge_overlapping(a)[, sum(end - start)]

parse_region <- function(s) {
  m <- regmatches(s, regexec("^(\\w+):([0-9,]+)-([0-9,]+)$", s))[[1]]
  if (!length(m)) stop("region must be chr:start-end, got ", s, call. = FALSE)
  list(chrom = norm_chrom(m[2]), start = as.numeric(gsub(",", "", m[3])), end = as.numeric(gsub(",", "", m[4])))
}

# ============================================================================
# SVs / fusions
# ============================================================================

# SV class that joins a 5' gene to a 3' gene given both genes' positions and strands
sv_class <- function(chrom5, pos5, strand5, chrom3, pos3, strand3) {
  up <- fifelse(strand5 == "+", pos5 < pos3, pos5 > pos3)
  fcase(is.na(chrom5) | is.na(pos5) | is.na(strand5) | is.na(chrom3) | is.na(pos3) | is.na(strand3), NA_character_,
        chrom5 != chrom3, "TRA", strand5 != strand3, "INV", up, "DEL", default = "DUP")
}

# a FLAG only - never a predicted mouse fusion
fusion_feasibility <- function(n5, n3, human_cls, mouse_cls) {
  fcase(n5 == 0 | n3 == 0, "not_possible_no_ortholog", is.na(mouse_cls), "not_assessable",
        !is.na(human_cls) & human_cls == mouse_cls, "feasible_same_sv_class",
        default = paste0("feasible_needs_", mouse_cls))
}

# ============================================================================
# SNVs
# ============================================================================

PAVE_DEFAULT <- c("gene", "transcript", "canonical_effect", "canonical_coding_effect", "splice_region",
                  "hgvs_coding", "hgvs_protein", "other_effects", "worst_coding_effect", "genes_affected")
PAVE_ALIASES <- c(gene = "gene", transcript = "transcript", canonicaleffect = "canonical_effect",
                  canonicalcodingeffect = "canonical_coding_effect", spliceregion = "splice_region",
                  hgvscodingimpact = "hgvs_coding", hgvsproteinimpact = "hgvs_protein",
                  otherreportableeffects = "other_effects", worstcodingeffect = "worst_coding_effect",
                  genesaffected = "genes_affected")

pave_fields <- function(header_line) {
  m <- regmatches(header_line, regexpr("\\[[^]]*\\]", header_line))
  if (!length(m)) return(PAVE_DEFAULT)
  n <- gsub("[^a-z]", "", tolower(strsplit(gsub("^\\[|\\]$", "", m), ",")[[1]]))
  fcoalesce(unname(PAVE_ALIASES[n]), n)
}

# PASS SNVs/indels on chr1-22/X/Y from PURPLE's somatic VCF with PAVE's canonical impact (bcftools)
read_purple_somatic <- function(vcf, log) {
  hdr <- system2("bcftools", c("view", "-h", shQuote(vcf)), stdout = TRUE)
  info_ids <- sub("^##INFO=<ID=([^,]+),.*", "\\1", grep("^##INFO=<ID=", hdr, value = TRUE))
  imp <- grep("^##INFO=<ID=IMPACT,", hdr, value = TRUE)
  fields <- if (length(imp)) pave_fields(imp) else { log("no INFO/IMPACT in %s", vcf, level = "WARNING"); PAVE_DEFAULT }
  log("PAVE IMPACT fields: %s", paste(fields, collapse = ","))
  tags <- intersect(c("IMPACT", "REPORTED", "TIER", "HOTSPOT", "PURPLE_AF", "PURPLE_CN", "BIALLELIC", "SUBCL"), info_ids)
  fmt <- paste0("%CHROM\\t%POS\\t%REF\\t%ALT", paste0("\\t%INFO/", tags, collapse = ""), "\\n")
  cmd <- sprintf("bcftools query -i 'FILTER=\"PASS\"' -f '%s' %s", fmt, shQuote(vcf))
  v <- fread(cmd = cmd, sep = "\t", header = FALSE, quote = "", na.strings = ".", colClasses = "character",
             col.names = c("chrom", "pos", "ref", "alt", tags))
  v <- v[chrom %chin% STD_CHROMS]
  v[, pos := as.integer(pos)]
  if (any(grepl(",", v$alt))) log("%d records with several ALTs, using the first", sum(grepl(",", v$alt)), level = "WARNING")
  v[, alt := sub(",.*", "", alt)]
  for (t in setdiff(c("IMPACT", "REPORTED", "TIER", "HOTSPOT", "PURPLE_AF", "PURPLE_CN", "BIALLELIC", "SUBCL"), tags))
    set(v, j = t, value = NA_character_)
  parts <- tstrsplit(fcoalesce(v$IMPACT, ""), ",", fixed = TRUE, fill = NA_character_)
  get_f <- function(f) { k <- match(f, fields); if (is.na(k) || k > length(parts)) NA_character_ else parts[[k]] }
  v[, `:=`(gene_symbol_h = get_f("gene"), transcript_h = get_f("transcript"),
           pave_effect = get_f("canonical_effect"), pave_coding_effect = get_f("canonical_coding_effect"),
           splice_region = truthy(get_f("splice_region")),
           hgvsc_h = get_f("hgvs_coding"), hgvsp_h = get_f("hgvs_protein"),
           reported = truthy(REPORTED), tier = TIER, hotspot = TIER %chin% "HOTSPOT" | truthy(HOTSPOT),
           purple_af = as.numeric(PURPLE_AF), purple_cn = as.numeric(PURPLE_CN),
           biallelic = truthy(BIALLELIC), subclonal = as.numeric(SUBCL))]
  for (c in c("gene_symbol_h", "transcript_h", "hgvsc_h", "hgvsp_h")) set(v, i = which(v[[c]] == ""), j = c, value = NA_character_)
  v[, c("IMPACT", "REPORTED", "TIER", "HOTSPOT", "PURPLE_AF", "PURPLE_CN", "BIALLELIC", "SUBCL") := NULL]
  v[, variant_id := paste(chrom, pos, ref, alt, sep = ":")]
  v[, vtype := fcase(nchar(ref) == 1 & nchar(alt) == 1, "SNV", nchar(ref) == nchar(alt), "MNV", default = "INDEL")]
  log("PURPLE somatic PASS: %d records", nrow(v))
  v[]
}

# VCF-style (anchored indels) -> H2M/MAF style, as bioh2m's vcf_to_maf: indels drop the anchor
vcf_to_maf <- function(pos, ref, alt) {
  lr <- nchar(ref); la <- nchar(alt)
  snp <- lr == la
  ins1 <- !snp & lr == 1 & la >= 2
  del1 <- !snp & lr >= 2 & la == 1
  cplx <- !snp & !ins1 & !del1
  data.table(
    start_h = fcase(snp | ins1, pos, default = pos + 1L),
    end_h = fcase(snp, pos + lr - 1L, ins1, pos + 1L, del1, pos + lr - 1L, default = pos + lr - 1L),
    ref_seq_h = fcase(snp, ref, ins1, "-", default = substring(ref, 2)),
    alt_seq_h = fcase(snp, alt, del1, "-", default = substring(alt, 2)),
    type_h = fcase(snp & lr == 1, "SNP", snp & lr == 2, "DNP", snp & lr == 3, "TNP", snp, "ONP",
                   ins1 | (cplx & la > lr), "INS", default = "DEL"))
}

is_frameshift <- function(ref, alt) (nchar(ref) - nchar(alt)) %% 3 != 0

AA3 <- c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q", Glu = "E", Gly = "G", His = "H",
         Ile = "I", Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P", Ser = "S", Thr = "T", Trp = "W",
         Tyr = "Y", Val = "V", Ter = "*", Sec = "U")

# 'p.Arg172His' / 'p.R172H' / 'R172H' -> 'R172H'
short_protein <- function(x) {
  s <- gsub("^\\(|\\)$", "", sub("^.*p\\.", "", as.character(x)))
  for (k in names(AA3)) s <- gsub(k, AA3[[k]], s, fixed = TRUE)
  s[is.na(x) | x == ""] <- NA_character_
  s
}

same_protein_change <- function(observed, expected) {
  o <- short_protein(observed)
  !is.na(o) & o == short_protein(expected)
}

# ============================================================================
# trinucleotide contexts / spectrum rescaling
# ============================================================================

BASES <- c("A", "C", "G", "T")
SUBS <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
SBS96 <- unlist(lapply(SUBS, function(s) as.vector(outer(BASES, BASES, function(a, b) paste0(a, "[", s, "]", b)))))
SBS96 <- SBS96[order(match(substr(SBS96, 3, 5), SUBS), substr(SBS96, 1, 1), substr(SBS96, 7, 7))]
TRI32 <- unlist(lapply(c("C", "T"), function(m) as.vector(t(outer(BASES, BASES, function(a, b) paste0(a, m, b))))))
TRI32 <- TRI32[order(substr(TRI32, 2, 2), substr(TRI32, 1, 1), substr(TRI32, 3, 3))]

revcomp <- function(s) vapply(strsplit(chartr("ACGT", "TGCA", s), ""), function(x) paste(rev(x), collapse = ""), "")

# 64 trinucleotide counts (named) -> 32 pyrimidine-centred
collapse_trinuc <- function(c64) {
  k <- names(c64)
  key <- fifelse(substr(k, 2, 2) %chin% c("C", "T"), k, revcomp(k))
  r <- tapply(as.numeric(c64), key, sum)
  out <- setNames(as.numeric(r[TRI32]), TRI32)
  out[is.na(out)] <- 0
  out
}

# ctx = 5'-ref-3' on the + strand -> 'A[C>T]G' on the pyrimidine strand
sbs96_label <- function(ctx, ref, alt) {
  flip <- ref %chin% c("A", "G")
  ctx[flip] <- revcomp(ctx[flip])
  ref[flip] <- chartr("ACGT", "TGCA", ref[flip])
  alt[flip] <- chartr("ACGT", "TGCA", alt[flip])
  paste0(substr(ctx, 1, 1), "[", ref, ">", alt, "]", substr(ctx, 3, 3))
}

# human SBS96 counts -> expected counts in mouse opportunity: channel x (mouse / human trinucleotide
# frequency), scaled back to the human total
rescale_spectrum <- function(counts96, tri_h, tri_m) {
  c96 <- setNames(as.numeric(counts96[SBS96]), SBS96)
  c96[is.na(c96)] <- 0
  fh <- tri_h / sum(tri_h); fm <- tri_m / sum(tri_m)
  tri <- paste0(substr(SBS96, 1, 1), substr(SBS96, 3, 3), substr(SBS96, 7, 7))
  ratio <- unname(fm[tri] / fh[tri])
  raw <- c96 * ratio
  scaled <- if (sum(raw) > 0) raw * sum(c96) / sum(raw) else raw
  data.table(MutationType = SBS96, trinucleotide = tri, human_count = unname(c96),
             human_prop = if (sum(c96)) unname(c96 / sum(c96)) else 0,
             human_tri_freq = unname(fh[tri]), mouse_tri_freq = unname(fm[tri]), ratio_mouse_human = ratio,
             expected_mouse_count = unname(scaled),
             expected_mouse_prop = if (sum(scaled)) unname(scaled / sum(scaled)) else 0)
}

# 32 pyrimidine-centred trinucleotide counts of the autosomes of a FASTA (Biostrings counts
# only A/C/G/T triplets), cached as TSV
genome_trinuc <- function(fasta, n_autosomes, cache, log) {
  if (file.exists(cache)) { d <- fread(cache); return(setNames(d$count, d$trinucleotide)) }
  fa <- Rsamtools::FaFile(fasta)
  idx <- Rsamtools::seqinfo(fa)
  tot <- NULL
  for (i in seq_len(n_autosomes)) {
    nm <- intersect(c(paste0("chr", i), as.character(i)), GenomeInfoDb::seqnames(idx))[1]
    if (is.na(nm)) stop(fasta, ": chromosome ", i, " missing", call. = FALSE)
    gr <- GenomicRanges::GRanges(nm, IRanges::IRanges(1, GenomeInfoDb::seqlengths(idx)[[nm]]))
    s <- Rsamtools::scanFa(fa, gr)[[1]]
    f <- Biostrings::trinucleotideFrequency(s)
    tot <- if (is.null(tot)) f else tot + f
    log("trinucleotides %s %s", basename(fasta), nm)
  }
  t32 <- collapse_trinuc(tot)
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(trinucleotide = names(t32), count = t32), cache, sep = "\t")
  t32
}

# ============================================================================
# liftOver (UCSC binary)
# ============================================================================

# bed: chrom start(0-based) end name [strand]; names unique.
# Returns list(mapped = chrom/start/end/name/strand, status = name/status/n_hits/reason)
run_liftover <- function(bed, chain, min_match, workdir, tag, multiple = FALSE) {
  if (!nzchar(Sys.which("liftOver"))) stop("liftOver not on PATH - activate the cross-species-h2m env", call. = FALSE)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  inp <- file.path(workdir, paste0(tag, ".in.bed"))
  out <- file.path(workdir, paste0(tag, ".out.bed"))
  un <- file.path(workdir, paste0(tag, ".unmapped.bed"))
  b <- bed[, .(chrom, start = as.integer(start), end = as.integer(end), name = as.character(name), score = 0L,
               strand = if ("strand" %chin% names(bed)) strand else "+")]
  fwrite(b, inp, sep = "\t", col.names = FALSE, quote = FALSE)
  args <- c(paste0("-minMatch=", min_match), if (multiple) "-multiple", inp, chain, out, un)
  rc <- system2("liftOver", args, stdout = FALSE, stderr = FALSE)
  if (rc != 0) stop("liftOver failed (", tag, ")", call. = FALSE)
  mapped <- if (file.size(out) > 0) fread(out, header = FALSE, select = 1:6, colClasses = list(character = 4),
                                         col.names = c("chrom", "start", "end", "name", "score", "strand")) else
    data.table(chrom = character(), start = integer(), end = integer(), name = character(), score = integer(), strand = character())
  mapped[, score := NULL]
  ul <- readLines(un)
  reason <- character(); last <- NA_character_
  for (l in ul) {
    if (startsWith(l, "#")) last <- sub("^#", "", l)
    else if (nzchar(l)) reason[strsplit(l, "\t", fixed = TRUE)[[1]][4]] <- last
  }
  n <- mapped[, .N, by = name]
  st <- data.table(name = b$name)
  st[n, on = "name", n_hits := i.N]
  st[is.na(n_hits), n_hits := 0L]
  st[, status := fcase(n_hits == 1L, "mapped", n_hits > 1L, "multi_mapped", default = "unmapped")]
  st[, reason := unname(reason[name])]
  list(mapped = mapped, status = st)
}

# ============================================================================
# gene-level events from the per-layer outputs (05b, 06, 06b)
# ============================================================================

NONCODING_GENE_CLASSES <- c("TERT_promoter", "promoter", "UTR5", "UTR3", "splice_region")
INACTIVATING <- c("NONSENSE_OR_FRAMESHIFT", "SPLICE", "DEEP_DEL", "LOSS", "LOH", "HOM_DISRUPTION", "DISRUPTION")
ACTIVATING <- c("AMP", "GAIN")
ROLE_DEPENDENT <- c("MISSENSE", "TERT_promoter", "promoter")
MATRIX_CODE <- c(AMP = "AMP", GAIN = "GAIN", LOSS = "DEL", DEEP_DEL = "DEL", LOH = "LOH",
                 HOM_DISRUPTION = "SV", DISRUPTION = "SV", FUSION_PARTNER = "FUSION")

read_layer <- function(P, name) {
  f <- file.path(P$hp, name)
  if (file.exists(f)) fread(f, na.strings = "NA") else NULL
}

# long table gene_id_h symbol_h layer event matrix_code confidence detail event_key (one row per human
# event x gene). confidence of the human call: high = PURPLE-reported / hotspot / driver-gene SNV,
# whole-gene CN, LOH, reported SV, TERT promoter; medium = other coding SNV, partial-gene CN,
# unreported disruption; low = frameshift / splice / H2M class 2, other non-coding.
collect_gene_events <- function(P) {
  out <- list()
  snv <- read_layer(P, "expected_mouse_snv.tsv")
  if (!is.null(snv) && nrow(snv)) {
    s <- unique(snv, by = "variant_id")
    out$snv <- s[, .(gene_id_h, symbol_h = gene_symbol_h, layer = "snv_coding", event = pave_coding_effect,
                     matrix_code = fifelse(vtype == "INDEL", "INDEL", "SNV"),
                     confidence = fifelse(truthy(uncertain_effect), "low",
                                          fifelse(truthy(reported) | truthy(hotspot) | truthy(purple_driver_gene), "high", "medium")),
                     detail = fcoalesce(hgvsp_h, hgvsc_h, ""), event_key = variant_id)]
  }
  nc <- read_layer(P, "expected_mouse_noncoding.tsv")
  if (!is.null(nc) && nrow(nc))
    out$nc <- nc[class %chin% NONCODING_GENE_CLASSES & !is.na(gene_id_h),
                 .(gene_id_h, symbol_h = gene_symbol_h, layer = "snv_noncoding", event = class,
                   matrix_code = fifelse(vtype == "INDEL", "INDEL", "SNV"),
                   confidence = fifelse(class == "TERT_promoter", "high", "low"), detail = variant_id, event_key = variant_id)]
  cn <- read_layer(P, "expected_mouse_cn_genes.tsv")
  if (!is.null(cn) && nrow(cn)) {
    out$cn <- cn[state %chin% c(GAIN_SIDE, LOSS_SIDE),
                 .(gene_id_h, symbol_h = gene, layer = "cn", event = state, matrix_code = unname(MATRIX_CODE[state]),
                   confidence = fifelse(truthy(partial), "medium", "high"),
                   detail = sprintf("CN %s-%s; segment %s", minCopyNumber, maxCopyNumber, segment),
                   event_key = paste("CN", segment, state, sep = ":"))]
    out$loh <- cn[truthy(loh), .(gene_id_h, symbol_h = gene, layer = "loh", event = "LOH", matrix_code = "LOH", confidence = "high",
                                 detail = paste("minor CN", minMinorAlleleCopyNumber), event_key = paste("LOH", segment, sep = ":"))]
  }
  sv <- read_layer(P, "expected_mouse_sv_genes.tsv")
  if (!is.null(sv) && nrow(sv))
    out$sv <- sv[event != "NON_DISRUPTIVE", .(gene_id_h, symbol_h = gene, layer = "sv", event,
                                              matrix_code = unname(MATRIX_CODE[event]),
                                              confidence = fifelse(truthy(reported), "high", "medium"),
                                              detail = paste("SV", sv_ids), event_key = paste("SV", sv_ids, sep = ":"))]
  ev <- rbindlist(out, use.names = TRUE)
  if (!nrow(ev)) ev <- data.table(gene_id_h = character(), symbol_h = character(), layer = character(), event = character(),
                                  matrix_code = character(), confidence = character(), detail = character(), event_key = character())
  ev[]
}

# ============================================================================
# frozen predictions (07b writes them; 10-12 read only from there)
# ============================================================================

FROZEN_MANIFEST <- "MANIFEST.sha256"

# stop unless every file in FROZEN/MANIFEST.sha256 exists with the same checksum and nothing was
# added; returns the FROZEN directory
verify_frozen <- function(P, log = NULL) {
  d <- P$frozen
  mf <- file.path(d, FROZEN_MANIFEST)
  if (!file.exists(mf)) stop("no frozen predictions (", mf, ") - run 07b_freeze.R first", call. = FALSE)
  x <- readLines(mf)
  m <- data.table(sha256 = sub(" .*", "", x), file = sub("^[0-9a-f]+ [ *]", "", x))
  now <- vapply(file.path(d, m$file), function(f) if (file.exists(f)) digest::digest(file = f, algo = "sha256") else NA_character_, "")
  bad <- m$file[is.na(now) | now != m$sha256]
  extra <- setdiff(list.files(d, recursive = TRUE), c(m$file, FROZEN_MANIFEST, "FROZEN_AT.txt"))
  if (length(bad) || length(extra))
    stop("FROZEN predictions changed - refusing to run. Changed/missing: ", paste(bad, collapse = ", "),
         if (length(extra)) paste0("; not in manifest: ", paste(extra, collapse = ", ")), call. = FALSE)
  if (!is.null(log)) log("FROZEN predictions verified: %d files (%s)", nrow(m), readLines(file.path(d, "FROZEN_AT.txt"), n = 1))
  d
}

CONFIDENCE_ORDER <- c("high", "medium", "low")
best_confidence <- function(x) CONFIDENCE_ORDER[min(match(x, CONFIDENCE_ORDER, nomatch = 3L))]

cgc_roles <- function(role) if (is.null(role) || is.na(role)) character() else trimws(strsplit(role, ",", fixed = TRUE)[[1]])

# expected functional direction of one event given the gene's CGC role string
event_direction <- function(event, role) {
  if (event %chin% INACTIVATING) return("inactivating")
  if (event %chin% ACTIVATING) return("activating")
  if (event %chin% ROLE_DEPENDENT) {
    r <- cgc_roles(role); onc <- "oncogene" %chin% r; tsg <- "TSG" %chin% r
    if (onc && !tsg) return("activating")
    if (tsg && !onc) return("inactivating")
  }
  "unknown"
}

role_consistency <- function(direction, role) {
  r <- cgc_roles(role); onc <- "oncogene" %chin% r; tsg <- "TSG" %chin% r
  if (!length(r) || direction %chin% c("unknown", "mixed")) return("")
  if ((direction == "inactivating" && tsg) || (direction == "activating" && onc)) return("consistent")
  if ((direction == "inactivating" && onc && !tsg) || (direction == "activating" && tsg && !onc)) return("discordant")
  ""
}

combine_directions <- function(dirs) {
  d <- unique(dirs[dirs != "unknown"])
  if (!length(d)) "unknown" else if (length(d) == 1) d else "mixed"
}
