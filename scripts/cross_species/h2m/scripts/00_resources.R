# ---------------------------------------------------------------------------
# 00_resources.R - download and pin the references for one Ensembl release.
#
#   Rscript scripts/00_resources.R                     # release from config (110)
#   Rscript scripts/00_resources.R --release latest    # current Ensembl, for a second run
#   Rscript scripts/00_resources.R --dry-run           # list files and sizes only
#
# Into <resources_dir>/ (files already there are skipped):
#   ensembl<R>/orthologs_hsapiens_mmusculus.tsv  human -> mouse orthologs from the release's Ensembl Compara
#                                                homology dumps, homo_sapiens + mus_musculus folders (each pair
#                                                is stored in one of them), protein + ncRNA, kept: gene IDs, symbols, homology type, %id both ways,
#                                                high-confidence flag, coords from the release GTFs (GRCm39).
#                                                Not BioMart: the archive BioMart returned HTTP 500 (2026-10-28).
#   ensembl<R>/Homo_sapiens.GRCh38.<R>.gtf.gz    + human_genes.tsv.gz, human_transcripts.tsv.gz
#   ensembl<R>/Mus_musculus.GRCm39.<R>.gtf.gz    + mouse_GRCm39_{genes,transcripts}.tsv.gz
#   ensembl102/Mus_musculus.GRCm38.102.gtf.gz    + mouse_GRCm38_{genes,transcripts}.tsv.gz
#   fasta/Mus_musculus.<build>.dna.primary_assembly.fa (+ .fai)   H2M (GRCm39) and 05 contexts
#   chains/hg38ToMm39, hg38ToMm10, mm39ToMm10 .over.chain.gz      UCSC liftOver
#   reactome/v<N>/Ensembl2Reactome.txt, ReactomePathways.txt, ReactomePathwaysRelation.txt
#                                                current Reactome release N (05b); pin with reactome.version
#   versions.json                                URL, bytes, sha256, download date per file
# With ensembl_release_check set, also that release's Compara files, human / GRCm39 GTF tables,
# ortholog table and
#   <out_dir>/human_predictions/ortholog_release_diff_<R>_vs_<R2>.tsv
#
# Stops before downloading more than max_download_gb (config) unless --allow-large.
# The human GRCh38 FASTA is not downloaded: config human_fasta (iGenomes, already on the server).
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

ENSEMBL_FTP <- "https://ftp.ensembl.org/pub/release-%s"
UCSC <- "https://hgdownload.soe.ucsc.edu/goldenPath/%s/liftOver/%s"
REACTOME_BASE <- "https://reactome.org/download/current"
REACTOME_VERSION <- "https://reactome.org/ContentService/data/database/version"
REACTOME_FILES <- c("Ensembl2Reactome.txt", "ReactomePathways.txt", "ReactomePathwaysRelation.txt")

remote_size <- function(url) {
  r <- tryCatch(httr::HEAD(url, httr::timeout(60)), error = function(e) NULL)
  if (is.null(r)) return(NA_real_)
  as.numeric(httr::headers(r)[["content-length"]] %||% NA)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

download <- function(url, dest, log) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(dest, ".part")
  log("downloading %s", url)
  old <- options(timeout = 36000); on.exit(options(old))
  if (download.file(url, tmp, mode = "wb", quiet = TRUE) != 0) stop("download failed: ", url, call. = FALSE)
  file.rename(tmp, dest)
}

main <- function() {
  a <- parse_cli(c("allow-large", "dry-run"))
  cfg <- load_config(a$config, a$release, query_latest = TRUE)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("00_resources", P)
  R <- cfg$ensembl_release
  log("Ensembl release %s (%s)", R, cfg$ensembl_release_tag)

  # ---- files to fetch ------------------------------------------------------------------
  files <- list(list(url = sprintf("%s/gtf/homo_sapiens/Homo_sapiens.GRCh38.%s.gtf.gz", sprintf(ENSEMBL_FTP, R), R),
                     dest = P$human_gtf))
  for (build in names(P$mouse)) {
    M <- P$mouse[[build]]; gr <- M$gtf_release; ftp <- sprintf(ENSEMBL_FTP, gr)
    files <- c(files, list(
      list(url = sprintf("%s/gtf/mus_musculus/Mus_musculus.%s.%s.gtf.gz", ftp, build, gr), dest = M$gtf),
      list(url = sprintf("%s/fasta/mus_musculus/dna/Mus_musculus.%s.dna.primary_assembly.fa.gz", ftp, build),
           dest = paste0(M$fasta, ".gz"), final = M$fasta),
      list(url = sprintf(UCSC, "hg38", basename(M$chain)), dest = M$chain)))
  }
  files <- c(files, list(list(url = sprintf(UCSC, "mm39", basename(P$mm39_to_mm10)), dest = P$mm39_to_mm10)))
  cf <- compara_files(P, R)
  files <- c(files, lapply(seq_along(cf), function(i) list(url = compara_url(R, COMPARA_SETS$dir[i], COMPARA_SETS$set[i]), dest = cf[i])))
  todo <- list()
  for (f in files) {
    done <- file.exists(f$dest) || (!is.null(f$final) && file.exists(f$final))
    size <- if (done) NA_real_ else remote_size(f$url)
    log("%-8s %8s  %s", if (done) "present" else "fetch", if (done) "" else sprintf("%.2f GB", fcoalesce(size, 0) / 1e9), basename(f$dest))
    if (!done) todo <- c(todo, list(c(f, size = size)))
  }
  total_gb <- sum(vapply(todo, function(f) fcoalesce(f$size, 0), 0)) / 1e9
  log("to download: %d files, %.2f GB", length(todo), total_gb)
  if (isTRUE(a$dry_run)) return(invisible())
  if (total_gb > cfg$max_download_gb && !isTRUE(a$allow_large)) {
    log("%.2f GB > max_download_gb %s - rerun with --allow-large if that is OK", total_gb, cfg$max_download_gb, level = "ERROR")
    quit(status = 2)
  }

  vf <- file.path(P$resources, "versions.json")
  versions <- if (file.exists(vf)) fromJSON(vf, simplifyVector = FALSE) else list()
  stamp <- function(dest, url, ...) versions[[sub(paste0("^", P$resources, "/"), "", dest)]] <<-
    c(list(url = url, bytes = file.size(dest), sha256 = sha256(dest, limit_gb = 5),
           downloaded = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), list(...))
  for (f in todo) { download(f$url, f$dest, log); stamp(f$dest, f$url) }

  # ---- FASTA: plain .fa + .fai (Rsamtools / pysam cannot index plain gzip) ---------------------
  for (build in names(P$mouse)) {
    M <- P$mouse[[build]]; gz <- paste0(M$fasta, ".gz")
    if (!file.exists(M$fasta) && file.exists(gz)) {
      log("decompressing %s", basename(gz))
      if (system2("gzip", c("-d", shQuote(gz))) != 0) stop("gzip -d failed: ", gz, call. = FALSE)
    }
    if (!file.exists(paste0(M$fasta, ".fai"))) Rsamtools::indexFa(M$fasta)
  }

  # ---- GTF -> gene / transcript tables ------------------------------------------------------
  tabs <- c(list(list(gtf = P$human_gtf, genes = P$human_genes, tx = P$human_tx)),
            lapply(P$mouse, function(M) list(gtf = M$gtf, genes = M$genes, tx = M$tx)))
  for (t in tabs) {
    if (file.exists(t$genes) && file.exists(t$tx)) next
    log("parsing %s", basename(t$gtf))
    g <- read_gtf_tables(t$gtf)
    write_tsv(g$genes, t$genes, log)
    write_tsv(g$tx, t$tx, log)
  }

  # ---- orthologs (Ensembl Compara homology dumps of the release) ------------------------------------
  fetch <- function(url, dest) if (!file.exists(dest)) { download(url, dest, log); stamp(dest, url) }
  tables <- function(gtf, gf, tf) if (!(file.exists(gf) && file.exists(tf))) {
    log("parsing %s", basename(gtf)); g <- read_gtf_tables(gtf); write_tsv(g$genes, gf, log); write_tsv(g$tx, tf, log)
  }
  orth_for <- function(release) {
    rd <- file.path(P$resources, paste0("ensembl", release))
    dest <- file.path(rd, "orthologs_hsapiens_mmusculus.tsv")
    if (file.exists(dest)) return(read_orthologs(dest))
    ftp <- sprintf(ENSEMBL_FTP, release)
    hgtf <- file.path(rd, sprintf("Homo_sapiens.GRCh38.%s.gtf.gz", release))
    mgtf <- file.path(rd, sprintf("Mus_musculus.GRCm39.%s.gtf.gz", release))
    fetch(sprintf("%s/gtf/homo_sapiens/Homo_sapiens.GRCh38.%s.gtf.gz", ftp, release), hgtf)
    fetch(sprintf("%s/gtf/mus_musculus/Mus_musculus.GRCm39.%s.gtf.gz", ftp, release), mgtf)
    tables(hgtf, file.path(rd, "human_genes.tsv.gz"), file.path(rd, "human_transcripts.tsv.gz"))
    tables(mgtf, file.path(rd, "mouse_GRCm39_genes.tsv.gz"), file.path(rd, "mouse_GRCm39_transcripts.tsv.gz"))
    cf <- compara_files(P, release)
    for (i in seq_along(cf)) fetch(compara_url(release, COMPARA_SETS$dir[i], COMPARA_SETS$set[i]), cf[i])
    log("building human -> mouse ortholog table from %s", paste(basename(cf), collapse = ", "))
    o <- orthologs_from_compara(setNames(cf, COMPARA_SETS$dir), load_genes(file.path(rd, "human_genes.tsv.gz")),
                                load_genes(file.path(rd, "mouse_GRCm39_genes.tsv.gz")), release)
    write_tsv(o, dest, log)
    read_orthologs(dest)
  }
  orth <- orth_for(R)
  cls <- classify_orthologs(orth)
  log("orthologs (release %s): %s", R, paste(cls[, .N, by = homology][, paste0(homology, "=", N)], collapse = ", "))

  outputs <- P$orthologs
  check <- cfg$ensembl_release_check
  if (!is.null(check)) {
    R2 <- resolve_release(check, cfg$resources_dir, query = TRUE)
    if (R2 != R) {
      o2 <- orth_for(R2)
      d <- ortholog_release_diff(orth, o2, R, R2)
      f <- write_tsv(d, file.path(P$hp, sprintf("ortholog_release_diff_%s_vs_%s.tsv", R, R2)), log)
      log("orthologs differing between %s and %s: %s", R, R2, paste(d[, .N, by = change][, paste0(change, "=", N)], collapse = ", "))
      outputs <- c(outputs, f)
    } else log("ensembl_release_check %s resolves to the same release - no diff", check)
  }

  # ---- Reactome (current release, kept in a versioned folder) ---------------------------------------
  want <- cfg$reactome$version
  rdir <- if (!is.null(want)) reactome_dir(P, want) else NULL
  if (is.null(rdir) || !file.exists(file.path(rdir, "Ensembl2Reactome.txt"))) {
    cur <- as.integer(trimws(httr::content(httr::GET(REACTOME_VERSION, httr::timeout(60)), as = "text", encoding = "UTF-8")))
    if (!is.null(want) && as.integer(want) != cur)
      stop("config reactome.version ", want, " is not on disk and Reactome now serves ", cur, call. = FALSE)
    rdir <- reactome_dir(P, cur)
    for (fn in REACTOME_FILES) {
      dest <- file.path(rdir, fn)
      if (!file.exists(dest)) { download(file.path(REACTOME_BASE, fn), dest, log); stamp(dest, file.path(REACTOME_BASE, fn), reactome_release = cur) }
    }
    log("Reactome release %d in %s", cur, rdir)
  }

  writeLines(toJSON(versions, auto_unbox = TRUE, pretty = TRUE), vf)
  res <- list(orthologs = P$orthologs, human_gtf = P$human_gtf, mm39_to_mm10 = P$mm39_to_mm10)
  for (b in names(P$mouse)) res[paste0(b, c("_gtf", "_chain", "_fasta_fai"))] <- list(P$mouse[[b]]$gtf, P$mouse[[b]]$chain, paste0(P$mouse[[b]]$fasta, ".fai"))
  write_provenance(P, cfg, "00_resources", resources = res, outputs = outputs,
                   extra = list(resource_versions = versions, human_fasta = cfg$human_fasta))
  log("done")
}

main()
