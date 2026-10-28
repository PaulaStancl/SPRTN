# ---------------------------------------------------------------------------
# 07b_freeze.R - fix the predictions before any mouse data are looked at.
#
#   Rscript scripts/07b_freeze.R [--config ...] [--release ...]           # freeze
#   Rscript scripts/07b_freeze.R --verify                                  # check an existing freeze
#
# Copies every prediction file of human_predictions/ (tables, BED/BEDPE, summary/, report.html)
# plus provenance.json and the config into human_predictions/FROZEN/, writes MANIFEST.sha256
# (sha256sum format, paths relative to FROZEN/) and FROZEN_AT.txt (date, Ensembl release, git
# commit of these scripts), and makes the copies read-only. An existing FROZEN/ is never
# overwritten: move it away yourself if the predictions really have to change (and say why).
# The mouse modules (10-12) read predictions only from FROZEN/ and stop if verify_frozen()
# (common.R) finds a changed, missing or extra file.
#
# Output  human_predictions/FROZEN/{<prediction files>, provenance.json, config.yaml,
#                                   MANIFEST.sha256, FROZEN_AT.txt}
# ---------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

main <- function() {
  a <- parse_cli("verify")
  cfg <- load_config(a$config, a$release)
  P <- paths(cfg)
  make_dirs(P)
  log <- make_logger("07b_freeze", P)
  if (isTRUE(a$verify)) { verify_frozen(P, log); return(invisible()) }
  if (dir.exists(P$frozen)) stop(P$frozen, " exists - predictions are already frozen. Move it away first if they must change.", call. = FALSE)

  files <- list.files(P$hp, recursive = TRUE, full.names = FALSE)
  files <- files[!startsWith(files, "FROZEN/") & grepl("\\.(tsv|bed|bedpe|html)$", files)]
  need <- c("expected_mouse_snv.tsv", "expected_mouse_noncoding.tsv", "expected_mouse_cn_genes.tsv", "expected_mouse_sv_genes.tsv",
            "expected_mouse_sbs96.tsv", "gene_pathways.tsv", "gene_matrix.tsv", "event_table.tsv", "report.html")
  miss <- setdiff(need, files)
  if (length(miss)) stop("not frozen - missing outputs: ", paste(miss, collapse = ", "), " (run 00-07 first)", call. = FALSE)

  dir.create(P$frozen)
  for (f in files) {
    dir.create(dirname(file.path(P$frozen, f)), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(P$hp, f), file.path(P$frozen, f), copy.date = TRUE)
  }
  file.copy(P$provenance, file.path(P$frozen, "provenance.json"))
  file.copy(cfg$config_path, file.path(P$frozen, "config.yaml"))
  all <- sort(list.files(P$frozen, recursive = TRUE))
  sums <- vapply(file.path(P$frozen, all), function(f) digest::digest(file = f, algo = "sha256"), "")
  writeLines(paste0(sums, "  ", all), file.path(P$frozen, FROZEN_MANIFEST))
  commit <- tryCatch(system2("git", c("-C", shQuote(H2M_SCRIPTS), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
                     error = function(e) "unknown")
  dirty <- tryCatch(length(system2("git", c("-C", shQuote(H2M_SCRIPTS), "status", "--porcelain", "."), stdout = TRUE, stderr = FALSE)) > 0,
                    error = function(e) NA)
  writeLines(c(sprintf("frozen %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               sprintf("ensembl_release %s (%s)", cfg$ensembl_release, cfg$ensembl_release_tag),
               sprintf("scripts git commit %s%s", fcoalesce(commit[1], "unknown"), if (isTRUE(dirty)) " (uncommitted changes)" else ""),
               sprintf("files %d", length(all))), file.path(P$frozen, "FROZEN_AT.txt"))
  Sys.chmod(file.path(P$frozen, c(all, FROZEN_MANIFEST, "FROZEN_AT.txt")), mode = "0444")
  log("froze %d files into %s", length(all), P$frozen)
  verify_frozen(P, log)
  write_provenance(P, cfg, "07b_freeze", outputs = file.path(P$frozen, FROZEN_MANIFEST),
                   extra = list(frozen_dir = P$frozen, n_files = length(all), scripts_commit = commit[1]))
}

main()
