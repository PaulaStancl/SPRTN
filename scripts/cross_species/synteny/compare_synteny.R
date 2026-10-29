# ---------------------------------------------------------------------------
# compare_synteny.R - our gene-order synteny blocks (h2m step 03) vs public human-mouse synteny datasets.
# Human GRCh38 vs mouse GRCm39 (the assemblies of Ensembl Compara 110 and UCSC mm39).
#
# Run the whole script:   ./run_synteny.sh            (or: Rscript compare_synteny.R --config config.yaml)
# Run it step by step:     open it in RStudio, set SYNTENY_DIR in section 0, then run section by section.
#                          Every section leaves its result in a named object you can look at.
#
# Datasets (downloaded into <resources_dir>/synteny/ when missing, ~0.3 GB)
#   ours              h2m human_predictions/synteny_blocks.GRCm39.tsv: runs of 1:1 orthologs in conserved order
#   ensembl_compara   Ensembl Compara "Hsap-Mmus synteny" (MySQL dump tables synteny_region, dnafrag_region,
#                     dnafrag, method_link, method_link_species_set, genome_db); release 110: 367 regions
#   ucsc_chains       UCSC liftOver chains hg38ToMm39 (from h2m 00_resources), human span >= min_block_kb
#   ucsc_synnet       UCSC hg38.mm39.syn.net.gz, top-level fills, human span >= min_block_kb (mouse coordinates of
#                     '-' fills are on the + strand; checked against Ensembl, see section 3d)
#   external          files listed in config.yaml `external` (Synteny Portal, Genomicus, halSynteny, MCScanX, ...)
#
# Measures, for every pair of datasets A and B
#   support_frac       fraction of A's human bp covered by B blocks on the same mouse chromosome that also
#                      overlap A's mouse interval
#   orientation_agree  of those bp, the fraction where A and B have the same orientation (+ / -)
#
# Output (config out_dir, default <results_base>/cross_species/synteny/)
#   synteny_datasets.tsv  synteny_dataset_summary.tsv  synteny_pairwise.tsv
#   synteny_compare.html  logs/compare_synteny.log  provenance.json
# ---------------------------------------------------------------------------


# ==== 0. settings ==========================================================================================
# Rscript sets SYNTENY_DIR itself; when running interactively, set it to the folder of this script.
args <- commandArgs(FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
SYNTENY_DIR <- if (length(script_file)) dirname(normalizePath(script_file)) else
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/cross_species/synteny"

cli <- commandArgs(TRUE)
CONFIG <- if ("--config" %in% cli) cli[which(cli == "--config") + 1] else file.path(SYNTENY_DIR, "config.yaml")

source(file.path(SYNTENY_DIR, "..", "h2m", "scripts", "common.R"))   # data.table, yaml, shared helpers
suppressPackageStartupMessages(library(ggplot2))

BUILD <- "GRCm39"
UCSC_SYNNET_URL <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm39/hg38.mm39.syn.net.gz"
COMPARA_URL <- "https://ftp.ensembl.org/pub/release-%s/mysql/ensembl_compara_%s/%s.txt.gz"
COMPARA_TABLES <- c("method_link", "method_link_species_set", "genome_db", "synteny_region", "dnafrag_region", "dnafrag")
MOUSE_CHROMS <- c(paste0("chr", 1:19), "chrX", "chrY")
BLOCK_COLUMNS <- c("h_chrom", "h_start", "h_end", "m_chrom", "m_start", "m_end", "strand")


# ==== 1. helper functions ==================================================================================

# download url to dest unless dest already exists
fetch <- function(url, dest) {
  if (file.exists(dest)) return(dest)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  log("downloading %s", url)
  old <- options(timeout = 36000)
  on.exit(options(old))
  download.file(url, paste0(dest, ".part"), mode = "wb", quiet = TRUE)
  file.rename(paste0(dest, ".part"), dest)
  dest
}

# any table with BLOCK_COLUMNS -> the common block table (1-based coordinates, standard chromosomes only)
as_blocks <- function(d, dataset) {
  d <- d[h_chrom %chin% STD_CHROMS & m_chrom %chin% MOUSE_CHROMS]
  data.table(dataset = dataset,
             block_id = paste0(dataset, "_", seq_len(nrow(d))),
             h_chrom = d$h_chrom, h_start = as.numeric(d$h_start), h_end = as.numeric(d$h_end),
             m_chrom = d$m_chrom, m_start = as.numeric(d$m_start), m_end = as.numeric(d$m_end),
             strand = d$strand,
             n_genes = if ("n_genes" %chin% names(d)) as.integer(d$n_genes) else NA_integer_)
}

# one table of an Ensembl Compara MySQL dump (tab-separated, no header); cols = names of the first columns
read_compara_table <- function(dir, table, cols) {
  f <- file.path(dir, paste0(table, ".txt.gz"))
  fread(cmd = paste("gzip -dc", shQuote(f)), header = FALSE, sep = "\t", quote = "",
        select = seq_along(cols), col.names = cols)
}

# UCSC chain headers -> one row per chain, mouse coordinates converted to the + strand
read_chain_headers <- function(f) {
  ch <- fread(cmd = sprintf("gzip -dc %s | grep '^chain'", shQuote(f)), header = FALSE, sep = " ",
              col.names = c("x", "score", "tName", "tSize", "tStrand", "tStart", "tEnd",
                            "qName", "qSize", "qStrand", "qStart", "qEnd", "id"))
  ch[, .(h_chrom = tName, h_start = tStart + 1, h_end = tEnd,
         m_chrom = qName,
         m_start = fifelse(qStrand == "-", qSize - qEnd + 1, qStart + 1),
         m_end   = fifelse(qStrand == "-", qSize - qStart, qEnd),
         strand  = qStrand)]
}

# top-level fills of a UCSC net file (lines "net <chrom> ..." and " fill ..." with one leading space)
read_net_fills <- function(f) {
  awk <- "awk '/^net /{c=$2} /^ fill /{print c\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}'"
  fread(cmd = paste("gzip -dc", shQuote(f), "|", awk), header = FALSE, sep = "\t",
        col.names = c("h_chrom", "tStart", "tSize", "m_chrom", "strand", "qStart", "qSize"))
}

# bp of A supported by B (same mouse chromosome, overlapping mouse interval) and orientation agreement
support_and_orientation <- function(A, B) {
  b <- B[, .(h_chrom, s = h_start, e = h_end, b_m_chrom = m_chrom, b_m_start = m_start, b_m_end = m_end, b_strand = strand)]
  setkey(b, h_chrom, s, e)
  a <- A[, .(block_id, h_chrom, s = h_start, e = h_end, m_chrom, m_start, m_end, strand)]
  hits <- foverlaps(a, b, by.x = c("h_chrom", "s", "e"), nomatch = NULL)        # human overlaps
  hits <- hits[b_m_chrom == m_chrom & b_m_start <= m_end & b_m_end >= m_start]  # ... that agree in mouse
  hits[, bp := pmin(e, i.e) - pmax(s, i.s) + 1]
  per_block <- hits[, .(supported_bp = min(sum(bp), i.e[1] - i.s[1] + 1),       # capped at A's block length
                        same_strand_bp = sum(bp[b_strand == strand]),
                        all_bp = sum(bp)), by = block_id]
  list(support_frac = sum(per_block$supported_bp) / sum(A$h_end - A$h_start + 1),
       orientation_agree = if (sum(per_block$all_bp) > 0) sum(per_block$same_strand_bp) / sum(per_block$all_bp) else NA_real_)
}

compare_pair <- function(A, B) {
  so <- support_and_orientation(A, B)
  data.table(support_frac = round(so$support_frac, 3), orientation_agree = round(so$orientation_agree, 3))
}

html_table <- function(d) {
  head <- paste0("<tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>")
  rows <- apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table>", head, paste(rows, collapse = ""), "</table>")
}

png_tag <- function(p, width, height) {
  f <- tempfile(fileext = ".png")
  ggsave(f, p, width = width, height = height, dpi = 110, bg = "white")
  sprintf("<img style='max-width:100%%' src='%s'/>", base64enc::dataURI(file = f, mime = "image/png"))
}


# ==== 2. configuration and paths ===========================================================================
scfg <- read_yaml(CONFIG)
h2m_config <- if (startsWith(scfg$h2m_config, "/")) scfg$h2m_config else file.path(dirname(normalizePath(CONFIG)), scfg$h2m_config)
cfg <- load_config(normalizePath(h2m_config))                  # shared paths: resources, orthologs, chains
P <- paths(cfg)
M <- P$mouse[[BUILD]]
blocks_file <- file.path(P$hp, sprintf("synteny_blocks.%s.tsv", BUILD))   # written by h2m step 03

scalars <- lapply(Filter(function(v) length(v) == 1 && (is.character(v) || is.numeric(v)), cfg), as.character)
OUT <- subst(scfg$out_dir, scalars)
if (grepl(CLOUD_RE, normalizePath(OUT, mustWork = FALSE), ignore.case = TRUE)) stop("refusing a cloud-synced out_dir: ", OUT)
dir.create(file.path(OUT, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "tmp"), recursive = TRUE, showWarnings = FALSE)
P_out <- list(logs = file.path(OUT, "logs"), provenance = file.path(OUT, "provenance.json"))
log <- make_logger("compare_synteny", P_out)

MIN_BP <- as.numeric(scfg$min_block_kb) * 1e3
SYNTENY_RES <- file.path(P$resources, "synteny")
log("out_dir %s", OUT)
log("h2m config %s; minimum UCSC block %g kb", h2m_config, scfg$min_block_kb)


# ==== 3. load the datasets =================================================================================
datasets <- list()

# 3a. ours - gene-order blocks from h2m step 03
if (!file.exists(blocks_file)) stop("missing ", blocks_file, " - run h2m step 03 first")
ours <- fread(blocks_file)
setnames(ours, "orientation", "strand")
datasets$ours <- as_blocks(ours, "ours")

# 3b. Ensembl Compara synteny regions (human and mouse side of each region in dnafrag_region)
compara_dir <- file.path(SYNTENY_RES, paste0("ensembl_compara_", cfg$ensembl_release))
for (t in COMPARA_TABLES) fetch(sprintf(COMPARA_URL, cfg$ensembl_release, cfg$ensembl_release, t), file.path(compara_dir, paste0(t, ".txt.gz")))

method_link <- read_compara_table(compara_dir, "method_link", c("method_link_id", "type"))
mlss        <- read_compara_table(compara_dir, "method_link_species_set", c("mlss_id", "method_link_id", "species_set_id", "name"))
genome_db   <- read_compara_table(compara_dir, "genome_db", c("genome_db_id", "taxon_id", "name", "assembly"))
synteny_id  <- mlss[method_link_id %in% method_link[type == "SYNTENY", method_link_id] & grepl("Hsap", name) & grepl("Mmus", name)]
stopifnot(nrow(synteny_id) == 1)
genomes     <- genome_db[name %chin% c("homo_sapiens", "mus_musculus")]
log("Ensembl Compara %s: %s (mlss %s); %s", cfg$ensembl_release, synteny_id$name, synteny_id$mlss_id,
    paste(genomes$name, genomes$assembly, collapse = ", "))

regions <- read_compara_table(compara_dir, "synteny_region", c("synteny_region_id", "mlss_id"))[mlss_id == synteny_id$mlss_id]
pieces  <- read_compara_table(compara_dir, "dnafrag_region", c("synteny_region_id", "dnafrag_id", "start", "end", "strand"))
pieces  <- pieces[synteny_region_id %in% regions$synteny_region_id]
dnafrag <- fread(cmd = sprintf("gzip -dc %s | awk -F'\\t' '$4==%s || $4==%s'", shQuote(file.path(compara_dir, "dnafrag.txt.gz")),
                               genomes$genome_db_id[1], genomes$genome_db_id[2]),
                 header = FALSE, sep = "\t", quote = "", select = 1:4,
                 col.names = c("dnafrag_id", "length", "name", "genome_db_id"), colClasses = list(character = 3))
pieces[dnafrag, on = "dnafrag_id", `:=`(chrom = norm_chrom(i.name), genome_db_id = i.genome_db_id)]
pieces[genomes, on = "genome_db_id", species := i.name]

human_side <- pieces[species == "homo_sapiens", .(synteny_region_id, h_chrom = chrom, h_start = start, h_end = end, h_strand = strand)]
mouse_side <- pieces[species == "mus_musculus", .(synteny_region_id, m_chrom = chrom, m_start = start, m_end = end, m_strand = strand)]
ensembl <- human_side[mouse_side, on = "synteny_region_id", nomatch = NULL]
ensembl[, strand := fifelse(h_strand * m_strand < 0, "-", "+")]
datasets$ensembl_compara <- as_blocks(ensembl, "ensembl_compara")

# 3c. UCSC liftOver chains
chains <- read_chain_headers(M$chain)
datasets$ucsc_chains <- as_blocks(chains[h_end - h_start + 1 >= MIN_BP], "ucsc_chains")

# 3d. UCSC synteny-filtered net: top-level fills. The mouse start of a '-' fill is on the + strand: checked
#     against Ensembl Compara 110 (2026-10-29), 96.8% of '-' fills supported when read that way, 15.6% when
#     converted from the - strand.
synnet_file <- fetch(UCSC_SYNNET_URL, file.path(SYNTENY_RES, basename(UCSC_SYNNET_URL)))
fills <- read_net_fills(synnet_file)
fills <- fills[tSize >= MIN_BP]
fills[, `:=`(h_start = tStart + 1, h_end = tStart + tSize, m_start = qStart + 1, m_end = qStart + qSize)]
datasets$ucsc_synnet <- as_blocks(fills, "ucsc_synnet")

# 3e. external datasets from config.yaml
for (x in scfg$external) {
  d <- fread(x$file)
  need_cols(d, BLOCK_COLUMNS, paste("external", x$name))
  d[, `:=`(h_chrom = norm_chrom(h_chrom), m_chrom = norm_chrom(m_chrom))]
  datasets[[x$name]] <- as_blocks(d, x$name)
}


# ==== 4. dataset summary ===================================================================================
all_blocks <- rbindlist(datasets)
summary_tab <- all_blocks[, .(n_blocks = .N,
                              human_Mb = round(sum(h_end - h_start + 1) / 1e6, 1),
                              median_block_Mb = round(median(h_end - h_start + 1) / 1e6, 2),
                              human_chroms = uniqueN(h_chrom)), by = dataset]
summary_tab
write_tsv(all_blocks, file.path(OUT, "synteny_datasets.tsv"), log)
write_tsv(summary_tab, file.path(OUT, "synteny_dataset_summary.tsv"), log)


# ==== 5. pairwise comparison ===============================================================================
pairs <- CJ(A = names(datasets), B = names(datasets))[A != B]
pairwise <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  cbind(pairs[i], compare_pair(datasets[[pairs$A[i]]], datasets[[pairs$B[i]]]))
}))
pairwise[A == "ours"]
write_tsv(pairwise, file.path(OUT, "synteny_pairwise.tsv"), log)
for (i in which(pairwise$A == "ours"))
  log("ours vs %-16s support %.2f, orientation %.2f",
      pairwise$B[i], pairwise$support_frac[i], pairwise$orientation_agree[i])


# ==== 6. plots =============================================================================================
heatmap <- ggplot(pairwise, aes(B, A, fill = support_frac)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", support_frac)), size = 3.5) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2471a3", limits = c(0, 1)) +
  labs(x = "reference B", y = "dataset A", fill = "A bp\nsupported\nby B",
       title = "Fraction of A's human span covered by consistent B blocks") +
  theme_minimal(base_size = 11)
heatmap

colours <- c(ours = "black", ensembl_compara = "#c0392b", ucsc_chains = "#2471a3", ucsc_synnet = "#27ae60")
extra <- setdiff(names(datasets), names(colours))
colours <- c(colours, setNames(rep("#8e44ad", length(extra)), extra))

dot_plot <- function(human_chrom) {
  mouse_chroms <- all_blocks[h_chrom == human_chrom & dataset == "ours", unique(m_chrom)]
  d <- all_blocks[h_chrom == human_chrom & m_chrom %chin% mouse_chroms]
  d[, `:=`(y = fifelse(strand == "-", m_end, m_start) / 1e6,
           yend = fifelse(strand == "-", m_start, m_end) / 1e6,
           width = fifelse(dataset == "ours", 0.5, 1.3))]
  ggplot(d, aes(x = h_start / 1e6, xend = h_end / 1e6, y = y, yend = yend, colour = dataset)) +
    geom_segment(aes(linewidth = width), alpha = 0.7) +
    scale_linewidth_identity() +
    scale_colour_manual(values = colours) +
    facet_wrap(~ m_chrom, scales = "free_y") +
    labs(x = paste(human_chrom, "(human, Mb)"), y = "mouse position (Mb)",
         title = paste(human_chrom, "- all datasets, one panel per mouse chromosome")) +
    theme_minimal(base_size = 9) + theme(legend.position = "top")
}
dot_plots <- lapply(setNames(nm = intersect(unlist(scfg$plot_chroms), all_blocks$h_chrom)), dot_plot)
dot_plots[[1]]


# ==== 7. HTML report and provenance =======================================================================
html_file <- file.path(OUT, "synteny_compare.html")
writeLines(c(
  "<!doctype html><html><head><meta charset='utf-8'><title>Synteny comparison</title>",
  "<style>body{font:14px system-ui,sans-serif;max-width:1150px;margin:auto;padding:16px}",
  "td,th{border:1px solid #ccc;padding:3px 8px}table{border-collapse:collapse;margin:8px 0}</style></head><body>",
  "<h1>Human-mouse synteny: our gene-order blocks vs public datasets (GRCh38 / GRCm39)</h1>",
  sprintf("<p>UCSC chains and synNet fills shorter than %g kb on the human side are left out. ", scfg$min_block_kb),
  "Support = fraction of A's human span covered by B blocks on the same mouse chromosome that overlap A's mouse interval.</p>",
  "<h2>Datasets</h2>", html_table(summary_tab),
  "<h2>Pairwise</h2>", png_tag(heatmap, 7, 5), html_table(pairwise),
  "<h2>Dot plots</h2>", vapply(dot_plots, png_tag, "", width = 10, height = 6.5),
  "</body></html>"), html_file)
log("wrote %s", html_file)

cfg$config_path <- normalizePath(CONFIG)
write_provenance(P_out, cfg, "compare_synteny",
                 resources = list(chain = M$chain, synnet = synnet_file),
                 outputs = list.files(OUT, pattern = "^synteny_", full.names = TRUE),
                 extra = list(min_block_kb = scfg$min_block_kb, datasets = names(datasets)))
log("done")
