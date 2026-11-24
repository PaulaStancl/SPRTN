#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24
#
# 00_input_inspection.R - Task 0: what is in the mouse WGS callsets, before any analysis.
#
# Read-only. The input BCF/VCFs are never changed and NOT re-normalized: the small-variant
# files were already split (bcftools norm -m -any, no -f, so indels were not left-aligned)
# and, for setting A, merged across callers with INFO/CALLER. There is no mm10 FASTA here.
#
# Settings (folders under data/mouse/wgs):
#   A  tumor_vs_normal/snv_indel/   tumour vs matched non-adjacent KI normal; Strelka2 + Mutect2 + FreeBayes merged
#   B  tumor_vs_joint_wt/           tumour vs joint WT normals, Mutect2 (exploratory)
#   C  tumor_only_pon/              tumour-only + WT panel of normals, Mutect2
#   SV tumor_vs_normal/sv/calls/ (Delly, SvABA, Manta) and sv/merged/ (OctopuSV)
#
# What it does, block by block:
#   0  EDIT THIS BLOCK: paths, samples
#   1  setup and the bcftools helper
#   2  provenance files (processing notes, md5 lists) and md5 check
#   3  file inventory: setting, tumour, caller per file
#   4  headers: samples, contigs, command lines, annotation, INFO/FORMAT/FILTER fields
#   5  which normals each setting actually used (B first)
#   6  small-variant records: FILTER, CALLER, variant type, multiallelics, duplicate IDs
#   7  SV records: SVTYPE, FILTER, caller support (SOURCES)
#   8  manifest (one row per sample x file)
#   9  issues table
#   10 figure, outputs, summary.md, RUN_LOG.md
#
# Run in RStudio line by line, or: Rscript scripts/mouse_wgs/00_input_inspection.R
# Env on ssi: r-mutation (needs data.table, ggplot2 and bcftools on the PATH).
# qcVCF is not used here (headers and record counts only); it starts in Task 1.
# Output: results/mouse_wgs/00_input_inspection/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
ROOT     <- file.path(REPO, "data/mouse/wgs")
OUT      <- file.path(REPO, "results/mouse_wgs/00_input_inspection")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/00_input_inspection.R")
# bcftools: RStudio often does not see the micromamba env's PATH, so Sys.which() can be "".
# Then the env's own bin/bcftools is used. Check with: system2(BCFTOOLS, "--version")
BCFTOOLS <- Sys.which("bcftools")
if (!nzchar(BCFTOOLS)) BCFTOOLS <- "/common/WORK/pstancl/envs/r-mutation/bin/bcftools"

# Tumours and their matched KI normals (same mouse); names as in the file names
TUMOURS <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1",
             HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
             HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")
KI_NORMALS <- c(HCC_1 = "sfb_dikic_2025_05_05_KI1",
                HCC_2 = "sfb_dikic_2025_05_06_KI2",
                HCC_3 = "sfb_dikic_2025_05_07_KI3")

CHECK_MD5  <- TRUE     # recompute md5 of every file (a few minutes for the BCFs)

# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# Console as wide as the terminal (falls back to 80 when there is no terminal, e.g. batch Rscript)
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()

# bcftools helper: returns stdout lines, stops on a non-zero exit.
# Only needed later, to read the BCF files (headers, records); not required to start.
bcf <- function(args) {
  if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS, ": install it (micromamba install -c bioconda bcftools) or set BCFTOOLS in block 0")
  err <- tempfile()
  out <- tryCatch(suppressWarnings(system2(BCFTOOLS, args, stdout = TRUE, stderr = err)),
                  error = function(e) structure(character(), status = -1L, msg = conditionMessage(e)))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0)
    stop("bcftools failed: bcftools ", paste(args, collapse = " "), "\n",
         paste(c(attr(out, "msg"), readLines(err, warn = FALSE)), collapse = "\n"))
  out
}
bcftools_version <- tryCatch(bcf("--version")[1:2], error = function(e) c("bcftools not available", ""))

# ---- 2. Provenance -----------------------------------------------------------
# 2a. Processing notes from the data provider: read them before anything else
processed_txt <- readLines(file.path(ROOT, "patrik_mouse_wgs_processed.txt"), warn = FALSE)
cat(processed_txt, sep = "\n")

# 2b. md5 lists ("<md5>  <path>"), one row per listed file
md5_files <- list.files(ROOT, pattern = "^md5sums.*\\.txt$", full.names = TRUE)
md5_listed <- rbindlist(lapply(md5_files, function(f) {
  l <- trimws(readLines(f, warn = FALSE))
  l <- l[nzchar(l)]
  data.table(md5_list = basename(f),
             md5_expected = sub("\\s.*$", "", l),
             listed_path = sub("^\\S+\\s+\\*?", "", l))
}))
md5_listed[, file := basename(listed_path)]
md5_listed

# ---- 3. File inventory -------------------------------------------------------
files <- data.table(path = list.files(ROOT, recursive = TRUE, full.names = TRUE))
files[, rel := substring(path, nchar(ROOT) + 2)]
files[, file := basename(path)]
files[, size_mb := round(file.size(path) / 1e6, 2)]
files[, setting := fcase(
  startsWith(rel, "tumor_vs_normal/snv_indel/"), "A",
  startsWith(rel, "tumor_vs_joint_wt/"),         "B",
  startsWith(rel, "tumor_only_pon/"),            "C",
  startsWith(rel, "tumor_vs_normal/sv/calls/"),  "SV_raw",
  startsWith(rel, "tumor_vs_normal/sv/merged/"), "SV_merged",
  default = "other")]
files[, is_index := grepl("\\.(csi|tbi)$", file)]
files[, is_vcf := !is_index & grepl("\\.(bcf|vcf|vcf\\.gz|svcf|svcf\\.gz)$", file)]

# Tumour named in the file name (HCC_1/2/3)
files[, tumour := NA_character_]
for (k in names(TUMOURS)) files[grepl(TUMOURS[[k]], file, fixed = TRUE), tumour := k]

# Caller named in the file name (SV files; setting A is a merge, its callers come from INFO/CALLER)
files[, caller_in_name := fcase(
  grepl("delly", file, ignore.case = TRUE),    "delly",
  grepl("svaba", file, ignore.case = TRUE),    "svaba",
  grepl("manta", file, ignore.case = TRUE),    "manta",
  grepl("octopus", file, ignore.case = TRUE),  "octopusv",
  grepl("mutect2", file, ignore.case = TRUE),  "mutect2",
  grepl("strelka", file, ignore.case = TRUE),  "strelka",
  grepl("freebayes", file, ignore.case = TRUE), "freebayes",
  default = NA_character_)]

# Is there an index next to each variant file? (needed only for region queries; we use -t)
files[, has_index := is_vcf & (file.exists(paste0(path, ".csi")) | file.exists(paste0(path, ".tbi")))]

# 3b. md5 check against the lists
if (CHECK_MD5) files[, md5 := unname(tools::md5sum(path))] else files[, md5 := NA_character_]
files <- md5_listed[, .(file, md5_expected, md5_list)][files, on = "file"]
files[, md5_status := fcase(
  is.na(md5_expected), "not listed",
  is.na(md5),          "not checked",
  md5 == md5_expected, "ok",
  default = "MISMATCH")]
md5_missing <- md5_listed[!file %chin% files$file]     # listed, but not on disk

files[, .N, by = .(setting, is_vcf)]
files[is_vcf == TRUE, .(setting, tumour, caller_in_name, file, size_mb, has_index, md5_status)]
files[, .N, by = md5_status]
md5_missing

# ---- 4. Headers --------------------------------------------------------------
vcfs <- files[is_vcf == TRUE]

# Header lines. bcftools for BCF/VCF; plain text for anything bcftools cannot open (e.g. SVCF)
read_header <- function(f) {
  h <- tryCatch(bcf(c("view", "-h", shQuote(f))), error = function(e) NULL)
  if (is.null(h))
    h <- system(paste("zcat -f", shQuote(f), "| awk '/^#/{print; next} {exit}'"), intern = TRUE)
  h
}
headers <- setNames(lapply(vcfs$path, read_header), vcfs$rel)
lengths(headers)

# 4a. Sample columns (from the #CHROM line: columns 10 onwards)
sample_cols <- rbindlist(lapply(names(headers), function(r) {
  chrom_line <- tail(grep("^#CHROM", headers[[r]], value = TRUE), 1)
  s <- if (length(chrom_line)) strsplit(chrom_line, "\t", fixed = TRUE)[[1]] else character(0)
  if (length(s) < 10) return(data.table(rel = r, col = NA_integer_, sample = NA_character_,
                                        has_chrom_line = length(chrom_line) > 0))
  data.table(rel = r, col = seq_along(s[-(1:9)]), sample = s[-(1:9)], has_chrom_line = TRUE)
}))
# Files with no #CHROM line: header could not be read (e.g. a BCF without a working bcftools)
sample_cols[has_chrom_line == FALSE, rel]
sample_cols[, .(samples = paste(sample, collapse = ", ")), by = rel]

# 4b. Attribute of a header line, e.g. ID / Number / Type / Description / Version
hdr_attr <- function(x, key) {
  pat <- paste0('.*[<,]', key, '=("([^"]*)"|([^,>]*)).*')
  fifelse(grepl(paste0('[<,]', key, '='), x), sub(pat, "\\2\\3", x), NA_character_)
}

# INFO / FORMAT / FILTER definitions, one row per field per file
fields <- rbindlist(lapply(names(headers), function(r) {
  rbindlist(lapply(c("INFO", "FORMAT", "FILTER"), function(kind) {
    x <- grep(paste0("^##", kind, "=<"), headers[[r]], value = TRUE)
    if (!length(x)) return(NULL)
    data.table(rel = r, kind = kind, ID = hdr_attr(x, "ID"), Number = hdr_attr(x, "Number"),
               Type = hdr_attr(x, "Type"), Description = hdr_attr(x, "Description"))
  }))
}))
dcast(fields, rel ~ kind, value.var = "ID", fun.aggregate = length)

# Fields this project depends on: present in which file?
key_fields <- c(FORMAT_AD = "AD", FORMAT_AF = "AF", FORMAT_DP = "DP",
                INFO_CALLER = "CALLER", INFO_ANN = "ANN", INFO_CSQ = "CSQ",
                INFO_SVTYPE = "SVTYPE", INFO_SOURCES = "SOURCES")
has_field <- vcfs[, .(rel, setting)]
for (k in names(key_fields)) {
  knd <- sub("_.*", "", k)                     # INFO / FORMAT
  hit <- fields[kind == knd & ID == key_fields[[k]], unique(rel)]
  has_field[, (k) := rel %chin% hit]
}
has_field

# 4c. Contigs: naming style, count, lengths identical across files?
contigs <- rbindlist(lapply(names(headers), function(r) {
  x <- grep("^##contig=<", headers[[r]], value = TRUE)
  if (!length(x)) return(data.table(rel = r, contig = NA_character_, length = NA_real_))
  data.table(rel = r, contig = hdr_attr(x, "ID"), length = as.numeric(hdr_attr(x, "length")))
}))
contig_summary <- contigs[, .(n_contigs = sum(!is.na(contig)),
                              chr_prefix = any(startsWith(contig, "chr"), na.rm = TRUE),
                              contig_set = paste(contig, length, sep = ":", collapse = ",")),
                          by = rel]
contig_summary[, contig_set_id := as.integer(factor(contig_set))]   # same number = identical contig list
contig_summary[, contig_set := NULL]
contig_summary

# 4d. Other meta lines (##key=value): reference, source, annotation, command lines
meta <- rbindlist(lapply(names(headers), function(r) {
  x <- grep("^##", headers[[r]], value = TRUE)
  x <- x[!grepl("^##(INFO|FORMAT|FILTER|contig|ALT)=<", x)]
  if (!length(x)) return(NULL)
  data.table(rel = r, key = sub("^##([^=]+)=.*", "\\1", x), value = sub("^##[^=]+=", "", x))
}))
meta[, .N, by = key][order(-N)]
meta[key %chin% c("reference", "source", "fileformat")]

# GATK command lines (Mutect2, FilterMutectCalls, ...): tool, version, command
gatk_cmd <- meta[key == "GATKCommandLine",
                 .(rel, tool = hdr_attr(value, "ID"), version = hdr_attr(value, "Version"),
                   date = hdr_attr(value, "Date"), command = hdr_attr(value, "CommandLine"))]
gatk_cmd[, .(rel, tool, version)]

# Values of one or more options in a command line, e.g. --normal-sample / -normal
opt_values <- function(cmd, opts) {
  pat <- paste0("(^|\\s)(", paste(opts, collapse = "|"), ")[ =]+(\\S+)")
  m <- regmatches(cmd, gregexpr(pat, cmd))
  vapply(m, function(v) if (length(v)) paste(unique(sub(pat, "\\3", v)), collapse = ", ") else NA_character_, "")
}
gatk_cmd[, `:=`(
  normals_cmd   = opt_values(command, c("--normal-sample", "-normal")),
  tumour_cmd    = opt_values(command, c("--tumor-sample", "-tumor")),
  inputs_cmd    = opt_values(command, c("--input", "-I")),
  pon_cmd       = opt_values(command, c("--panel-of-normals", "-pon")),
  germline_cmd  = opt_values(command, c("--germline-resource")),
  ob_priors_cmd = opt_values(command, c("--orientation-bias-artifact-priors", "-ob-priors")),
  intervals_cmd = opt_values(command, c("--intervals", "-L")))]
gatk_cmd[tool == "Mutect2", .(rel, version, normals_cmd, tumour_cmd, pon_cmd, germline_cmd)]
gatk_cmd[tool == "FilterMutectCalls", .(rel, version, ob_priors_cmd)]

# Other tools' command lines and versions (bcftools, SnpEff, VEP, Strelka, Delly, Manta, OctopuSV ...)
tool_lines <- meta[grepl("Command|Version|source|VEP|SnpEff|octopus|cmdline", key, ignore.case = TRUE) &
                     key != "GATKCommandLine"]
tool_lines[, .N, by = key]

# 4e. Annotation type and version per file
annotation <- has_field[, .(rel, annotation = fcase(INFO_ANN, "SnpEff ANN", INFO_CSQ, "VEP CSQ", default = "none"))]
annotation[rel %chin% fields[kind == "INFO" & grepl("refGene|ANNOVAR", ID), rel], annotation := "ANNOVAR"]
ann_version <- meta[key %chin% c("SnpEffVersion", "VEP", "ANNOVAR"), .(ann_version = paste(value, collapse = " | ")), by = rel]
annotation <- ann_version[annotation, on = "rel"]
annotation[, .N, by = .(annotation, ann_version)]

# ---- 5. Which normals did each setting actually use? -------------------------
# From three places: the sample columns, ##normal_sample / ##tumor_sample lines (Mutect2),
# and the Mutect2 command line. Setting B's file names mention the KI normal (e.g. _vs_..._KI1):
# check whether the KI normal really was a normal there, or only the WT livers.
all_tumour_names <- unname(TUMOURS)
normals_used <- vcfs[setting %chin% c("A", "B", "C"), .(rel, setting, tumour)]
normals_used[, ki_normal_expected := KI_NORMALS[tumour]]
normals_used[, ki_normal_in_name := mapply(function(f, n) grepl(n, f, fixed = TRUE), basename(rel), ki_normal_expected)]
normals_used[sample_cols[, .(cols = paste(sample, collapse = ", "),
                             non_tumour_cols = paste(setdiff(sample, all_tumour_names), collapse = ", ")), by = rel],
             on = "rel", `:=`(sample_columns = i.cols, non_tumour_columns = i.non_tumour_cols)]
normals_used[meta[key == "normal_sample", .(v = paste(value, collapse = ", ")), by = rel],
             on = "rel", normal_sample_lines := i.v]
normals_used[meta[key == "tumor_sample", .(v = paste(value, collapse = ", ")), by = rel],
             on = "rel", tumor_sample_lines := i.v]
normals_used[gatk_cmd[tool == "Mutect2", .(n = paste(na.omit(normals_cmd), collapse = ", "),
                                           p = paste(na.omit(pon_cmd), collapse = ", ")), by = rel],
             on = "rel", `:=`(normals_mutect2_cmd = i.n, pon_mutect2_cmd = i.p)]
# Is the matched KI normal among the normals actually used?
normals_used[, ki_normal_used := mapply(function(n, a, b, c) any(grepl(n, c(a, b, c), fixed = TRUE)),
                                        ki_normal_expected, non_tumour_columns,
                                        fcoalesce(normal_sample_lines, ""), fcoalesce(normals_mutect2_cmd, ""))]
setorder(normals_used, setting, tumour)
normals_used[setting == "B"]     # report this first
normals_used

# WT normal names: every normal that is neither a tumour nor a matched KI normal
wt_names <- setdiff(unique(c(unlist(strsplit(na.omit(normals_used$non_tumour_columns), ", ")),
                             unlist(strsplit(na.omit(normals_used$normal_sample_lines), ", ")),
                             unlist(strsplit(na.omit(normals_used$normals_mutect2_cmd), ", ")))),
                    c(unname(KI_NORMALS), all_tumour_names, ""))
wt_names

# ---- 6. Small-variant records ------------------------------------------------
# One row per record: CHROM POS REF ALT FILTER (+ CALLER when the file has it).
read_small <- function(path, r) {
  has_caller <- has_field[rel == r, INFO_CALLER]
  fmt <- paste0("%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER", if (has_caller) "\\t%INFO/CALLER", "\\n")
  cmd <- paste(shQuote(BCFTOOLS), "query",
               "-f", shQuote(fmt), shQuote(path))
  cn <- c("CHROM", "POS", "REF", "ALT", "FILTER", if (has_caller) "CALLER")
  x <- tryCatch(fread(cmd = cmd, sep = "\t", header = FALSE, col.names = cn,
                      colClasses = list(character = cn[cn != "POS"], integer = "POS")),
                error = function(e) NULL)
  if (is.null(x) || !nrow(x)) x <- setnames(data.table(matrix(character(0), 0, length(cn))), cn)
  if (!has_caller) x[, CALLER := "mutect2"]          # B and C are Mutect2 only
  x[, POS := as.integer(POS)]
  x
}
small_files <- vcfs[setting %chin% c("A", "B", "C")]
small <- setNames(Map(read_small, small_files$path, small_files$rel), small_files$rel)
sapply(small, nrow)          # row counts per file

small_dt <- rbindlist(small, idcol = "rel")
small_dt <- small_files[, .(rel, setting, tumour)][small_dt, on = "rel"]
small_dt[, ID := paste(CHROM, POS, REF, ALT, sep = ":")]
small_dt[, TYPE := fcase(
  grepl(",", ALT, fixed = TRUE),                        "multiallelic",
  ALT %chin% c("*", "."),                               "star_or_missing",
  nchar(REF) == 1 & nchar(ALT) == 1,                    "SNV",
  nchar(REF) == nchar(ALT),                             "MNV",
  nchar(REF) < nchar(ALT) & startsWith(ALT, substr(REF, 1, 1)) & nchar(REF) == 1, "INS",
  nchar(REF) > nchar(ALT) & startsWith(REF, substr(ALT, 1, 1)) & nchar(ALT) == 1, "DEL",
  default = "complex")]
small_dt[, PASS := FILTER == "PASS"]

# 6a. Records, PASS and type per file
records <- small_dt[, .(records = .N, PASS = sum(PASS), nonPASS = sum(!PASS),
                        SNV = sum(TYPE == "SNV"), MNV = sum(TYPE == "MNV"), INS = sum(TYPE == "INS"),
                        DEL = sum(TYPE == "DEL"), complex = sum(TYPE == "complex"),
                        multiallelic = sum(TYPE == "multiallelic"), star = sum(TYPE == "star_or_missing")),
                    by = .(setting, tumour, rel)]
setorder(records, setting, tumour)
records

# 6b. FILTER reasons (a record can have several, separated by ";")
filter_reasons <- small_dt[, .(reason = unlist(strsplit(FILTER, ";", fixed = TRUE))), by = .(setting, tumour, CALLER)][
  , .N, by = .(setting, tumour, CALLER, reason)]
setorder(filter_reasons, setting, tumour, CALLER, -N)
filter_reasons

# Filtered? FilterMutectCalls in the header, and non-PASS records kept in the file
filter_state <- small_files[, .(rel, setting, tumour)]
filter_state[, filtermutectcalls := rel %chin% gatk_cmd[tool == "FilterMutectCalls", rel]]
filter_state[records, on = "rel", `:=`(nonPASS_kept = i.nonPASS > 0, records = i.records)]
filter_state

# 6c. Setting A: which callers per record (INFO/CALLER, e.g. "mutect2,strelka")
caller_combos <- small_dt[setting == "A", .N, by = .(tumour, CALLER, PASS)]
setorder(caller_combos, tumour, -N)
caller_combos

# 6d. Duplicate chrom:pos:ref:alt IDs within a file (should be 0 after the merge)
dups <- small_dt[, .(n = .N, callers = paste(sort(unique(CALLER)), collapse = "|")), by = .(rel, ID)][n > 1]
dups[, .N, by = rel]
head(dups)

# ---- 7. SV records -----------------------------------------------------------
# Raw caller files and OctopuSV merged files. Read as text (SVCF may not open in bcftools);
# only the 8 fixed columns are used. No normalization for SVs.
read_sv <- function(path) {
  src <- if (grepl("\\.bcf$", path)) paste(shQuote(BCFTOOLS), "view", shQuote(path)) else paste("zcat -f", shQuote(path))
  cmd <- paste(src, "| grep -v '^##'")
  x <- tryCatch(fread(cmd = cmd, sep = "\t", header = TRUE, colClasses = "character",
                      select = 1:8), error = function(e) NULL)
  if (is.null(x)) return(data.table())
  setnames(x, 1:8, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
  info_val <- function(key) fifelse(grepl(paste0("(^|;)", key, "="), x$INFO),
                                    sub(paste0("^(.*;)?", key, "=([^;]*).*$"), "\\2", x$INFO), NA_character_)
  x[, `:=`(SVTYPE = info_val("SVTYPE"), SVLEN = info_val("SVLEN"), END = info_val("END"),
           SOURCES = info_val("SOURCES"), SUPPORT = info_val("SUPPORT"),
           IMPRECISE = grepl("(^|;)IMPRECISE(;|$)", INFO))]
  x[, INFO := NULL]
  x
}
sv_files <- vcfs[setting %chin% c("SV_raw", "SV_merged")]
sv <- setNames(lapply(sv_files$path, read_sv), sv_files$rel)
sapply(sv, nrow)

sv_dt <- rbindlist(sv, idcol = "rel", fill = TRUE)
sv_dt <- sv_files[, .(rel, setting, tumour, caller_in_name)][sv_dt, on = "rel"]

# 7a. Records per file, SVTYPE and FILTER
sv_counts <- sv_dt[, .N, by = .(setting, tumour, caller_in_name, SVTYPE, PASS = FILTER == "PASS")]
setorder(sv_counts, setting, tumour, caller_in_name, SVTYPE)
dcast(sv_counts, setting + tumour + caller_in_name + PASS ~ SVTYPE, value.var = "N", fill = 0)

# 7b. OctopuSV caller support: SOURCES values and number of callers per event
if ("SOURCES" %in% names(sv_dt)) {
  sv_dt[!is.na(SOURCES), n_callers := lengths(strsplit(SOURCES, ",", fixed = TRUE))]
  print(sv_dt[setting == "SV_merged", .N, by = .(tumour, SOURCES, n_callers)][order(tumour, -N)])
}

# 7c. SV inventory: one row per SV file, with format, samples, versions
sv_inventory <- sv_files[, .(rel, setting, tumour, caller_in_name, size_mb, md5_status)]
sv_inventory[, fileformat := sapply(rel, function(r) sub("^##fileformat=", "", grep("^##fileformat=", headers[[r]], value = TRUE)[1]))]
sv_inventory[sample_cols[, .(s = paste(sample, collapse = ", ")), by = rel], on = "rel", samples := i.s]
sv_inventory[tool_lines[, .(t = paste(paste0(key, "=", substr(value, 1, 200)), collapse = " | ")), by = rel],
             on = "rel", tool_header_lines := i.t]
sv_inventory[sv_dt[, .(records = .N, PASS = sum(FILTER == "PASS")), by = rel], on = "rel",
             `:=`(records = i.records, PASS = i.PASS)]
sv_inventory[, has_SOURCES := rel %chin% has_field[INFO_SOURCES == TRUE, rel]]
sv_inventory

# ---- 8. Manifest -------------------------------------------------------------
# One row per sample column per file: who the sample is and what role it has in that file
manifest <- vcfs[, .(rel, setting, tumour, file, md5, md5_status)][sample_cols, on = "rel"]
manifest[, sample_class := fcase(
  sample %chin% unname(TUMOURS),    "KI tumour",
  sample %chin% unname(KI_NORMALS), "KI normal (non-adjacent liver)",
  sample %chin% wt_names,           "WT normal liver",
  sample %chin% c("TUMOR", "TUMOUR"), "tumour (caller placeholder name)",
  sample %chin% c("NORMAL"),          "normal (caller placeholder name)",
  default = "unknown - check")]
manifest[, animal := fcase(
  sample %chin% unname(TUMOURS),    names(TUMOURS)[match(sample, TUMOURS)],
  sample %chin% unname(KI_NORMALS), names(KI_NORMALS)[match(sample, KI_NORMALS)],
  default = sample)]
manifest[, animal := sub("^HCC_", "KI mouse ", animal)]
manifest[, tissue := fifelse(grepl("tumour|TUMOR", sample_class), "liver tumour", "liver")]
manifest[, genotype := fifelse(sample_class == "WT normal liver", "WT", fifelse(grepl("^KI", sample_class), "KI", NA_character_))]
manifest[, role_in_file := fifelse(grepl("tumour", sample_class), "tumour", "normal")]
manifest[, pairing := fifelse(is.na(tumour), NA_character_, paste0(TUMOURS[tumour], " | ", KI_NORMALS[tumour]))]
manifest[gatk_cmd[tool == "Mutect2", .(v = paste(unique(version), collapse = ",")), by = rel],
         on = "rel", mutect2_version := i.v]
manifest[small_dt[, .(c = paste(sort(unique(unlist(strsplit(CALLER, ",", fixed = TRUE)))), collapse = ",")), by = rel],
         on = "rel", callers := i.c]
manifest[is.na(callers), callers := vcfs$caller_in_name[match(rel, vcfs$rel)]]
manifest[, path := file.path(ROOT, rel)]
setcolorder(manifest, c("sample", "animal", "tissue", "genotype", "sample_class", "role_in_file", "pairing",
                        "setting", "tumour", "callers", "mutect2_version", "file", "md5", "md5_status"))
manifest[, .N, by = .(setting, sample_class)]
manifest

# ---- 9. Issues ---------------------------------------------------------------
issue <- function(what, detail) if (length(detail) && any(nzchar(detail))) data.table(issue = what, detail = detail[nzchar(detail)])
issues <- rbindlist(list(
  issue("md5 mismatch", files[md5_status == "MISMATCH", rel]),
  issue("file not in any md5 list", files[md5_status == "not listed" & !is_index, rel]),
  issue("listed in md5 file, not on disk", md5_missing[, paste(md5_list, listed_path)]),
  issue("tumour missing from a setting",
        CJ(setting = c("A", "B", "C"), tumour = names(TUMOURS))[!vcfs, on = .(setting, tumour), paste(setting, tumour)]),
  issue("sample column not recognised", unique(manifest[sample_class == "unknown - check", paste(rel, sample)])),
  issue("contig list differs between files",
        if (uniqueN(contig_summary$contig_set_id) > 1) contig_summary[, paste(rel, "-> contig set", contig_set_id)] else character()),
  issue("more than one Mutect2 version", if (uniqueN(gatk_cmd[tool == "Mutect2", version]) > 1)
        gatk_cmd[tool == "Mutect2", paste(rel, version)] else character()),
  issue("setting B: matched KI normal used as a normal", normals_used[setting == "B" & ki_normal_used == TRUE, rel]),
  issue("setting B: file name names the KI normal, but it is not in the file",
        normals_used[setting == "B" & ki_normal_in_name & !ki_normal_used, rel]),
  issue("multiallelic records left", records[multiallelic > 0, paste(rel, multiallelic)]),
  issue("duplicate chrom:pos:ref:alt within a file", dups[, .N, by = rel][, paste(rel, N)]),
  issue("setting A file without INFO/CALLER", has_field[setting == "A" & !INFO_CALLER, rel]),
  issue("small-variant file without FORMAT AD or DP", has_field[setting %chin% c("A", "B", "C") & !(FORMAT_AD & FORMAT_DP), rel]),
  issue("small-variant file without FORMAT AF", has_field[setting %chin% c("A", "B", "C") & !FORMAT_AF, rel]),
  issue("no annotation (ANN/CSQ)", annotation[annotation == "none" & rel %chin% small_files$rel, rel]),
  issue("Mutect2 file without FilterMutectCalls in the header", filter_state[setting %chin% c("B", "C") & !filtermutectcalls, rel]),
  issue("OctopuSV merged file without SOURCES", sv_inventory[setting == "SV_merged" & !has_SOURCES, rel]),
  issue("records with complex REF/ALT (not left-aligned; no FASTA to fix)", records[complex > 0, paste(rel, complex)])
), fill = TRUE)
issues

# ---- 10. Figure and outputs --------------------------------------------------
# Records per file, PASS vs non-PASS, by variant type
plot_dt <- small_dt[, .N, by = .(setting, tumour, TYPE, FILTER_class = fifelse(PASS, "PASS", "non-PASS"))]
p_records <- ggplot(plot_dt, aes(x = tumour, y = N, fill = TYPE)) +
  geom_col() +
  facet_grid(FILTER_class ~ setting, scales = "free_y") +
  labs(x = NULL, y = "records", fill = "type",
       title = paste("Small-variant records per file", "(genome-wide)"),
       subtitle = "A: vs matched KI normal (3 callers merged)  B: vs joint WT (Mutect2)  C: tumour-only + WT PoN (Mutect2)") +
  theme_bw(base_size = 11)
p_records
ggsave(file.path(OUT, "records_per_file.pdf"), p_records, width = 10, height = 6)
ggsave(file.path(OUT, "records_per_file.png"), p_records, width = 10, height = 6, dpi = 150)

# Tables
fwrite(files[, !"path"],      file.path(OUT, "file_inventory.tsv"), sep = "\t")
fwrite(manifest,              file.path(OUT, "manifest.tsv"), sep = "\t")
fwrite(fields,                file.path(OUT, "header_fields.tsv"), sep = "\t")
fwrite(has_field,             file.path(OUT, "key_fields_present.tsv"), sep = "\t")
fwrite(contig_summary,        file.path(OUT, "contigs.tsv"), sep = "\t")
fwrite(meta,                  file.path(OUT, "header_meta_lines.tsv"), sep = "\t")
fwrite(gatk_cmd,              file.path(OUT, "gatk_command_lines.tsv"), sep = "\t")
fwrite(annotation,            file.path(OUT, "annotation.tsv"), sep = "\t")
fwrite(normals_used,          file.path(OUT, "normals_used.tsv"), sep = "\t")
fwrite(records,               file.path(OUT, paste0("records", ".tsv")), sep = "\t")
fwrite(filter_reasons,        file.path(OUT, paste0("filter_reasons", ".tsv")), sep = "\t")
fwrite(filter_state,          file.path(OUT, paste0("filter_state", ".tsv")), sep = "\t")
fwrite(caller_combos,         file.path(OUT, paste0("settingA_caller_combinations", ".tsv")), sep = "\t")
fwrite(dups,                  file.path(OUT, paste0("duplicate_ids", ".tsv")), sep = "\t")
fwrite(sv_counts,             file.path(OUT, paste0("sv_counts", ".tsv")), sep = "\t")
fwrite(sv_inventory,          file.path(OUT, paste0("sv_inventory", ".tsv")), sep = "\t")
fwrite(issues,                file.path(OUT, paste0("issues", ".tsv")), sep = "\t")
fwrite(files[, .(rel, md5, md5_expected, md5_status)], file.path(OUT, "input_md5.tsv"), sep = "\t")

# One RDS with everything (small_dt is the per-record key table Task 1 can start from)
saveRDS(list(files = files, manifest = manifest, headers = headers, sample_cols = sample_cols,
             fields = fields, has_field = has_field, contigs = contigs, meta = meta, gatk_cmd = gatk_cmd,
             annotation = annotation, normals_used = normals_used, wt_names = wt_names,
             small_dt = small_dt, records = records, filter_reasons = filter_reasons,
             caller_combos = caller_combos, dups = dups, sv_dt = sv_dt, sv_counts = sv_counts,
             sv_inventory = sv_inventory, issues = issues, processed_txt = processed_txt),
        file.path(OUT, paste0("inspection", ".rds")))

# Parameters, versions, session
params <- data.table(parameter = c("ROOT", "OUT", "BCFTOOLS", "CHECK_MD5",
                                   paste0("TUMOUR_", names(TUMOURS)), paste0("KI_NORMAL_", names(KI_NORMALS))),
                     value = c(ROOT, OUT, BCFTOOLS, CHECK_MD5, TUMOURS, KI_NORMALS))
fwrite(params, file.path(OUT, "parameters.tsv"), sep = "\t")
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
git_dirty <- tryCatch(length(system2("git", c("-C", shQuote(REPO), "status", "--porcelain", shQuote(SCRIPT)), stdout = TRUE)) > 0,
                      error = function(e) NA)
writeLines(c(R.version.string, bcftools_version,
             paste("data.table", packageVersion("data.table")), paste("ggplot2", packageVersion("ggplot2")),
             paste("git", git_hash, if (isTRUE(git_dirty)) "(script modified)" else "")),
           file.path(OUT, "versions.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

# summary.md: the key facts in a few lines
md_table <- function(x) {
  x <- as.data.table(x)
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
summary_md <- c(
  "# 00_input_inspection summary", "",
  paste0("Run ", format(started, "%Y-%m-%d %H:%M"), "; records from ",
         "the whole genome", "."), "",
  "## Normals actually used (setting B first)", "",
  md_table(normals_used[, .(setting, tumour, sample_columns, normal_sample_lines, normals_mutect2_cmd, pon_mutect2_cmd, ki_normal_used)]), "",
  paste("WT normal names found:", paste(wt_names, collapse = ", ")), "",
  "## Records per file", "", md_table(records[, !"rel"]), "",
  "## Filter state", "", md_table(filter_state[, !"rel"]), "",
  "## Annotation", "", md_table(annotation[, .N, by = .(annotation, ann_version)]), "",
  "## Mutect2 / FilterMutectCalls versions", "", md_table(unique(gatk_cmd[, .(tool, version)])), "",
  "## SV files", "", md_table(sv_inventory[, .(setting, tumour, caller_in_name, records, PASS, has_SOURCES)]), "",
  "## Issues", "", if (nrow(issues)) md_table(issues) else "None found.", "",
  "## Limitations (fixed by the data, not checked here)", "",
  "- No BAMs: no pileups, callable regions, mutations per Mb, SV read support or copy number.",
  "- No mm10 FASTA: indels were split (bcftools norm -m -any) but not left-aligned, so the same indel",
  "  can be written differently by different callers/settings; REF was not checked against the reference.",
  "- SnpEff ANN after splitting lists every allele of the original record; Task 1 must keep the entries for each ALT.",
  "- FreeBayes PASS is the provider's crude filter (MQM>=40, normal 0/0, tumour 0/1, DP>=30 both, QUAL/AO>=5, SAF>0, SAR>0), not a somatic model.")
writeLines(summary_md, file.path(OUT, "summary.md"))

# RUN_LOG.md: newest entry on top, older entries kept
log_entry <- c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - 00_input_inspection"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ", ", bcftools_version[1], ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash, if (isTRUE(git_dirty)) " +dirty" else ""),
  paste0("- Parameters: CHECK_MD5=", CHECK_MD5),
  paste0("- Inputs:     ", ROOT, " (", nrow(vcfs), " variant files; already split with bcftools norm -m -any, not re-normalized)"),
  paste0("- Outputs:    ", OUT),
  paste0("- Counts/QC:  small-variant records ", nrow(small_dt), ", SV records ", nrow(sv_dt), ", issues ", nrow(issues)),
  "")
log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(log_entry, old_log), log_file)

message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
