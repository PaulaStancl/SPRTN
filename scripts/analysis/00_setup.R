# ---------------------------------------------------------------------------
# 00_setup.R - shared setup for the RJALS mutation analysis. Sourced by 01 / 02.
#
# Patient data: the results live on the server and must NOT end up in OneDrive.
# Run on the server (Rscript / RStudio there), or point SPRTN_RESULTS at a
# private local copy. Environment variables override the defaults:
#   SPRTN_RESULTS   results/wgs folder     SPRTN_OUT   where tables/plots go
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

# ---- paths ------------------------------------------------------------------
RESULTS <- Sys.getenv("SPRTN_RESULTS", "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/results/wgs")
OUT     <- Sys.getenv("SPRTN_OUT", file.path(RESULTS, "analysis"))
for (p in c(RESULTS, OUT)) {
  if (grepl("OneDrive|CloudStorage|Dropbox|iCloud|Google Drive", normalizePath(p, mustWork = FALSE), ignore.case = TRUE))
    stop("Refusing to use a cloud-synced folder for patient data: ", p, call. = FALSE)
}

PATIENT <- "RJALS"
TUMOUR  <- "RJALS_Tm"
NORMAL  <- "RJALS_N"
PAIR    <- paste0(TUMOUR, "_vs_", NORMAL)           # sarek's name for the paired calls

SAREK <- file.path(RESULTS, "sarek", PATIENT)
ONCO  <- file.path(RESULTS, "oncoanalyser", PATIENT, PATIENT)   # .../oncoanalyser/RJALS/RJALS/purple etc.
TEVO  <- file.path(RESULTS, "tumourevo", PATIENT)

# ---- genome (GRCh38 primary assembly) ----------------------------------------
CHR_LEN <- c(chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
             chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
             chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
             chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
             chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
             chr21 = 46709983, chr22 = 50818468, chrX = 156040895, chrY = 57227415)
STD_CHR    <- names(CHR_LEN)
CHR_OFFSET <- setNames(c(0, head(cumsum(CHR_LEN), -1)), STD_CHR)   # for genome-wide plots

# ---- qcVCF -------------------------------------------------------------------
HAVE_QCVCF <- requireNamespace("qcVCF", quietly = TRUE)
if (!HAVE_QCVCF) message("qcVCF is not installed - the qcVCF sections will be skipped")

theme_set(theme_bw(base_size = 11))

# ---- helpers -----------------------------------------------------------------
# Exactly one file under `dir` matching `pattern`; stops with the candidates if not.
find_one <- function(dir, pattern, required = TRUE) {
  hits <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 1) return(hits)
  if (length(hits) == 0 && !required) { message("not found (skipped): ", pattern, " under ", dir); return(NA_character_) }
  stop(sprintf("Expected 1 file matching '%s' under %s, found %d%s", pattern, dir, length(hits),
               if (length(hits)) paste0(":\n  ", paste(hits, collapse = "\n  ")) else ""), call. = FALSE)
}

need_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(sprintf("%s: missing column(s) %s. Columns are: %s", what,
                                 paste(miss, collapse = ", "), paste(names(df), collapse = ", ")), call. = FALSE)
}

norm_chr <- function(x) { x <- as.character(x); fifelse(grepl("^chr", x), x, paste0("chr", x)) }

save_plot <- function(p, name, dir, w = 7, h = 4) ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h)

# VCF data lines as a data.table via fread - far faster than VariantAnnotation::readVcf;
# INFO / FORMAT are split into columns afterwards by vcf_metrics(). .gz is decompressed
# with gzip (fread's own .gz support needs R.utils).
read_vcf_dt <- function(path) {
  src <- if (grepl("\\.gz$", path)) list(cmd = paste("gzip -dc", shQuote(path))) else list(file = path)
  v <- do.call(fread, c(src, list(skip = "#CHROM", sep = "\t", quote = "", showProgress = FALSE,
                                  colClasses = list(character = c(1, 3, 4, 5, 7)))))
  setnames(v, 1, "CHROM")
  v
}

# Value of FORMAT field `key` in one sample column, per record (NA where absent).
format_field <- function(fmt, values, key) {
  out <- rep(NA_character_, length(fmt))
  for (f in unique(fmt)) {
    k <- match(key, strsplit(f, ":", fixed = TRUE)[[1]])
    if (is.na(k)) next
    i <- which(fmt == f)
    out[i] <- tstrsplit(values[i], ":", fixed = TRUE, fill = NA_character_, keep = k)[[1]]
  }
  out
}

# Column holding a sample: its name if present (Mutect2 RJALS_RJALS_Tm, SAGE RJALS_Tm),
# else the generic name Strelka uses (TUMOR / NORMAL). NA if neither is there.
sample_col <- function(v, id, generic) {
  hit <- grep(id, names(v), fixed = TRUE, value = TRUE)
  if (length(hit) == 1) return(hit)
  if (generic %in% names(v)) return(generic)
  NA_character_
}

# Numbers where a column is all-numeric, TRUE for flags; "12,30"-style lists stay text.
convert_types <- function(d) {
  for (j in names(d)) set(d, j = j, value = type.convert(d[[j]], as.is = TRUE, na.strings = c(".", "NA")))
  d
}

# Every INFO key as its own column, prefixed info_ (e.g. info_DP, info_TLOD, info_SomaticEVS,
# info_TIER). Flags without a value (SOMATIC) become TRUE. NULL if the VCF has no INFO.
info_fields <- function(info) {
  n     <- length(info)
  parts <- strsplit(info, ";", fixed = TRUE)
  l <- data.table(i = rep(seq_len(n), lengths(parts)), kv = unlist(parts, use.names = FALSE))[kv != "."]
  if (!nrow(l)) return(NULL)
  l[, `:=`(key   = sub("=.*", "", kv),
           value = fifelse(grepl("=", kv, fixed = TRUE), sub("^[^=]*=", "", kv), "TRUE"))]
  l <- unique(l, by = c("i", "key"))                              # a repeated key: keep the first
  w <- dcast(l, i ~ key, value.var = "value")                     # no fun.aggregate: ~30x faster
  w <- w[data.table(i = seq_len(n)), on = "i"][, i := NULL]      # keep records with no INFO
  setnames(w, paste0("info_", names(w)))
  convert_types(w)
}

# Every FORMAT key of one sample column as its own column, prefixed (t_DP, t_AD, n_DP ...).
sample_fields <- function(fmt, values, prefix) {
  keys <- unique(unlist(strsplit(unique(fmt), ":", fixed = TRUE)))
  w <- as.data.table(setNames(lapply(keys, function(k) format_field(fmt, values, k)), paste0(prefix, keys)))
  convert_types(w)
}

# INFO + tumour (t_) + normal (n_) FORMAT fields for every record of a VCF read by read_vcf_dt().
vcf_metrics <- function(v) {
  tcol <- sample_col(v, TUMOUR, "TUMOR"); ncl <- sample_col(v, NORMAL, "NORMAL")
  fmt  <- "FORMAT" %in% names(v)
  parts <- list(
    if ("INFO" %in% names(v)) info_fields(v$INFO),
    if (fmt && !is.na(tcol)) sample_fields(v$FORMAT, v[[tcol]], "t_"),
    if (fmt && !is.na(ncl))  sample_fields(v$FORMAT, v[[ncl]],  "n_"))
  parts <- Filter(Negate(is.null), parts)
  if (length(parts)) do.call(cbind, parts) else NULL
}

# Strelka writes no AF. Tumour VAF from its tier-1 counts, as its documentation recommends:
# SNVs alt / (ref + alt) from AU/CU/GU/TU, indels TIR / (TAR + TIR).
strelka_vaf <- function(d) {
  if (!nrow(d)) return(numeric(0))
  t1 <- function(x) suppressWarnings(as.numeric(sub(",.*", "", x)))
  if (all(paste0("t_", c("A", "C", "G", "T"), "U") %in% names(d))) {
    cnt <- vapply(c("A", "C", "G", "T"), function(b) t1(d[[paste0("t_", b, "U")]]), numeric(nrow(d)))
    if (is.null(dim(cnt))) cnt <- matrix(cnt, nrow = 1, dimnames = list(NULL, c("A", "C", "G", "T")))
    r <- cnt[cbind(seq_len(nrow(d)), match(d$REF,  colnames(cnt)))]
    a <- cnt[cbind(seq_len(nrow(d)), match(d$ALT1, colnames(cnt)))]
    return(a / (r + a))
  }
  if (all(c("t_TAR", "t_TIR") %in% names(d))) { r <- t1(d$t_TAR); a <- t1(d$t_TIR); return(a / (r + a)) }
  rep(NA_real_, nrow(d))
}

# One row per VCF record: the core columns (VCF names, upper case), then every INFO (info_*), tumour (t_*) and normal
# (n_*) FORMAT field - depth, quality, strand and position metrics for later artefact work.
#   ALT1      first ALT allele. Mutect2 filters sites with >1 ALT as `multiallelic` (never
#             PASS); Strelka and SAGE write one ALT per record - see ALT_COUNT.
#   TYPE      SNV, MNV (same-length multi-base, e.g. SAGE), INDEL, or OTHER (symbolic / *)
#   VAF       tumour VAF: FORMAT/AF (Mutect2, SAGE), or from Strelka's read counts
read_vcf_table <- function(path) {
  v   <- read_vcf_dt(path)
  out <- v[, .(CHROM, POS, REF, ALT, ALT1 = sub(",.*", "", ALT),
               ALT_COUNT = lengths(strsplit(ALT, ",", fixed = TRUE)),
               FILTER, QUAL = suppressWarnings(as.numeric(QUAL)))]
  out[, TYPE := fcase(grepl("^[<*.]", ALT1),               "OTHER",
                      nchar(REF) == 1 & nchar(ALT1) == 1, "SNV",
                      nchar(REF) == nchar(ALT1),          "MNV",
                      default = "INDEL")]
  m <- vcf_metrics(v)
  if (!is.null(m)) out <- cbind(out, m)
  tvaf <- if ("t_AF" %in% names(out)) suppressWarnings(as.numeric(sub(",.*", "", out$t_AF))) else strelka_vaf(out)
  out[, VAF := tvaf]
  setcolorder(out, c("CHROM", "POS", "REF", "ALT", "ALT1", "ALT_COUNT", "TYPE", "FILTER", "QUAL", "VAF"))
  out
}

# SV VCF (Manta, ESVEE): one row per record plus all INFO / FORMAT metrics, as above.
# NB a translocation / inversion is two BND records (the two breakends), so BND counts
# are ~2x the number of events.
read_sv_vcf <- function(path) {
  v   <- read_vcf_dt(path)
  out <- v[, .(ID, CHROM, POS, FILTER,
               SVTYPE = fifelse(grepl("(^|;)SVTYPE=", INFO), sub(".*(^|;)SVTYPE=([^;]+).*", "\\2", INFO), NA_character_))]
  m <- vcf_metrics(v)
  if (!is.null(m)) out <- cbind(out, m)
  out
}

# ---------------------------------------------------------------------------
# Shared analysis steps - used by 01_sarek.R and 02_oncoanalyser.R
# ---------------------------------------------------------------------------

# Read the VCFs in files[keys] into one table tagged with `label` (NULL if none found).
load_calls <- function(files, keys, label) {
  keys <- intersect(keys, names(files))
  if (!length(keys)) return(NULL)
  rbindlist(lapply(files[keys], read_vcf_table), fill = TRUE)[, CALLER := label]   # callers differ in INFO/FORMAT keys
}

# calls: table from read_vcf_table() with a `CALLER` column. Writes counts, caller overlap
# (when there is more than one caller), the VAF histogram, the substitution spectrum and
# the PASS calls into od. Returns the PASS calls.
snv_indel_summary <- function(calls, od) {
  calls      <- calls[CHROM %chin% STD_CHR]
  calls_pass <- calls[FILTER == "PASS"]

  counts <- calls[, .(total = .N, pass = sum(FILTER == "PASS")), by = .(CALLER, TYPE)]
  print(counts); fwrite(counts, file.path(od, "snv_indel_counts.csv"))

  if (uniqueN(calls_pass$CALLER) > 1) {
    # NB indels can be written differently by different callers (normalise with
    # bcftools norm before trusting the indel overlap).
    concord <- unique(calls_pass[, .(CALLER, TYPE, key = paste(CHROM, POS, REF, ALT1, sep = ":"))])
    concord <- concord[, .(callers = paste(sort(CALLER), collapse = "+")), by = .(key, TYPE)]
    concord <- concord[, .(n = .N), by = .(TYPE, callers)]
    print(concord); fwrite(concord, file.path(od, "snv_indel_concordance.csv"))
    save_plot(ggplot(concord, aes(reorder(callers, n), n, fill = TYPE)) +
                geom_col(position = position_dodge(width = 0.9)) +
                geom_text(aes(label = n), position = position_dodge(width = 0.9), hjust = -0.15, size = 3.5) +
                scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +      # room for the labels
                coord_flip() + labs(x = NULL, y = "PASS calls", title = "Overlap between callers"),
              "snv_indel_concordance", od)
  }

  vaf_tbl <- calls_pass[!is.na(VAF) & TYPE == "SNV"]
  if (nrow(vaf_tbl)) save_plot(ggplot(vaf_tbl, aes(VAF)) + geom_histogram(bins = 50) + facet_wrap(~CALLER) +
                                 labs(x = "tumour VAF", title = "PASS SNV allele fractions"), "snv_vaf_hist", od)

  # substitution spectrum: 6 classes, pyrimidine reference
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  snv  <- calls_pass[TYPE == "SNV" & REF %chin% names(comp) & ALT1 %chin% names(comp)]
  snv[, pyr := REF %chin% c("C", "T")]
  snv[, class := paste0(fifelse(pyr, REF,  unname(comp[REF])), ">",
                        fifelse(pyr, ALT1, unname(comp[ALT1])))]
  spec <- snv[, .(n = .N), by = .(CALLER, class)]
  fwrite(spec, file.path(od, "snv_substitution_spectrum.csv"))
  save_plot(ggplot(spec, aes(class, n, fill = CALLER)) + geom_col(position = "dodge") +
              labs(x = NULL, y = "PASS SNVs", title = "Substitution spectrum"), "snv_spectrum", od)

  fwrite(calls_pass, file.path(od, "snv_indel_pass.csv"))
  invisible(calls_pass)
}

# sv: table from read_sv_vcf() with a `CALLER` column.
sv_summary <- function(sv, od) {
  sv <- copy(sv)[CHROM %chin% STD_CHR]
  sv[, SVTYPE := fcoalesce(SVTYPE, "unknown")]
  sv_counts <- sv[, .(total = .N, pass = sum(FILTER == "PASS")), by = .(CALLER, SVTYPE)]
  print(sv_counts); fwrite(sv_counts, file.path(od, "sv_counts.csv"))
  save_plot(ggplot(sv[FILTER == "PASS"], aes(SVTYPE, fill = CALLER)) + geom_bar(position = "dodge") +
              labs(x = NULL, y = "PASS records", title = "SVs by type (BND = two records per event)"), "sv_types", od)
  invisible(sv)
}

# Copy-number segments as one table: source, chr, start, end, cn (total), minor.
ascat_segments <- function(path) {
  asc <- fread(path)
  need_cols(asc, c("chr", "startpos", "endpos", "nMajor", "nMinor"), "ASCAT segments")
  asc[, .(source = "ASCAT", chr = norm_chr(chr), start = startpos, end = endpos,
          cn = nMajor + nMinor, minor = nMinor)]
}
purple_segments <- function(path) {
  pur <- fread(path)
  need_cols(pur, c("chromosome", "start", "end", "copyNumber", "minorAlleleCopyNumber"), "PURPLE cnv")
  pur[, .(source = "PURPLE", chr = norm_chr(chromosome), start = start, end = end,
          cn = copyNumber, minor = minorAlleleCopyNumber)]
}

# Genome-wide copy-number plot, fraction of the autosomal genome per copy-number state, LOH.
cn_summary <- function(seg, od) {
  seg <- seg[chr %chin% STD_CHR]
  seg[, `:=`(x0 = start + unname(CHR_OFFSET[chr]), x1 = end + unname(CHR_OFFSET[chr]))]
  fwrite(seg, file.path(od, "cn_segments.csv"))
  save_plot(
    ggplot(seg, aes(x = x0, xend = x1, y = cn, yend = cn)) +
      geom_vline(xintercept = CHR_OFFSET, colour = "grey85", linewidth = 0.2) +
      geom_segment(linewidth = 1) + facet_grid(source ~ .) +
      scale_x_continuous(breaks = CHR_OFFSET + CHR_LEN / 2, labels = sub("chr", "", STD_CHR), expand = c(0, 0)) +
      coord_cartesian(ylim = c(0, 8)) + labs(x = "chromosome", y = "total copy number", title = "Copy number"),
    "cn_genome", od, w = 11, h = 5)

  auto <- seg[!chr %chin% c("chrX", "chrY")]
  auto[, `:=`(w = end - start + 1, state = pmin(round(cn), 6), loh = round(minor) == 0)]
  state_frac <- auto[, .(w = sum(w)), by = .(source, state)][, fraction := w / sum(w), by = source][, w := NULL]
  cn_frac    <- auto[, .(LOH = sum(w[which(loh)]) / sum(w)), by = source]
  print(state_frac); print(cn_frac)
  fwrite(state_frac, file.path(od, "cn_state_fraction.csv"))
  fwrite(cn_frac,    file.path(od, "cn_loh_fraction.csv"))
  invisible(seg)
}
