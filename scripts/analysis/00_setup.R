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
# GRCh38 FASTA sarek aligned to (iGenomes GATK.GRCh38, with its .fai) - for the SNVs'
# trinucleotide context, read with samtools faidx (samtools from this env or wgs-tools).
FASTA <- Sys.getenv("SPRTN_FASTA", file.path("/common/WORK/pstancl/references/igenomes/Homo_sapiens/GATK/GRCh38",
                                             "Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"))
SAMTOOLS <- Filter(nzchar, c(Sys.which("samtools"), Sys.getenv("SPRTN_SAMTOOLS"),
                             if (file.exists("/common/WORK/pstancl/envs/wgs-tools/bin/samtools"))
                               "/common/WORK/pstancl/envs/wgs-tools/bin/samtools"))[1]

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

# Native Strelka2 somatic output has no AF (some wrappers add one - if so, read_vcf_table()
# uses it and this is never called). Tumour VAF from the tier-1 allele-support counts, i.e.
# the first value of each comma pair: SNVs ALT / (REF + ALT) from AU/CU/GU/TU, indels
# TIR / (TAR + TIR).
# Strelka tier-1 allele-support counts for one sample (prefix "t_" tumour, "n_" normal),
# as list(ref, alt): SNVs from AU/CU/GU/TU, indels TAR (ref) / TIR (alt). NULL if absent.
strelka_tier1 <- function(d, prefix = "t_") {
  if (!nrow(d)) return(NULL)
  t1 <- function(x) suppressWarnings(as.numeric(sub(",.*", "", x)))
  base_cols <- paste0(prefix, c("A", "C", "G", "T"), "U")
  if (all(base_cols %in% names(d))) {
    cnt <- vapply(base_cols, function(b) t1(d[[b]]), numeric(nrow(d)))
    if (is.null(dim(cnt))) cnt <- matrix(cnt, nrow = 1)
    colnames(cnt) <- c("A", "C", "G", "T")
    return(list(ref = cnt[cbind(seq_len(nrow(d)), match(d$REF,  colnames(cnt)))],
                alt = cnt[cbind(seq_len(nrow(d)), match(d$ALT1, colnames(cnt)))]))
  }
  if (all(paste0(prefix, c("TAR", "TIR")) %in% names(d)))
    return(list(ref = t1(d[[paste0(prefix, "TAR")]]), alt = t1(d[[paste0(prefix, "TIR")]])))
  NULL
}
strelka_vaf <- function(d, prefix = "t_") {
  k <- strelka_tier1(d, prefix)
  if (is.null(k)) rep(NA_real_, nrow(d)) else k$alt / (k$ref + k$alt)
}

# One row per VCF record: the core columns (VCF names, upper case), then every INFO (info_*), tumour (t_*) and normal
# (n_*) FORMAT field - depth, quality, strand and position metrics for later artefact work.
#   ALT1      first ALT allele. Mutect2 filters sites with >1 ALT as `multiallelic` (never
#             PASS); Strelka and SAGE write one ALT per record - see ALT_COUNT.
#   TYPE      SNV, MNV (same-length multi-base, e.g. SAGE), INDEL, or OTHER (symbolic / *)
#   VAF       tumour VAF: FORMAT/AF (Mutect2, SAGE), or Strelka's tier-1 allele-support counts
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
  out[, MUTTYPE := fcase(TYPE == "INDEL" & nchar(ALT1) > nchar(REF), "INS",
                         TYPE == "INDEL", "DEL", default = TYPE)]
  # VAF / N_VAF: FORMAT/AF of the tumour / normal (Mutect2, SAGE), else Strelka's tier-1 counts,
  # else FORMAT/AD as ALT / (REF + ALT) (MuSE, which has no AF).
  # T_ALT: tumour reads supporting ALT - FORMAT/AD's second value, else Strelka's tier-1 count.
  num1   <- function(x) suppressWarnings(as.numeric(sub(",.*", "", x)))
  second <- function(x) { x <- as.character(x)
    suppressWarnings(as.numeric(fifelse(grepl(",", x, fixed = TRUE), sub("^[^,]*,([^,]*).*$", "\\1", x), NA_character_))) }
  tk <- strelka_tier1(out, "t_")
  ad_vaf <- function(x) { a <- second(x); a / (num1(x) + a) }       # REF,ALT[,...] -> ALT / (REF + ALT)
  af <- function(p) if (paste0(p, "AF") %in% names(out)) num1(out[[paste0(p, "AF")]])
                    else if (!is.null(strelka_tier1(out, p))) strelka_vaf(out, p)
                    else if (paste0(p, "AD") %in% names(out)) ad_vaf(out[[paste0(p, "AD")]])
                    else rep(NA_real_, nrow(out))
  out[, `:=`(VAF   = af("t_"),
             N_VAF = af("n_"),
             T_ALT = if ("t_AD" %in% names(out)) second(out$t_AD) else if (!is.null(tk)) tk$alt else NA_real_)]
  setcolorder(out, c("CHROM", "POS", "REF", "ALT", "ALT1", "ALT_COUNT", "TYPE", "MUTTYPE",
                     "FILTER", "QUAL", "VAF", "N_VAF", "T_ALT"))
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

# Somatic SNV / indel VCFs of every caller that ran, one or more files per caller:
#   sarek main run + extra-callers run (RJALS_vc): mutect2, strelka (snvs + indels), muse,
#   freebayes, lofreq;  oncoanalyser: sage (PURPLE's final VCF).
# The bcftools-normalised copy (scripts/wgs/06_normalize_vcfs.sh) is used when it exists.
# Returns a data.table: CALLER, file, normalised.
find_caller_vcfs <- function() {
  vc_dirs <- file.path(RESULTS, "sarek", paste0(PATIENT, c("", "_vc")), "variant_calling")
  norm_s  <- file.path(SAREK, "normalized_bcftools", PAIR)
  norm_o  <- file.path(RESULTS, "oncoanalyser", PATIENT, "normalized_bcftools")
  known <- data.table(CALLER  = c("mutect2", "strelka", "strelka", "muse", "freebayes", "lofreq"),
                      pattern = paste0("^", PAIR, c("\\.mutect2\\.filtered", "\\.strelka\\.somatic_snvs",
                                                    "\\.strelka\\.somatic_indels", "\\.muse", "\\.freebayes", ".*lofreq.*")))
  hits <- rbindlist(lapply(seq_len(nrow(known)), function(i) {
    raw <- unlist(lapply(vc_dirs[dir.exists(vc_dirs)], list.files, pattern = paste0(known$pattern[i], "\\.vcf\\.gz$"),
                         recursive = TRUE, full.names = TRUE))
    raw <- raw[grepl(paste0("/", PAIR, "/"), raw)]
    raw <- if (length(raw)) raw[1] else NA_character_                # main run first, then RJALS_vc
    nrm <- if (dir.exists(norm_s)) list.files(norm_s, paste0(known$pattern[i], "\\.norm\\.vcf\\.gz$"), full.names = TRUE)[1] else NA
    f <- if (!is.na(nrm)) nrm else raw
    if (is.na(f)) NULL else data.table(CALLER = known$CALLER[i], file = f, normalised = !is.na(nrm))
  }))
  pv <- if (dir.exists(norm_o)) list.files(norm_o, "\\.purple\\.somatic\\.norm\\.vcf\\.gz$", full.names = TRUE)[1] else NA
  pr <- list.files(file.path(ONCO, "purple"), "\\.purple\\.somatic\\.vcf\\.gz$", full.names = TRUE)[1]
  if (!is.na(pv) || !is.na(pr))
    hits <- rbind(hits, data.table(CALLER = "sage", file = if (!is.na(pv)) pv else pr, normalised = !is.na(pv)))
  hits
}

# Read the VCFs in files[keys] into one table tagged with `label` (NULL if none found).
load_calls <- function(files, keys, label) {
  keys <- intersect(keys, names(files))
  if (!length(keys)) return(NULL)
  rbindlist(lapply(files[keys], read_vcf_table), fill = TRUE)[, CALLER := label]   # callers differ in INFO/FORMAT keys
}

# PASS calls with each MNV split into its single-base substitutions (as bcftools norm --atomize
# does), so Mutect2's / SAGE's MNVs can match Strelka's SNVs: Strelka2's somatic output has no
# MNVs, it writes a CC>TT as two C>T. Split rows keep the MNV's metrics and VAF, are flagged
# FROM_MNV, and MNV_KEY names the record they came from. Everything else is unchanged.
atomize_mnv <- function(calls) {
  calls <- copy(calls)[, `:=`(FROM_MNV = FALSE, MNV_KEY = NA_character_)]
  mnv <- calls[MUTTYPE == "MNV"]
  if (!nrow(mnv)) return(calls)
  mnv[, MNV_KEY := paste0(CHROM, ":", POS, ":", REF, ">", ALT1)]
  ex <- mnv[rep(seq_len(nrow(mnv)), nchar(REF))]
  ex[, k := seq_len(.N), by = .(CALLER, MNV_KEY)]
  ex[, `:=`(POS = POS + k - 1L, REF = substr(REF, k, k), ALT1 = substr(ALT1, k, k))]
  ex <- ex[REF != ALT1]
  ex[, `:=`(ALT = ALT1, ALT_COUNT = 1L, TYPE = "SNV", MUTTYPE = "SNV", FROM_MNV = TRUE)][, k := NULL]
  rbindlist(list(calls[MUTTYPE != "MNV"], ex), use.names = TRUE)
}

# Trinucleotide context on the + strand (base before, REF, base after, from the reference)
# for each SNV - turn it into the pyrimidine-strand class with sbs96(). One samtools faidx call for all positions.
# Returns NC_3 in the rows' order; NA if FASTA or samtools are missing.
trinuc_context <- function(d, fasta = FASTA) {
  if (!nrow(d)) return(character(0))
  if (is.na(SAMTOOLS) || !file.exists(fasta) || !file.exists(paste0(fasta, ".fai"))) {
    message("no NC_3: need samtools and ", fasta, " (+ .fai) - set SPRTN_FASTA / SPRTN_SAMTOOLS")
    return(rep(NA_character_, nrow(d)))
  }
  reg <- unique(d[, .(region = paste0(CHROM, ":", POS - 1L, "-", POS + 1L))])
  rf  <- tempfile(fileext = ".txt"); on.exit(unlink(rf))
  writeLines(reg$region, rf)
  fa  <- system2(SAMTOOLS, c("faidx", "-r", shQuote(rf), shQuote(fasta)), stdout = TRUE)
  hdr <- startsWith(fa, ">")                                   # one 3-base sequence line per region
  seqs <- data.table(region = sub("^>", "", fa[hdr]), NC_3 = toupper(fa[which(hdr) + 1L]))
  seqs[d[, .(region = paste0(CHROM, ":", POS - 1L, "-", POS + 1L))], on = "region"]$NC_3
}

# SBS96 class in SigProfiler / COSMIC notation, e.g. A[C>T]G. nc3 is the + strand context from
# trinuc_context(); the class is written on the pyrimidine strand, so for a G or A reference
# the context and the change are reverse-complemented (G>A in CGT -> A[C>T]G), as
# SigProfilerMatrixGenerator and palimpsest do. NA where nc3 is NA.
sbs96 <- function(ref, alt, nc3) {
  rc  <- function(x) chartr("ACGT", "TGCA", x)
  pur <- ref %chin% c("A", "G")
  b5  <- fifelse(pur, rc(substr(nc3, 3, 3)), substr(nc3, 1, 1))
  b3  <- fifelse(pur, rc(substr(nc3, 1, 1)), substr(nc3, 3, 3))
  out <- paste0(b5, "[", fifelse(pur, rc(ref), ref), ">", fifelse(pur, rc(alt), alt), "]", b3)
  out[is.na(nc3)] <- NA_character_
  out
}
SBS96_TYPES <- sort(CJ(b5 = c("A", "C", "G", "T"), sub = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                       b3 = c("A", "C", "G", "T"))[, paste0(b5, "[", sub, "]", b3)], method = "radix")

# PASS SNVs of `atom` (MNVs split, from snv_indel_summary()) get their trinucleotide context
# from the reference FASTA: NC_3 on the + strand, SBS96 the pyrimidine-strand class (A[C>T]G)
# used by the 96-context plots and SigProfiler. Modifies atom in place.
add_sbs96 <- function(atom) {
  atom[MUTTYPE == "SNV", NC_3 := trinuc_context(.SD)]
  bad <- atom[!is.na(NC_3) & substr(NC_3, 2, 2) != REF, .N]      # middle base must be REF
  if (bad) warning(bad, " SNVs whose REF is not the reference base - is FASTA the genome the caller used?")
  atom[MUTTYPE == "SNV" & substr(NC_3, 2, 2) == REF, SBS96 := sbs96(REF, ALT1, NC_3)]
  invisible(atom)
}

# qcVCF 96-context plot of the PASS SNVs in atom (split MNVs included): one column per caller;
# rows all SNVs plus, if `rowsplit` names a column (e.g. QC_SHARED), one row per value.
# plot96_matrix() wants NC_3 on the pyrimidine strand (as palimpsest writes it) and does not
# flip it itself, so it gets SBS96's bases, not the + strand NC_3. It returns the figure(s)
# without saving; orderplots / showperc / dropempty / dontshowall must be single values (their
# defaults are vectors, which its if() checks reject). Needs packages qcVCF does not declare.
# show_all = FALSE: no extra "all SNVs" row (for overlapping groups, where summing them is meaningless)
plot_96context <- function(atom, od, rowsplit = NULL, name = "snv_96context", show_all = TRUE, roworder = NULL) {
  cols <- c("CHROM", "POS", "REF", "ALT1", "CALLER", "SBS96", rowsplit)
  snv96 <- atom[!is.na(SBS96), ..cols]
  setnames(snv96, c("ALT1", "CALLER"), c("ALT", "tool"))
  snv96[, NC_3 := paste0(substr(SBS96, 1, 1), substr(SBS96, 3, 3), substr(SBS96, 7, 7))][, SBS96 := NULL]
  if (!nrow(snv96)) return(invisible(NULL))
  fwrite(snv96, file.path(od, paste0(name, ".csv")))
  miss <- Filter(function(pk) !requireNamespace(pk, quietly = TRUE), c("qcVCF", "cowplot", "stringr", "stringi", "ggtext"))
  if (length(miss)) { message("96-context plot skipped - install: ", paste(miss, collapse = ", ")); return(invisible(NULL)) }
  fig96 <- qcVCF::plot96_matrix(snv96, rowsplit = rowsplit, plotsplitcol = "tool", orderplots = "no", showperc = "yes",
                                dropempty = "no", dontshowall = if (is.null(rowsplit) || !show_all) "yes" else "no",
                                roworder = paste(c("all SNVs", roworder), collapse = ","))
  nrow96 <- if (is.null(rowsplit)) 1 else uniqueN(snv96[[rowsplit]]) + show_all
  for (k in seq_along(fig96))
    save_plot(fig96[[k]], paste0(name, if (k > 1) paste0("_", k)), od,
              w = 9 * uniqueN(snv96$tool), h = 2 + 1.7 * nrow96)
  invisible(snv96)
}

# How many of the n_all callers found a mutation, as a label for the 96-context rows:
# "unique" (1), "2 of 3 callers", ..., "all 3 callers". Unlike a shared/unique split, the
# "all n callers" row is the same set of mutations in every caller's column.
sharing_label <- function(n, n_all) fifelse(n <= 1L, "unique",
  fifelse(n >= n_all, sprintf("all %d callers", n_all), sprintf("%d of %d callers", n, n_all)))
sharing_levels <- function(n_all) c("unique", if (n_all > 2) sprintf("%d of %d callers", 2:(n_all - 1), n_all),
                                   sprintf("all %d callers", n_all))

# Signature sets (rows of atom with a SET column, MNVs split) -> signatures/input/ for the
# SigProfiler script (01b) and FitMS (01c): one minimal VCF per set in vcf/ (one "sample" each;
# the matrix generator rejoins adjacent SNVs into doublets itself), pass_sets.csv,
# pass_set_counts.csv, and the SBS96 counts computed here in SigProfiler's matrix format -
# NOT used for fitting, only to compare with the matrix generator's (sbs96_compare.csv).
write_signature_sets <- function(sets, sig_in, source_label) {
  dir.create(file.path(sig_in, "vcf"), recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(file.path(sig_in, "vcf"), "\\.vcf$", full.names = TRUE))   # no sets left from an earlier run
  sets <- sets[, .(SET, CHROM, POS, REF, ALT = ALT1, MUTTYPE, FROM_MNV, NC_3, SBS96)]
  sets <- unique(sets, by = c("SET", "CHROM", "POS", "REF", "ALT"))
  sets <- sets[order(SET, match(CHROM, STD_CHR), POS)]
  set_counts <- dcast(sets, SET ~ MUTTYPE, fun.aggregate = length, value.var = "POS")
  print(set_counts); fwrite(set_counts, file.path(sig_in, "pass_set_counts.csv"))
  fwrite(sets, file.path(sig_in, "pass_sets.csv"))
  for (st in unique(sets$SET)) {
    f <- file.path(sig_in, "vcf", paste0(st, ".vcf"))
    writeLines(c("##fileformat=VCFv4.2", paste0("##source=", source_label, " PASS set ", st),
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), f)
    fwrite(sets[SET == st, .(CHROM, POS, ID = ".", REF, ALT, QUAL = ".", FILTER = "PASS", INFO = ".")],
           f, sep = "\t", append = TRUE, col.names = FALSE)
  }
  m96 <- dcast(sets[!is.na(SBS96)], SBS96 ~ SET, fun.aggregate = length, value.var = "POS")
  m96 <- m96[data.table(SBS96 = SBS96_TYPES), on = "SBS96"]
  setnames(m96, "SBS96", "MutationType")
  for (cl in setdiff(names(m96), "MutationType")) set(m96, which(is.na(m96[[cl]])), cl, 0L)
  fwrite(m96, file.path(sig_in, paste0(PATIENT, ".SBS96.from_R.txt")), sep = "\t")
  invisible(sets)
}

# What counts as a call. FILTER == "PASS" for every caller except MuSE, whose FILTER is a
# confidence tier, not pass / fail: PASS (highest), Tier1 ... Tier5 (lowest). For WGS, MuSE's
# authors recommend keeping all of them (README, github.com/wwylab/MuSE; WES: all but Tier5).
# SPRTN_MUSE_TIERS=1,2,3,4 (say) keeps only those tiers besides PASS. See CALLERS.md.
MUSE_KEEP <- c("PASS", paste0("Tier", strsplit(Sys.getenv("SPRTN_MUSE_TIERS", "1,2,3,4,5"), ",", fixed = TRUE)[[1]]))
is_pass <- function(d) d$FILTER == "PASS" | (d$CALLER == "muse" & d$FILTER %chin% MUSE_KEEP)

# All records, PASS and filtered, per caller and mutation type (SNV / MNV / INS / DEL): how
# many passed, why the rest failed, and how depth, allele fraction and score compare between
# PASS and FAIL. Writes pass_fail_counts, filter_reasons, metrics_summary and metric_* to od.
raw_qc <- function(calls, od) {
  d <- calls[CHROM %chin% STD_CHR]
  d[, `:=`(STATUS  = factor(fifelse(is_pass(d), "PASS", "FAIL"), levels = c("PASS", "FAIL")),   # MuSE tiers: see is_pass()
           MUTTYPE = factor(MUTTYPE, levels = intersect(c("SNV", "MNV", "INS", "DEL", "OTHER"), unique(MUTTYPE))))]

  # 1. PASS vs FAIL
  counts <- d[, .(total = .N, pass = sum(STATUS == "PASS"), fail = sum(STATUS == "FAIL")), by = .(CALLER, MUTTYPE)]
  counts[, pass_pct := round(100 * pass / total, 1)]
  setorder(counts, CALLER, MUTTYPE)
  print(counts); fwrite(counts, file.path(od, "pass_fail_counts.csv"))
  if ("muse" %in% d$CALLER) {                                   # MuSE's own confidence tiers
    mt <- d[CALLER == "muse", .N, by = .(FILTER, MUTTYPE)][order(match(FILTER, c("PASS", paste0("Tier", 1:5)))), ]
    mt[, kept := FILTER %chin% MUSE_KEEP]
    print(mt); fwrite(mt, file.path(od, "muse_tiers.csv"))
  }
  pf <- d[, .N, by = .(CALLER, MUTTYPE, STATUS)]
  # share of PASS / FAIL per bar (counts differ by orders of magnitude between types), with the counts on it
  save_plot(ggplot(pf, aes(MUTTYPE, N, fill = STATUS)) +
              geom_col(position = "fill") +
              geom_text(aes(label = N), position = position_fill(vjust = 0.5), size = 3) +
              scale_y_continuous(labels = function(x) paste0(100 * x, "%"), expand = c(0, 0)) +
              facet_wrap(~CALLER, scales = "free_x") +
              labs(x = NULL, y = "share of records (numbers = records)", title = "FILTER outcome per caller and mutation type"),
            "pass_fail_counts", od, w = 9, h = 5)

  # 2. why records failed - a record failing several filters counts once for each
  f <- d[STATUS == "FAIL", .(CALLER, MUTTYPE, FILTER)]
  if (nrow(f)) {
    parts   <- strsplit(f$FILTER, ";", fixed = TRUE)
    reasons <- data.table(CALLER  = rep(f$CALLER,  lengths(parts)),
                          MUTTYPE = rep(f$MUTTYPE, lengths(parts)),
                          reason  = unlist(parts, use.names = FALSE))[, .N, by = .(CALLER, MUTTYPE, reason)]
    setorder(reasons, CALLER, -N)
    fwrite(reasons, file.path(od, "filter_reasons.csv"))
    save_plot(ggplot(reasons, aes(reorder(reason, N, sum), N, fill = MUTTYPE)) + geom_col() + coord_flip() +
                facet_wrap(~CALLER, scales = "free") +
                labs(x = NULL, y = "failed records", title = "Why records failed FILTER"),
              "filter_reasons", od, w = 10, h = 6)
  }

  # 3. metrics, PASS vs FAIL. SCORE is each caller's own confidence score - on different
  #    scales, so compare it within a caller only.
  num1  <- function(x) suppressWarnings(as.numeric(sub(",.*", "", x)))
  score <- rep(NA_real_, nrow(d))
  for (col in c("info_TLOD", "info_SomaticEVS")) if (col %in% names(d)) score <- fcoalesce(score, num1(d[[col]]))
  d[, SCORE := fcoalesce(score, QUAL)]
  metrics <- c(t_DP = "tumour depth (DP)", n_DP = "normal depth (DP)", VAF = "tumour VAF",
               N_VAF = "normal VAF", T_ALT = "tumour reads supporting ALT",
               SCORE = "caller score (Mutect2 TLOD / Strelka SomaticEVS / SAGE QUAL)")
  metrics <- metrics[names(metrics) %in% names(d)]
  for (m in names(metrics)) set(d, j = m, value = suppressWarnings(as.numeric(d[[m]])))
  long <- melt(d[, c("CALLER", "MUTTYPE", "STATUS", names(metrics)), with = FALSE],
               id.vars = c("CALLER", "MUTTYPE", "STATUS"), variable.name = "metric",
               value.name = "value", na.rm = TRUE)
  summ <- long[, .(n = .N, median = median(value), q25 = quantile(value, 0.25), q75 = quantile(value, 0.75)),
               by = .(metric, CALLER, MUTTYPE, STATUS)]
  setorder(summ, metric, CALLER, MUTTYPE, STATUS)
  print(summ); fwrite(summ, file.path(od, "metrics_summary.csv"))
  for (m in names(metrics)) {
    x <- long[metric == m]
    if (!nrow(x)) next
    x <- x[, .SD[value <= quantile(value, 0.99)], by = CALLER]      # display only: drop each caller's top 1%
    save_plot(ggplot(x, aes(MUTTYPE, value, fill = STATUS)) + geom_boxplot(outlier.shape = NA) +
                facet_wrap(~CALLER, scales = "free") +
                labs(x = NULL, y = metrics[[m]], title = paste0(metrics[[m]], ": PASS vs FAIL"),
                     subtitle = "box = median and IQR; outliers not drawn; each caller's top 1% left out"),
              paste0("metric_", m), od, w = 9, h = 5)
  }
  invisible(summ)
}

# calls: table from read_vcf_table() with a `CALLER` column. PASS calls only: caller overlap
# (after splitting MNVs, so callers compare like for like), the VAF histogram, the substitution
# spectrum (split MNVs included), and the PASS tables, written to od. Returns, invisibly,
# list(pass = PASS records as called, atom = the same with MNVs split - see atomize_mnv()).
snv_indel_summary <- function(calls, od) {
  calls_pass <- calls[CHROM %chin% STD_CHR & is_pass(calls)]          # MuSE: PASS + Tier1-5 (see is_pass)
  atom       <- atomize_mnv(calls_pass)

  if (uniqueN(atom$CALLER) > 1) {
    # NB indels can be written differently by different callers (normalise with
    # bcftools norm before trusting the indel overlap).
    concord <- unique(atom[, .(CALLER, TYPE, key = paste(CHROM, POS, REF, ALT1, sep = ":"))])
    concord <- concord[, .(callers = paste(sort(CALLER), collapse = "+")), by = .(key, TYPE)]
    concord <- concord[, .(n = .N), by = .(TYPE, callers)]
    print(concord); fwrite(concord, file.path(od, "snv_indel_concordance.csv"))
    save_plot(ggplot(concord, aes(reorder(callers, n), n, fill = TYPE)) +
                geom_col(position = position_dodge(width = 0.9, preserve = "single")) +
                geom_text(aes(label = n), position = position_dodge(width = 0.9, preserve = "single"), hjust = -0.15, size = 3.5) +
                scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +      # room for the labels
                coord_flip() + labs(x = NULL, y = "PASS calls (MNVs split into SNVs)", title = "Overlap between callers"),
              "snv_indel_concordance", od)
  }

  vaf_tbl <- calls_pass[!is.na(VAF) & TYPE == "SNV"]
  if (nrow(vaf_tbl)) save_plot(ggplot(vaf_tbl, aes(VAF)) + geom_histogram(bins = 50) + facet_wrap(~CALLER) +
                                 labs(x = "tumour VAF", title = "PASS SNV allele fractions"), "snv_vaf_hist", od)

  # substitution spectrum: 6 classes, pyrimidine reference; split MNVs included
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  snv  <- atom[TYPE == "SNV" & REF %chin% names(comp) & ALT1 %chin% names(comp)]
  snv[, pyr := REF %chin% c("C", "T")]
  snv[, class := paste0(fifelse(pyr, REF,  unname(comp[REF])), ">",
                        fifelse(pyr, ALT1, unname(comp[ALT1])))]
  spec <- snv[, .(n = .N), by = .(CALLER, class)]
  fwrite(spec, file.path(od, "snv_substitution_spectrum.csv"))
  save_plot(ggplot(spec, aes(class, n, fill = CALLER)) + geom_col(position = position_dodge(preserve = "single")) +
              labs(x = NULL, y = "PASS SNVs (incl. split MNVs)", title = "Substitution spectrum"), "snv_spectrum", od)

  fwrite(calls_pass, file.path(od, "snv_indel_pass.csv"))
  fwrite(atom,       file.path(od, "snv_indel_pass_atomized.csv"))
  invisible(list(pass = calls_pass, atom = atom))
}

# sv: table from read_sv_vcf() with a `CALLER` column.
sv_summary <- function(sv, od) {
  sv <- copy(sv)[CHROM %chin% STD_CHR]
  sv[, SVTYPE := fcoalesce(SVTYPE, "unknown")]
  sv_counts <- sv[, .(total = .N, pass = sum(FILTER == "PASS")), by = .(CALLER, SVTYPE)]
  print(sv_counts); fwrite(sv_counts, file.path(od, "sv_counts.csv"))
  save_plot(ggplot(sv[FILTER == "PASS"], aes(SVTYPE, fill = CALLER)) + geom_bar(position = position_dodge(preserve = "single")) +
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
