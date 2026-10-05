# ---------------------------------------------------------------------------
# 00_setup.R - shared setup for the RJALS mutation analysis. Sourced by 01 / 02.
#
# Patient data: the results live on the server and must NOT end up in OneDrive.
# Run on the server (Rscript / RStudio there), or point SPRTN_RESULTS at a
# private local copy. Environment variables override the defaults:
#   SPRTN_RESULTS   results/wgs folder     SPRTN_OUT   where tables/plots go
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(VariantAnnotation); library(GenomicRanges)
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

# ---- qcVCR -------------------------------------------------------------------
HAVE_QCVCR <- requireNamespace("qcVCR", quietly = TRUE)
if (!HAVE_QCVCR) message("qcVCR is not installed - the qcVCR sections will be skipped")

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

norm_chr <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))

save_plot <- function(p, name, dir, w = 7, h = 4) ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h)

# Tumour allele fraction from FORMAT/AF (Mutect2, SAGE/PURPLE); NA when absent (Strelka).
tumour_vaf <- function(vcf) {
  na <- rep(NA_real_, nrow(vcf))
  tryCatch({
    g <- geno(vcf)
    if (!"AF" %in% names(g)) return(na)
    col <- grep(TUMOUR, colnames(vcf), fixed = TRUE)
    if (length(col) != 1) return(na)
    m <- g[["AF"]]
    x <- if (length(dim(m)) == 3) m[, col, 1] else m[, col]
    if (is.list(x) || methods::is(x, "List")) x <- vapply(x, function(v) as.numeric(v)[1], numeric(1))
    as.numeric(x)
  }, error = function(e) na)
}

# One row per VCF record: chr, pos, ref, alt (first allele in alt1), FILTER, SNV/INDEL, VAF.
read_vcf_table <- function(path) {
  vcf  <- readVcf(path)
  rr   <- rowRanges(vcf)
  ref  <- as.character(rr$REF)
  alt  <- vapply(rr$ALT, function(a) paste(as.character(a), collapse = ","), character(1))
  alt1 <- sub(",.*", "", alt)
  tibble(chr = as.character(seqnames(rr)), pos = start(rr), ref = ref, alt = alt, alt1 = alt1,
         filter = as.character(rr$FILTER),
         type = if_else(nchar(ref) == 1 & nchar(alt1) == 1, "SNV", "INDEL"),
         vaf = tumour_vaf(vcf))
}

# SV VCF (Manta, ESVEE): one row per record. NB a translocation / inversion is two BND
# records (the two breakends), so BND counts are ~2x the number of events.
read_sv_vcf <- function(path) {
  vcf <- readVcf(path); rr <- rowRanges(vcf); inf <- info(vcf)
  tibble(id = names(rr), chr = as.character(seqnames(rr)), pos = start(rr),
         filter = as.character(rr$FILTER),
         svtype = if ("SVTYPE" %in% names(inf)) as.character(inf$SVTYPE) else NA_character_)
}
