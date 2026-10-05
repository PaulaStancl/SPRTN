# ---------------------------------------------------------------------------
# 01_snv_indel.R - somatic SNVs and indels, RJALS tumour vs normal.
# Callers: Mutect2 + Strelka2 (sarek) and SAGE/PURPLE (oncoanalyser).
# Run:  Rscript 01_snv_indel.R   (or line by line in R - any of the folders above works)
# ---------------------------------------------------------------------------
# Finds 00_setup.R from this folder, from scripts/, or from the project root.
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
od <- file.path(OUT, "snv_indel"); dir.create(od, recursive = TRUE, showWarnings = FALSE)

# ---- 1. input files ---------------------------------------------------------
files <- c(
  mutect2       = find_one(file.path(SAREK, "variant_calling/mutect2", PAIR), "\\.mutect2\\.filtered\\.vcf\\.gz$"),
  strelka_snv   = find_one(file.path(SAREK, "variant_calling/strelka", PAIR), "somatic_snvs\\.vcf\\.gz$",   required = FALSE),
  strelka_indel = find_one(file.path(SAREK, "variant_calling/strelka", PAIR), "somatic_indels\\.vcf\\.gz$", required = FALSE),
  purple        = find_one(file.path(ONCO, "purple"), "\\.purple\\.somatic\\.vcf\\.gz$", required = FALSE)
)
files <- files[!is.na(files)]
writeLines(paste(names(files), files, sep = "\t"), file.path(od, "inputs_used.tsv"))

# ---- 2. load into one table (one row per record, tagged with the caller) -----
load_caller <- function(keys, label) {
  keys <- intersect(keys, names(files))
  if (!length(keys)) return(NULL)
  rbindlist(lapply(files[keys], read_vcf_table))[, caller := label]
}
calls <- rbindlist(list(                      # rbindlist() skips NULL (a caller that was not found)
  load_caller("mutect2", "mutect2"),
  load_caller(c("strelka_snv", "strelka_indel"), "strelka"),
  load_caller("purple", "sage_purple")
))[chr %in% STD_CHR]
calls_pass <- calls[filter == "PASS"]

# ---- 3. counts per caller ----------------------------------------------------
counts <- calls[, .(total = .N, pass = sum(filter == "PASS")), by = .(caller, type)]
print(counts); fwrite(counts, file.path(od, "counts_per_caller.csv"))

# ---- 4. caller concordance (PASS calls) --------------------------------------
# NB indels can be written differently by different callers (normalise with
# bcftools norm before trusting the indel overlap).
concord <- unique(calls_pass[, .(caller, type, key = paste(chr, pos, ref, alt1, sep = ":"))])
concord <- concord[, .(callers = paste(sort(caller), collapse = "+")), by = .(key, type)]
concord <- concord[, .(n = .N), by = .(type, callers)]
print(concord); fwrite(concord, file.path(od, "caller_concordance.csv"))
save_plot(ggplot(concord, aes(reorder(callers, n), n, fill = type)) + geom_col(position = "dodge") +
            coord_flip() + labs(x = NULL, y = "PASS calls", title = "Overlap between callers"),
          "caller_concordance", od)

# ---- 5. tumour allele-fraction distribution ----------------------------------
vaf_tbl <- calls_pass[!is.na(vaf) & type == "SNV"]
if (nrow(vaf_tbl)) save_plot(ggplot(vaf_tbl, aes(vaf)) + geom_histogram(bins = 50) + facet_wrap(~caller) +
                               labs(x = "tumour VAF", title = "PASS SNV allele fractions"), "vaf_hist", od)

# ---- 6. substitution spectrum (6 classes, pyrimidine reference) ---------------
comp <- c(A = "T", C = "G", G = "C", T = "A")
snv <- calls_pass[type == "SNV" & ref %chin% names(comp) & alt1 %chin% names(comp)]
snv[, pyr := ref %chin% c("C", "T")]
snv[, class := paste0(fifelse(pyr, ref,  unname(comp[ref])), ">",
                      fifelse(pyr, alt1, unname(comp[alt1])))]
spec <- snv[, .(n = .N), by = .(caller, class)]
fwrite(spec, file.path(od, "substitution_spectrum.csv"))
save_plot(ggplot(spec, aes(class, n, fill = caller)) + geom_col(position = "dodge") +
            labs(x = NULL, y = "PASS SNVs", title = "Substitution spectrum"), "spectrum", od)

# ---- 7. drivers (tumourevo) ---------------------------------------------------
# TODO: readRDS() the tumourevo driver object
#   <TEVO>/driver_annotation/annotate_driver/RJALS/RJALS/RJALS_RJALS_Tm/*_driver.rds
# and tabulate is_driver / driver_label. Look at str(x, max.level = 2) first.

# ---- 8. qcVCR -----------------------------------------------------------------
# TODO(Paula): call your qcVCR functions here and tell me what they take (a VCF
# path? a data frame?) so I can wire them in.
if (HAVE_QCVCR) {
  library(qcVCR)
  # qc <- qcVCR::<function>(files[["mutect2"]])
}

fwrite(calls_pass, file.path(od, "calls_pass.csv"))
message("done: ", od)
