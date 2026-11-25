#!/usr/bin/env Rscript
# Author: Paula Breyer
# Project: SPRTN (Maria Boskovic), mouse WGS | Created: 2026-11-24
#
# mouse_wgs_snv_indel_callers.R - FIRST STEP: callers in tumour vs matched KI normal (setting A).
#
# Folder: data/mouse/wgs/tumor_vs_normal/snv_indel/
# Callers: Strelka2, Mutect2, FreeBayes. The provider split multiallelics (bcftools norm -m -any)
# and merged the callers into one file per tumour with INFO/CALLER (the "union").
# PASS as set by the provider:
#   strelka   nf-core/sarek defaults
#   mutect2   FilterMutectCalls with sarek's orientation-bias priors
#   freebayes crude prefilter: MQM>=40, normal 0/0, tumour 0/1, DP>=30 in both, QUAL/AO>=5, SAF>0, SAR>0
#
# MERGED FILE FILTER (checked 2026-11-24, all 6.8 M merged records, no exception): it is the bcftools
# merge default - PASS only if >=1 caller is PASS and no caller that found the variant has a named
# filter; FreeBayes writes "." when its prefilter fails, so it never blocks PASS. INFO/CALLER lists every
# caller that CALLED the variant, passed or not. So "merged PASS + >=2 callers" is NOT ">=2 callers PASS".
# The correct per-caller PASS sets and intersects come from the single-caller files (R3, "my intersect").
#
# Questions answered, per tumour:
#   1. how many SNVs and indels each caller finds
#   2. how many are shared between each pair of callers (and by 1, 2 or 3 callers)
#   3. how the merged union compares with the single callers
#   4. the SBS96 (trinucleotide context) profile of each caller, the union, >=2 callers and all 3
#
# Blocks:
#   0  EDIT THIS BLOCK
#   1  setup and helpers
#   2  files: per-caller files and/or the merged file
#   3  read PASS records; per-caller sets (from per-caller files if present, else from INFO/CALLER)
#   4  variant types; MNVs split into SNVs for the comparison
#   5  counts per caller and the union (+ qcVCF plot_mutation_counts)
#   6  pairwise sharing between callers (+ qcVCF plot_pairwise_shared_mutations), exclusive categories
#   7  union vs single callers
#   8  SBS96 context: profiles, plots (ggplot + qcVCF plot96_matrix), cosine similarity
#   9  outputs, summary.md, RUN_LOG.md
#
# Comparison unit: chrom:pos:ref:alt. MNVs are split into SNVs (Strelka writes a CC>TT as two C>T,
# Mutect2/FreeBayes as one MNV), as in scripts/analysis/01_sarek.R. Doublet halves (adjacent SNVs
# within one set) are left out of SBS96, as in scripts/analysis/00_setup.R.
# Open question (memory, 2026-10-07): split MNVs vs join doublets as the comparison unit - not changed here.
# Indels are not left-aligned (no -f in the provider's bcftools norm), so FreeBayes indels written
# differently from Mutect2/Strelka will not match even when they are the same event.
#
# Env on ssi: r-mutation (data.table, ggplot2, qcVCF, bcftools, samtools).
# Output: results/mouse_wgs/mouse_wgs_snv_indel_callers/

# ---- 0. EDIT THIS BLOCK ------------------------------------------------------
REPO     <- "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN"
DIR_A    <- file.path(REPO, "data/mouse/wgs/tumor_vs_normal/snv_indel")
OUT      <- file.path(REPO, "results/mouse_wgs/mouse_wgs_snv_indel_callers")
SCRIPT   <- file.path(REPO, "scripts/mouse_wgs/mouse_wgs_snv_indel_callers.R")

# GRCm38 FASTA (Ensembl 102 primary assembly, from scripts/cross_species/h2m 00_resources.R), + .fai
FASTA    <- "/common/WORK/pstancl/references/cross_species/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa"

# Tools: RStudio often does not see the env's PATH, so fall back to the env's bin/
ENV_BIN  <- "/common/WORK/pstancl/envs/r-mutation/bin"
BCFTOOLS <- Sys.which("bcftools"); if (!nzchar(BCFTOOLS)) BCFTOOLS <- file.path(ENV_BIN, "bcftools")
SAMTOOLS <- Sys.which("samtools"); if (!nzchar(SAMTOOLS)) SAMTOOLS <- file.path(ENV_BIN, "samtools")

TUMOURS <- c(HCC_1 = "sfb_dikic_2025_05_10_KI_HCC_1",
             HCC_2 = "sfb_dikic_2025_05_11_KI_HCC_2",
             HCC_3 = "sfb_dikic_2025_05_12_KI_HCC_3")
CALLERS <- c("strelka", "mutect2", "freebayes")


# ---- 1. Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# Console as wide as the terminal (falls back to 80 when there is no terminal, e.g. batch Rscript)
options(width = max(80, suppressWarnings(as.numeric(system("tput cols", intern = TRUE, ignore.stderr = TRUE)))[1], na.rm = TRUE))
HAVE_QCVCF <- requireNamespace("qcVCF", quietly = TRUE)
if (!file.exists(BCFTOOLS)) stop("bcftools not found at ", BCFTOOLS, " - install it in r-mutation or set BCFTOOLS")
HAVE_FASTA <- file.exists(SAMTOOLS) && file.exists(FASTA) && file.exists(paste0(FASTA, ".fai"))
if (!HAVE_FASTA) message("No samtools or FASTA(.fai): the SBS96 part (block 8) will be skipped")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()

# bcftools: stdout lines; on failure shows bcftools' own message
bcf <- function(args) {
  err <- tempfile()
  out <- tryCatch(suppressWarnings(system2(BCFTOOLS, args, stdout = TRUE, stderr = err)),
                  error = function(e) structure(character(), status = -1L, msg = conditionMessage(e)))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0)
    stop("bcftools failed: bcftools ", paste(args, collapse = " "), "\n",
         paste(c(attr(out, "msg"), readLines(err, warn = FALSE)), collapse = "\n"))
  out
}
save_plot <- function(p, name, w = 9, h = 6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 150, limitsize = FALSE)
  invisible(p)
}
# Callers in INFO/CALLER, matched like bcftools' INFO/CALLER~"mutect2" (name anywhere in the value,
# whatever the separator). Returns a TRUE/FALSE column per caller.
caller_flags <- function(x) setNames(lapply(CALLERS, function(cl) grepl(cl, x, ignore.case = TRUE)), CALLERS)

# ---- 2. Files ----------------------------------------------------------------
files <- data.table(path = list.files(DIR_A, recursive = TRUE, full.names = TRUE, pattern = "\\.(bcf|vcf|vcf\\.gz)$"))
files[, rel := substring(path, nchar(DIR_A) + 2)]
files[, tumour := NA_character_]
for (k in names(TUMOURS)) files[grepl(TUMOURS[[k]], rel, fixed = TRUE), tumour := k]
files[, caller_in_name := NA_character_]
for (cl in CALLERS) files[is.na(caller_in_name) & grepl(cl, rel, ignore.case = TRUE), caller_in_name := cl]   # file or folder name
hdr <- setNames(lapply(files$path, function(f) bcf(c("view", "-h", shQuote(f)))), files$rel)
files[, has_CALLER := sapply(rel, function(r) any(startsWith(hdr[[r]], "##INFO=<ID=CALLER,")))]
# merged = has INFO/CALLER; caller = file named after one caller; else unknown
files[, kind := fcase(has_CALLER, "merged", !is.na(caller_in_name), "caller", default = "unknown")]
files[, .(tumour, kind, caller_in_name, file = rel)]
stopifnot(!anyNA(files$tumour))
files[, .N, by = .(tumour, kind)]

# ---- 3. Read records ---------------------------------------------------------
read_sites <- function(i) {
  f <- files[i]
  fmt <- paste0("%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER", if (f$has_CALLER) "\\t%INFO/CALLER", "\\n")
  cn <- c("CHROM", "POS", "REF", "ALT", "FILTER", if (f$has_CALLER) "CALLER")
  cmd <- paste(shQuote(BCFTOOLS), "query", "-f", shQuote(fmt), shQuote(f$path))
  x <- tryCatch(fread(cmd = cmd, sep = "\t", header = FALSE, col.names = cn, quote = "",
                      colClasses = list(character = setdiff(cn, "POS"), integer = "POS")),
                error = function(e) { message(f$rel, ": no records read (", conditionMessage(e), ")"); data.table() })
  if (!nrow(x)) return(NULL)
  if (!f$has_CALLER) x[, CALLER := f$caller_in_name]
  x[, `:=`(tumour = f$tumour, kind = f$kind, rel = f$rel)]
  x
}
raw <- rbindlist(lapply(seq_len(nrow(files)), read_sites), use.names = TRUE, fill = TRUE)
raw[, PASS := FILTER == "PASS"]
raw[, .(records = .N, PASS = sum(PASS)), by = .(tumour, kind, rel)]        # row counts per file

# 3a. The merged file: one record per variant, INFO/CALLER = the callers that CALLED it (any FILTER).
# NOTE: its FILTER is the bcftools-merge union (see header and R5c), so INFO/CALLER on merged PASS
# records does not mean those callers passed the variant.
merged <- raw[kind == "merged"]
if (nrow(merged)) {
  merged[, (CALLERS) := caller_flags(CALLER)]
  merged[, N_CALLERS := rowSums(.SD), .SDcols = CALLERS]
  merged[, CALLERS_SET := apply(.SD, 1, function(v) paste(CALLERS[v], collapse = "+")), .SDcols = CALLERS]
  print(merged[, .N, by = .(CALLER)][order(-N)])           # raw INFO/CALLER values, to check the matching
  print(merged[N_CALLERS == 0, .N, by = CALLER])           # records with none of CALLERS (should be none)
  # FILTER per caller combination: what PASS means for a merged record
  merged_filter <- merged[, .(records = .N, PASS = sum(PASS)), by = .(tumour, CALLERS_SET)][order(tumour, -records)]
  print(merged_filter)
}

# 3b. Per-caller PASS sets from INFO/CALLER of the merged PASS records.
# Only if there is no merged file, fall back to per-caller files.
if (nrow(merged)) {
  per_caller <- rbindlist(lapply(CALLERS, function(cl)
    merged[PASS == TRUE & get(cl), .(tumour, SET = cl, CHROM, POS, REF, ALT, FILTER)]))
  PER_CALLER_SOURCE <- "merged file, INFO/CALLER (record FILTER = merged FILTER)"
} else {
  per_caller <- raw[kind == "caller" & PASS == TRUE, .(tumour, SET = CALLER, CHROM, POS, REF, ALT, FILTER)]
  PER_CALLER_SOURCE <- "per-caller files (no merged file found)"
}
per_caller <- unique(per_caller[SET %chin% CALLERS])
PER_CALLER_SOURCE
per_caller[, .N, by = .(tumour, SET)][order(tumour, SET)]

# 3c. The union = PASS records of the merged file
union_pass <- if (nrow(merged)) merged[PASS == TRUE, .(tumour, SET = "union", CHROM, POS, REF, ALT, FILTER)] else data.table()
union_pass[, .N, by = tumour]


### FROM RAW -> PASS -> INTERSECT
# Check: rebuild the intersects from the single-caller files and compare with the merged file.
# Needs the single-caller files (kind == "caller" in `files`); skipped if there are none.
rawc <- raw[kind == "caller"]
if (nrow(rawc) && nrow(merged)) {

  # R1. What the single-caller files hold: records per FILTER value (is PASS set there at all?
  #     FreeBayes may have FILTER "." if its crude filter was applied only when merging)
  print(rawc[, .N, by = .(tumour, CALLER, rel, FILTER)][order(tumour, CALLER, -N)])

  # R2. Split multiallelic ALTs (as the provider's bcftools norm -m -any did) and key the variants
  rawc <- rawc[, .(ALT = unlist(strsplit(ALT, ",", fixed = TRUE))), by = .(tumour, CALLER, CHROM, POS, REF, FILTER, PASS)]
  rawc[, ID := paste(CHROM, POS, REF, ALT, sep = ":")]
  merged[, ID := paste(CHROM, POS, REF, ALT, sep = ":")]

  # R3. RAW -> PASS -> INTERSECT: PASS set per caller, then which callers share each variant
  rp <- unique(rawc[PASS == TRUE, .(tumour, CALLER, ID)])
  rp_w <- dcast(rp, tumour + ID ~ CALLER, fun.aggregate = length, value.var = "ID")
  for (cl in CALLERS) if (!cl %in% names(rp_w)) rp_w[, (cl) := 0L]
  rp_w[, N_CALLERS := rowSums(.SD > 0), .SDcols = CALLERS]
  rp_w[, CALLERS_SET := apply(.SD > 0, 1, function(v) paste(CALLERS[v], collapse = "+")), .SDcols = CALLERS]

  # R4. Same question, both routes side by side, per caller combination
  #     from_raw_PASS:    PASS in each single-caller file, then intersected
  #     from_merged_PASS: merged records with FILTER PASS, combination from INFO/CALLER
  #     from_merged_all:  all merged records, any FILTER
  check_raw_vs_merged <- Reduce(function(a, b) merge(a, b, by = c("tumour", "CALLERS_SET"), all = TRUE), list(
    rp_w[, .(from_raw_PASS = .N), by = .(tumour, CALLERS_SET)],
    merged[PASS == TRUE, .(from_merged_PASS = .N), by = .(tumour, CALLERS_SET)],
    merged[, .(from_merged_all = .N), by = .(tumour, CALLERS_SET)]))
  print(check_raw_vs_merged[order(tumour, CALLERS_SET)])
  # Your number: >1 caller and PASS, both routes
  print(merge(rp_w[N_CALLERS > 1, .(raw_PASS_then_intersect = .N), by = tumour],
              merged[N_CALLERS > 1 & PASS == TRUE, .(merged_PASS = .N), by = tumour], by = "tumour", all = TRUE))

  # R5. Where does the merged PASS come from? For each merged record: its state in each caller's own
  #     file (PASS / filtered / absent), next to the merged FILTER
  st <- rawc[, .(state = if (any(PASS)) "PASS" else "filtered"), by = .(tumour, CALLER, ID)]
  st_w <- dcast(st, tumour + ID ~ CALLER, value.var = "state")
  mg <- st_w[merged[, .(tumour, ID, CALLER, merged_PASS = PASS)], on = .(tumour, ID)]
  for (cl in CALLERS) { if (!cl %in% names(mg)) mg[, (cl) := NA_character_]; set(mg, which(is.na(mg[[cl]])), cl, "absent") }
  print(mg[, .N, by = c("merged_PASS", CALLERS)][order(merged_PASS, -N)])
  # Candidate rules for the merged PASS: PASS in any caller, or PASS in every caller that found it
  mg[, any_PASS := Reduce(`|`, lapply(.SD, `==`, "PASS")), .SDcols = CALLERS]
  mg[, all_found_PASS := any_PASS & Reduce(`&`, lapply(.SD, `%chin%`, c("PASS", "absent"))), .SDcols = CALLERS]
  print(mg[, .N, by = .(merged_PASS, any_PASS, all_found_PASS)][order(merged_PASS, any_PASS)])
  # Does INFO/CALLER list every caller that found the variant, or only those where it PASSed?
  for (cl in CALLERS) print(mg[, .N, by = .(in_INFO_CALLER = grepl(cl, CALLER, ignore.case = TRUE), state_in_own_file = get(cl))][
    , caller := cl][order(-N)])

  # R5b. Which rule reproduces the merged FILTER? Agreement of merged PASS with each candidate rule.
  #      Neither "any" nor "all found" fits exactly, so test: FILTER taken from the first caller listed
  #      in INFO/CALLER, or from one fixed caller.
  mg[, first_caller := tolower(trimws(sub("[,|;&+].*$", "", CALLER)))]
  mg[, rule_first_listed := fcase(first_caller == "strelka", strelka == "PASS", first_caller == "mutect2", mutect2 == "PASS",
                                  first_caller == "freebayes", freebayes == "PASS", default = NA)]
  for (cl in CALLERS) mg[get(cl) != "absent", paste0("rule_", cl, "_if_found") := get(cl) == "PASS"]
  rule_cols <- c("any_PASS", "all_found_PASS", "rule_first_listed", paste0("rule_", CALLERS, "_if_found"))
  rule_check <- rbindlist(lapply(rule_cols, function(r) mg[!is.na(get(r)), .(rule = r, records = .N,
    agree = sum(merged_PASS == get(r)), merged_PASS_rule_FALSE = sum(merged_PASS & !get(r)),
    merged_notPASS_rule_TRUE = sum(!merged_PASS & get(r)))]))
  print(rule_check)
  # Order of callers in INFO/CALLER (is it always mutect2 first, then strelka, then freebayes?)
  print(mg[, .N, by = CALLER][order(-N)])

  # R5c. The rule that reproduces the merged FILTER exactly (bcftools merge default): PASS if >=1 caller
  #      is PASS and no caller has a named filter (FreeBayes "." does not count). Expect no TRUE/FALSE mix.
  merge_rule <- rawc[, .(fail = any(!FILTER %chin% c("PASS", ".")), pass = any(FILTER == "PASS")), by = .(tumour, ID)][
    mg, on = .(tumour, ID)][, .N, by = .(merged_PASS, bcftools_merge_rule = pass & !fail)]
  print(merge_rule)

  # R6. The variants where the two routes disagree (>1 caller): raw combination vs merged combination,
  #     and the state in each caller's own file. raw_set NA = not >1-caller PASS in the raw files;
  #     merged_set NA = not a >1-caller PASS record in the merged file; all states "absent" = the
  #     variant is written differently in the raw files (representation, e.g. not split/trimmed).
  diff_routes <- merge(rp_w[N_CALLERS > 1, .(tumour, ID, raw_set = CALLERS_SET)],
                       merged[N_CALLERS > 1 & PASS == TRUE, .(tumour, ID, merged_set = CALLERS_SET, merged_CALLER = CALLER)],
                       by = c("tumour", "ID"), all = TRUE)[is.na(raw_set) | is.na(merged_set) | raw_set != merged_set]
  diff_routes <- st_w[diff_routes, on = .(tumour, ID)]
  for (cl in CALLERS) { if (!cl %in% names(diff_routes)) diff_routes[, (cl) := NA_character_]
                        set(diff_routes, which(is.na(diff_routes[[cl]])), cl, "absent") }
  diff_routes[, merged_FILTER := merged[diff_routes, on = .(tumour, ID), x.FILTER, mult = "first"]]
  print(diff_routes[, .N, by = c("tumour", "raw_set", "merged_set", CALLERS, "merged_FILTER")][order(tumour, -N)])
  print(head(diff_routes[order(tumour)], 20))
  fwrite(diff_routes, file.path(OUT, "raw_vs_merged_disagreements.tsv"), sep = "\t")

  # R7. My intersect (PASS in each caller's own file, then overlapped; >=2 callers) vs the merged file
  #     (FILTER PASS, >=2 callers in INFO/CALLER). Fill = SNV / INDEL.
  #     + Mutect2 & Strelka: PASS in the Mutect2 file AND PASS in the Strelka file (FreeBayes ignored)
  int <- rbind(rp_w[N_CALLERS >= 2, .(tumour, ID, source = "my intersect")],
               rp_w[mutect2 > 0 & strelka > 0, .(tumour, ID, source = "Mutect2 & Strelka")],
               merged[N_CALLERS >= 2 & PASS == TRUE, .(tumour, ID, source = "merged file")])
  int <- unique(int)                     # each variant once per tumour and source (also if >1 merged file)
  int[, type := fifelse(nchar(tstrsplit(ID, ":")[[3]]) == nchar(tstrsplit(ID, ":")[[4]]), "SNV", "INDEL")]
  int[, .N, by = .(tumour, source, type)][order(tumour, source, type)]

  # R7b. Mutect2 & Strelka intersect: the variant list (chrom:pos:ref:alt as in the raw files, MNVs not
  #      split - a Strelka doublet written as two SNVs does not match a Mutect2 MNV) and its counts
  m2_strelka <- int[source == "Mutect2 & Strelka", .(tumour, ID, type)]
  m2_strelka[, c("CHROM", "POS", "REF", "ALT") := tstrsplit(ID, ":", fixed = TRUE, type.convert = FALSE)]
  m2_strelka[, POS := as.integer(POS)]
  m2_strelka_n <- m2_strelka[, .N, by = .(tumour, type)][order(tumour, type)]
  print(m2_strelka_n)
  fwrite(m2_strelka[order(tumour, CHROM, POS)], file.path(OUT, "intersect_mutect2_strelka_PASS.tsv"), sep = "\t")
  fwrite(m2_strelka_n, file.path(OUT, "intersect_mutect2_strelka_PASS_counts.tsv"), sep = "\t")

  # Venn per tumour, SNVs: my intersect vs the merged file (same PASS >=2-caller sets as the bar plot)
  v2 <- dcast(int[type == "SNV" & source != "Mutect2 & Strelka"], tumour + ID ~ source, fun.aggregate = length, value.var = "ID")
  setnames(v2, c("my intersect", "merged file"), c("mine", "mrg"), skip_absent = TRUE)
  venn2 <- v2[, .(only_mine = sum(mine > 0 & mrg == 0), shared = sum(mine > 0 & mrg > 0),
                  only_merged = sum(mine == 0 & mrg > 0), n_mine = sum(mine > 0), n_merged = sum(mrg > 0)), by = tumour]
  print(venn2)
  fwrite(venn2, file.path(OUT, "venn_snv_my_intersect_vs_merged.tsv"), sep = "\t")
  venn2_lab <- melt(venn2, id.vars = "tumour", measure.vars = c("only_mine", "shared", "only_merged"), value.name = "N")
  venn2_lab[, `:=`(x = c(only_mine = -1.0, shared = 0, only_merged = 1.0)[as.character(variable)], y = 0)]
  venn2_names <- melt(venn2, id.vars = "tumour", measure.vars = c("n_mine", "n_merged"), value.name = "n")
  venn2_names[, `:=`(x = fifelse(variable == "n_mine", -0.6, 0.6), y = 1.3,
                     label = paste0(fifelse(variable == "n_mine", "my intersect", "merged file"), "\n(n=", n, ")"))]
  circ2 <- data.table(set = c("my intersect", "merged file"), cx = c(-0.55, 0.55))[
    , .(x = cx + cos(seq(0, 2 * pi, length.out = 200)), y = sin(seq(0, 2 * pi, length.out = 200))), by = set]
  p_venn2 <- ggplot() +
    geom_polygon(data = circ2, aes(x, y, group = set, fill = set), alpha = 0.3, colour = "grey30") +
    geom_text(data = venn2_lab, aes(x, y, label = N), size = 4.5) +
    geom_text(data = venn2_names, aes(x, y, label = label), size = 3.5, fontface = "bold") +
    facet_wrap(~ tumour) + coord_equal(clip = "off") +
    scale_fill_manual(values = c("my intersect" = "#4C72B0", "merged file" = "#DD8452"), guide = "none") +
    labs(title = "SNVs PASS by >=2 callers: my intersect vs merged file") +
    theme_void() + theme(strip.text = element_text(face = "bold", size = 12), plot.margin = margin(10, 20, 10, 20))
  save_plot(p_venn2, "venn_snv_my_intersect_vs_merged", 12, 4.5)
  # Panels nested: tumour, then source (ggh4x::facet_nested; plain facet_grid if ggh4x is missing)
  int[, source := factor(source, levels = c("my intersect", "Mutect2 & Strelka", "merged file"))]
  p_int <- ggplot(int[, .N, by = .(tumour, source, type)], aes(type, N, fill = type)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = N), vjust = -0.3, size = 3.2) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    (if (requireNamespace("ggh4x", quietly = TRUE)) ggh4x::facet_nested(~ tumour + source) else facet_grid(~ tumour + source)) +
    labs(x = NULL, y = "PASS variants, >=2 callers (Mutect2 & Strelka: PASS in both)") + theme_bw() + theme(legend.position = "none")
  save_plot(p_int, "intersect_raw_vs_merged", 14, 5)
} else message("RAW -> PASS -> INTERSECT skipped: no single-caller files (kind == 'caller') or no merged file")



# ---- 4. Variant types; MNVs split into SNVs ------------------------------------
add_type <- function(x) {
  x[, TYPE := fcase(nchar(REF) == 1 & nchar(ALT) == 1, "SNV",
                    nchar(REF) == nchar(ALT),          "MNV",
                    nchar(REF) <  nchar(ALT) & substr(REF, 1, 1) == substr(ALT, 1, 1), "INS",
                    nchar(REF) >  nchar(ALT) & substr(REF, 1, 1) == substr(ALT, 1, 1), "DEL",
                    default = "complex")]
  x[, CLASS := fifelse(TYPE %chin% c("SNV", "MNV"), "SNV", "INDEL")]
  # SNV / INDEL from qcVCF::add_mutational_type() where it works (it knows only SNV / INDEL, so TYPE
  # with MNV / INS / DEL stays ours). Called with the table only; the column it adds is found by name.
  if (HAVE_QCVCF && nrow(x)) {
    y <- copy(x[, .(CHROM, POS, REF, ALT)])
    res <- tryCatch(qcVCF::add_mutational_type(y), error = function(e) { message("qcVCF add_mutational_type not used: ", conditionMessage(e)); NULL })
    if (!is.null(res)) {
      res <- if (is.data.frame(res)) as.data.table(res) else y
      new <- setdiff(names(res), c("CHROM", "POS", "REF", "ALT"))
      if (length(new) == 1 && nrow(res) == nrow(x)) {
        q <- toupper(res[[new]])
        print(table(qcVCF = q, ours = x$CLASS, useNA = "ifany"))   # MNVs: qcVCF may call them INDEL - check
        x[, CLASS_qcVCF := q]
        x[TYPE != "MNV" & CLASS_qcVCF %chin% c("SNV", "INDEL"), CLASS := CLASS_qcVCF]  # MNVs stay SNV (split later)
      } else message("qcVCF add_mutational_type: could not find its output column - CLASS from add_type()")
    }
  }
  x
}
sets <- add_type(rbind(per_caller, union_pass))
sets[, .N, by = .(tumour, SET, TYPE)][order(tumour, SET, TYPE)]   # records as called (MNVs whole)

# Split MNVs into their SNVs (positions where REF and ALT differ)
split_mnv <- function(x) {
  m <- x[TYPE == "MNV"]
  if (!nrow(m)) return(x[, FROM_MNV := FALSE][])
  n <- nchar(m$REF)
  e <- m[rep(seq_len(.N), n)]
  e[, k := sequence(n)]
  e[, `:=`(POS = POS + k - 1L, REF = substring(REF, k, k), ALT = substring(ALT, k, k))]
  e <- e[REF != ALT][, `:=`(TYPE = "SNV", FROM_MNV = TRUE, k = NULL)]
  rbind(x[TYPE != "MNV"][, FROM_MNV := FALSE], e)
}
atom <- split_mnv(copy(sets))
atom[, ID := paste(CHROM, POS, REF, ALT, sep = ":")]
atom <- unique(atom, by = c("tumour", "SET", "ID"))
atom[, .N, by = .(tumour, SET, CLASS)][order(tumour, SET, CLASS)]

SET_LEVELS <- c(CALLERS, "union")
atom[, SET := factor(SET, levels = SET_LEVELS)]

# ---- 5. Counts per caller and the union ---------------------------------------
counts <- sets[, .(records = .N, SNV = sum(TYPE == "SNV"), MNV = sum(TYPE == "MNV"),
                   INS = sum(TYPE == "INS"), DEL = sum(TYPE == "DEL"), complex = sum(TYPE == "complex")),
               by = .(tumour, SET)]
counts[atom[, .(SNV_after_MNV_split = sum(CLASS == "SNV"), INDEL = sum(CLASS == "INDEL")), by = .(tumour, SET)],
       on = .(tumour, SET), `:=`(SNV_after_MNV_split = i.SNV_after_MNV_split, INDEL = i.INDEL)]
counts[, SET := factor(SET, levels = SET_LEVELS)]
setorder(counts, tumour, SET)
counts

p_counts <- ggplot(atom[, .N, by = .(tumour, SET, CLASS)], aes(SET, N, fill = CLASS)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = N), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  facet_wrap(~ tumour, scales = "free_y") +
  labs(x = NULL, y = "PASS variants (MNVs split)", fill = NULL,
       title = paste("Tumour vs matched KI normal: SNVs and indels per caller",
                     "")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_counts, "counts_per_caller", 10, 5)

# qcVCF: counts per caller (MUTTYPE classes SNV, MNV, INS, DEL; records as called)
if (HAVE_QCVCF) {
  library(qcVCF)
  for (t in names(TUMOURS)) {
    q <- sets[tumour == t & SET %chin% CALLERS & TYPE != "complex",
              .(CHROM, POS, REF, ALT, subjectID = tumour, tool = SET, mutType = TYPE)]
    if (nrow(q)) save_plot(plot_mutation_counts(q, grp_col = "tool", cohort = t), paste0("qcVCF_mutation_counts_", t), 8, 6)
  }
}

# ---- 6. Pairwise sharing between callers -------------------------------------
cal <- atom[SET %chin% CALLERS]
cal[, SET := droplevels(SET)]

# 6a. Pairwise table: shared, % of the row caller, Jaccard
pairwise <- cal[, {
  s <- split(ID, SET)
  pr <- CJ(caller_1 = names(s), caller_2 = names(s))
  pr[, `:=`(n_1 = lengths(s[caller_1]), n_2 = lengths(s[caller_2]),
            shared = mapply(function(a, b) length(intersect(s[[a]], s[[b]])), caller_1, caller_2))]
  pr[, `:=`(pct_of_1 = round(100 * shared / n_1, 1), jaccard = round(shared / (n_1 + n_2 - shared), 3))]
  pr
}, by = .(tumour, CLASS)]
pairwise[caller_1 != caller_2]

p_pair <- ggplot(pairwise, aes(caller_2, caller_1, fill = pct_of_1)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%d\n%.0f%%", shared, pct_of_1)), size = 3) +
  facet_grid(CLASS ~ tumour) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100)) +
  labs(x = "found also by", y = "caller (row = 100%)", fill = "% of row",
       title = "Pairwise shared PASS variants between callers") +
  theme_bw()
save_plot(p_pair, "pairwise_shared", 10, 6)

# qcVCF pairwise plots, one for SNVs and one for indels: every tumour x caller set against every other
# (HCC_1_mutect2, HCC_2_mutect2, HCC_1_strelka, ...), so sharing between tumours shows too
qc_pairwise <- data.table()
if (HAVE_QCVCF) {
  qc_pairwise <- rbindlist(lapply(c("SNV", "INDEL"), function(cl) {
    q <- cal[CLASS == cl, .(CHROM, POS, REF, ALT, subjectID = tumour, tool = as.character(SET))]
    q[, sample_tool := paste(subjectID, tool, sep = "_")]
    # qcVCF's by = get() fails on very small classes; skip those instead of stopping
    pw <- tryCatch(plot_pairwise_shared_mutations(q, sample_col = "sample_tool", cohort = cl),
                   error = function(e) { message("qcVCF pairwise ", cl, " skipped: ", conditionMessage(e)); NULL })
    if (is.null(pw)) return(NULL)
    save_plot(pw$plot, paste0("qcVCF_pairwise_shared_", cl), 12, 10)
    as.data.table(pw$data)[, CLASS := cl]
  }), fill = TRUE)
  print(qc_pairwise)
}

# 6b. Exclusive categories: which callers found each variant; number of callers
memb <- dcast(cal, tumour + CLASS + ID ~ SET, fun.aggregate = length, value.var = "ID")
for (cl in CALLERS) if (!cl %in% names(memb)) memb[, (cl) := 0L]
memb[, N_CALLERS := rowSums(.SD > 0), .SDcols = CALLERS]
memb[, CALLERS_SET := apply(.SD > 0, 1, function(v) paste(CALLERS[v], collapse = "+")), .SDcols = CALLERS]
exclusive <- memb[, .N, by = .(tumour, CLASS, CALLERS_SET, N_CALLERS)]
exclusive[, pct := round(100 * N / sum(N), 1), by = .(tumour, CLASS)]
setorder(exclusive, tumour, CLASS, -N)
exclusive
by_n_callers <- memb[, .N, by = .(tumour, CLASS, N_CALLERS)][order(tumour, CLASS, N_CALLERS)]
by_n_callers

lvl <- c(CALLERS, combn(CALLERS, 2, paste, collapse = "+"), paste(CALLERS, collapse = "+"))
p_upset <- ggplot(exclusive, aes(factor(CALLERS_SET, levels = lvl), N, fill = factor(N_CALLERS))) +
  geom_col() + geom_text(aes(label = N), vjust = -0.3, size = 3) +
  facet_grid(CLASS ~ tumour, scales = "free_y") +
  labs(x = "found by (exclusive)", y = "PASS variants", fill = "callers",
       title = "Caller overlap, exclusive categories (UpSet-style)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_upset, "caller_overlap_exclusive", 12, 6)

# ---- 7. Union vs single callers ---------------------------------------------
# Union = PASS records of the merged file. Per caller: how much of its PASS set is in the union, and
# how much of the union it contributes. Since the caller sets come from the same file (INFO/CALLER),
# union_check should show only_merged = only_callers = 0; anything else flags a CALLER value none
# of CALLERS matches (see block 3a).
un <- atom[SET == "union", .(tumour, CLASS, ID)]
cal_union <- unique(cal[, .(tumour, CLASS, ID)])          # union of the callers' PASS sets
union_vs <- rbind(
  cal[, .(n_caller = .N, in_union = sum(ID %chin% un[tumour == .BY$tumour & CLASS == .BY$CLASS, ID])),
      by = .(tumour, CLASS, SET)][, pct_of_caller_in_union := round(100 * in_union / n_caller, 1)],
  fill = TRUE)
union_vs[un[, .(n_union = .N), by = .(tumour, CLASS)], on = .(tumour, CLASS), n_union := i.n_union]
union_vs[, pct_of_union_from_caller := round(100 * in_union / n_union, 1)]
setorder(union_vs, tumour, CLASS, SET)
union_vs

# Merged-union PASS vs union of the callers' PASS sets
union_check <- merge(un[, .(merged_union = .N), by = .(tumour, CLASS)],
                     cal_union[, .(callers_union = .N), by = .(tumour, CLASS)], by = c("tumour", "CLASS"), all = TRUE)
union_check[, both := mapply(function(t, cl) length(intersect(un[tumour == t & CLASS == cl, ID],
                                                              cal_union[tumour == t & CLASS == cl, ID])), tumour, CLASS)]
union_check[, `:=`(only_merged = merged_union - both, only_callers = callers_union - both)]
union_check

# Union by number of callers (from INFO/CALLER of the merged PASS records)
if (nrow(merged)) {
  union_n <- add_type(merged[PASS == TRUE])[, .N, by = .(tumour, CLASS, N_CALLERS)][order(tumour, CLASS, N_CALLERS)]
  union_n[, pct := round(100 * N / sum(N), 1), by = .(tumour, CLASS)]
  print(union_n)
  p_union <- ggplot(union_n, aes(factor(N_CALLERS), N)) + geom_col(fill = "grey30") +
    geom_text(aes(label = paste0(N, "\n(", pct, "%)")), vjust = -0.2, size = 3) +
    facet_grid(CLASS ~ tumour, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = "number of callers (INFO/CALLER)", y = "PASS records in the merged union",
         title = "Merged union: how many callers support each record") + theme_bw()
  save_plot(p_union, "union_by_n_callers", 9, 6)
} else union_n <- data.table()

# ---- 8. SBS96 context ----------------------------------------------------------
# Sets: each caller, the union, >=2 callers, all 3 callers. SNVs after MNV split; doublet halves out.
SBS96_TYPES <- sort(CJ(b5 = c("A", "C", "G", "T"), sub = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                       b3 = c("A", "C", "G", "T"))[, paste0(b5, "[", sub, "]", b3)], method = "radix")
# + strand trinucleotide from the FASTA (samtools faidx); contig names matched to the .fai
trinuc_context <- function(d) {
  fai <- fread(paste0(FASTA, ".fai"), header = FALSE, select = 1)[[1]]
  chrom <- fifelse(d$CHROM %chin% fai, d$CHROM, sub("^chr", "", d$CHROM))
  reg <- paste0(chrom, ":", d$POS - 1L, "-", d$POS + 1L)
  rf <- tempfile(fileext = ".txt"); on.exit(unlink(rf))
  writeLines(unique(reg), rf)
  fa <- system2(SAMTOOLS, c("faidx", "-r", shQuote(rf), shQuote(FASTA)), stdout = TRUE)
  h <- startsWith(fa, ">")
  seqs <- data.table(region = sub("^>", "", fa[h]), NC_3 = toupper(fa[which(h) + 1L]))
  seqs[data.table(region = reg), on = "region"]$NC_3
}
# SBS96 class on the pyrimidine strand, e.g. A[C>T]G (G/A references reverse-complemented)
sbs96 <- function(ref, alt, nc3) {
  rc <- function(x) chartr("ACGT", "TGCA", x)
  pur <- ref %chin% c("A", "G")
  b5 <- fifelse(pur, rc(substr(nc3, 3, 3)), substr(nc3, 1, 1))
  b3 <- fifelse(pur, rc(substr(nc3, 1, 1)), substr(nc3, 3, 3))
  out <- paste0(b5, "[", fifelse(pur, rc(ref), ref), ">", fifelse(pur, rc(alt), alt), "]", b3)
  out[is.na(nc3)] <- NA_character_
  out
}
# Doublet half: another SNV of the same set at the neighbouring position
in_doublet <- function(d) {
  here <- paste(d$tumour, d$SET, d$CHROM, d$POS)
  paste(d$tumour, d$SET, d$CHROM, d$POS - 1L) %chin% here | paste(d$tumour, d$SET, d$CHROM, d$POS + 1L) %chin% here
}

snv96 <- data.table(); m96 <- data.table(); cosine <- data.table(); sbs96_filter_n <- data.table()
if (HAVE_FASTA) {
  snv <- atom[CLASS == "SNV", .(tumour, SET = as.character(SET), CHROM, POS, REF, ALT, ID)]
  snv_cal <- memb[CLASS == "SNV", .(tumour, ID, N_CALLERS)]
  snv <- rbind(snv,
               snv[SET %chin% CALLERS][snv_cal[N_CALLERS >= 2], on = .(tumour, ID), nomatch = NULL][, SET := ">=2 callers"][, N_CALLERS := NULL],
               snv[SET %chin% CALLERS][snv_cal[N_CALLERS == 3], on = .(tumour, ID), nomatch = NULL][, SET := "all 3 callers"][, N_CALLERS := NULL])
  snv <- unique(snv, by = c("tumour", "SET", "ID"))
  snv[, NC_3 := trinuc_context(snv)]
  ref_mismatch <- snv[!is.na(NC_3) & substr(NC_3, 2, 2) != REF, .N]
  if (ref_mismatch) warning(ref_mismatch, " SNVs whose REF is not the FASTA base - is this the genome the callers used?")
  snv[, IN_DBS := in_doublet(snv)]
  snv[!IN_DBS & substr(NC_3, 2, 2) == REF, SBS96 := sbs96(REF, ALT, NC_3)]
  print(snv[, .(SNVs = .N, doublet_halves = sum(IN_DBS), with_SBS96 = sum(!is.na(SBS96))), by = .(tumour, SET)])

  SET96 <- c(CALLERS, "union", ">=2 callers", "all 3 callers")
  snv96 <- snv[!is.na(SBS96)]
  snv96[, SET := factor(SET, levels = SET96)]
  snv96[, SUB := substr(SBS96, 3, 5)]

  # Count matrix: rows = 96 types, columns = tumour x set (SigProfiler format)
  m96 <- dcast(snv96[, .N, by = .(SBS96, col = paste(tumour, SET, sep = "|"))], SBS96 ~ col, value.var = "N", fill = 0)
  m96 <- m96[data.table(SBS96 = SBS96_TYPES), on = "SBS96"]
  for (cl in setdiff(names(m96), "SBS96")) set(m96, which(is.na(m96[[cl]])), cl, 0L)
  setnames(m96, "SBS96", "MutationType")

  # 8b. SBS96 plots with ALL records, split by FILTER: within each set three rows - all SNVs, PASS, not PASS.
  #     (m96 and the cosine similarity above stay PASS-only.) Sets:
  #       callers         every record of the caller's own file; PASS = that file's FILTER
  #       union           every record of the merged file; PASS = merged FILTER
  #       >=2 / all 3     found by >=2 / all 3 callers (any FILTER); PASS = PASS in >=2 / all 3 callers' own files
  #     Multiallelics and MNVs split, doublet halves out.
  split_alt <- function(x) {                       # one row per ALT of multiallelic records
    m <- x[ALT %like% ","]
    if (!nrow(m)) return(x)
    s <- strsplit(m$ALT, ",", fixed = TRUE)
    rbind(x[!ALT %like% ","], m[rep(seq_len(.N), lengths(s))][, ALT := unlist(s)])
  }
  a <- rbind(raw[kind == "caller", .(tumour, SET = CALLER, CHROM, POS, REF, ALT, PASS)],
             raw[kind == "merged", .(tumour, SET = "union", CHROM, POS, REF, ALT, PASS)])
  a <- split_alt(a)[nchar(REF) == nchar(ALT) & REF %like% "^[ACGT]+$" & ALT %like% "^[ACGT]+$"]
  a[, TYPE := fifelse(nchar(REF) == 1L, "SNV", "MNV")]
  a <- split_mnv(a)[, .(PASS = any(PASS)), by = .(tumour, SET, CHROM, POS, REF, ALT)]
  cc <- a[SET %chin% CALLERS, .(n_found = .N, n_pass = sum(PASS)), by = .(tumour, CHROM, POS, REF, ALT)]
  a <- rbind(a,
             cc[n_found >= 2, .(tumour, SET = ">=2 callers", CHROM, POS, REF, ALT, PASS = n_pass >= 2)],
             cc[n_found == length(CALLERS), .(tumour, SET = "all 3 callers", CHROM, POS, REF, ALT, PASS = n_pass == length(CALLERS))])
  rm(cc)
  a[, NC_3 := trinuc_context(a)]
  a <- a[!in_doublet(a) & substr(NC_3, 2, 2) == REF]
  a[, SBS96 := sbs96(REF, ALT, NC_3)]
  a[, `:=`(NC_3 = paste0(substr(SBS96, 1, 1), substr(SBS96, 3, 3), substr(SBS96, 7, 7)),   # pyrimidine strand
           filter = fifelse(PASS, "PASS", "not PASS"), SET = factor(SET, levels = SET96))]
  sbs96_filter_n <- a[, .(SNVs = .N, PASS = sum(PASS), not_PASS = sum(!PASS)), by = .(tumour, SET)][order(tumour, SET)]
  print(sbs96_filter_n)

  # ggplot: rows = set x (all SNVs, PASS, not PASS), columns = substitution; % of each row's SNVs
  FILTER_LEVELS <- c("all SNVs", "PASS", "not PASS")
  prof <- rbind(a[, .N, by = .(tumour, SET, filter = "all SNVs", SBS96)], a[, .N, by = .(tumour, SET, filter, SBS96)])
  prof <- prof[CJ(tumour = unique(prof$tumour), SET = factor(SET96, levels = SET96), filter = FILTER_LEVELS,
                  SBS96 = SBS96_TYPES, unique = TRUE), on = .(tumour, SET, filter, SBS96)]
  prof[is.na(N), N := 0L][, SUB := substr(SBS96, 3, 5)]
  prof[, `:=`(n_row = sum(N), pct = if (sum(N)) 100 * N / sum(N) else 0), by = .(tumour, SET, filter)]
  prof <- prof[n_row > 0]
  prof[, filter := factor(filter, levels = FILTER_LEVELS)]
  setorder(prof, tumour, SET, filter)
  prof[, row := paste0(SET, "\n", filter, " (n=", n_row, ")")]
  prof[, row := factor(row, levels = unique(row))]
  prof[, SBS96 := factor(SBS96, levels = SBS96_TYPES[order(substr(SBS96_TYPES, 3, 5), SBS96_TYPES)])]
  sub_cols <- c("C>A" = "#1EBFF0", "C>G" = "#050708", "C>T" = "#E62725", "T>A" = "#CBCACB", "T>C" = "#A1CF64", "T>G" = "#EDC8C5")
  for (t in unique(prof$tumour)) {
    p96 <- ggplot(prof[tumour == t], aes(SBS96, pct, fill = SUB)) +
      geom_col(width = 0.8) +
      facet_grid(row ~ SUB, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = sub_cols, guide = "none") +
      labs(x = NULL, y = "% of SNVs", title = paste(t, "- SBS96 per set: all SNVs, PASS, not PASS (GRCm38 context)")) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5, family = "mono"),
                         strip.text.y = element_text(angle = 0), panel.spacing.x = unit(0.5, "mm"))
    save_plot(p96, paste0("sbs96_", t), 16, 2 + 1.4 * uniqueN(prof[tumour == t, row]))
  }

  # qcVCF 96-context plot per tumour: the sets stacked in one column, rows all SNVs / PASS / not PASS
  if (HAVE_QCVCF && all(sapply(c("cowplot", "stringr", "stringi", "ggtext"), requireNamespace, quietly = TRUE))) {
    for (t in names(TUMOURS)) {
      q96 <- a[tumour == t][order(SET), .(CHROM, POS, REF, ALT, NC_3, tool = as.character(SET), filter)]
      if (!nrow(q96)) next
      fig <- tryCatch(qcVCF::plot96_matrix(q96, rowsplit = "filter", plotsplitrow = "tool", orderplots = "no", showperc = "yes",
                                           dropempty = "yes", dontshowall = "no", roworder = "all SNVs,PASS,not PASS"),
                      error = function(e) { message("qcVCF plot96_matrix ", t, " skipped: ", conditionMessage(e)); NULL })
      # qcVCF's own layout: 4 cm per facet row + 3 cm per plot, 22 cm wide
      rows <- q96[, uniqueN(filter) + 1L, by = tool]$V1
      for (k in seq_along(fig)) save_plot(fig[[k]], paste0("qcVCF_sbs96_", t, if (k > 1) paste0("_", k)),
                                          22 / 2.54, (sum(rows) * 4 + length(rows) * 3) / 2.54)
    }
  }
  rm(a)

  # Cosine similarity of the SBS96 profiles, per tumour, between sets
  cosine <- rbindlist(lapply(names(TUMOURS), function(t) {
    cols <- grep(paste0("^", t, "\\|"), names(m96), value = TRUE)
    if (length(cols) < 2) return(NULL)
    M <- as.matrix(m96[, ..cols]); colnames(M) <- sub("^[^|]*\\|", "", cols)
    nrm <- sqrt(colSums(M^2)); C <- crossprod(M) / outer(nrm, nrm)
    data.table(tumour = t, set_1 = rownames(C)[row(C)], set_2 = colnames(C)[col(C)], cosine = round(as.vector(C), 3))
  }))
  print(dcast(cosine, tumour + set_1 ~ set_2, value.var = "cosine"))
  p_cos <- ggplot(cosine, aes(set_2, set_1, fill = cosine)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", cosine)), size = 3) + facet_wrap(~ tumour) +
    scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +
    labs(x = NULL, y = NULL, title = "Cosine similarity of SBS96 profiles") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_cos, "sbs96_cosine", 13, 5)
}

# ---- 9. Outputs --------------------------------------------------------------
w <- function(x, name) if (nrow(x)) fwrite(x, file.path(OUT, paste0(name, ".tsv")), sep = "\t")
w(files[, !"path"], "files")
if (nrow(merged)) w(merged_filter, "merged_filter_by_caller_combination")
w(counts, "counts_per_caller")
w(pairwise, "pairwise_shared"); w(qc_pairwise, "qcVCF_pairwise_shared")
w(exclusive, "caller_overlap_exclusive"); w(by_n_callers, "caller_overlap_by_n_callers")
w(union_vs, "union_vs_callers"); w(union_check, "union_merged_vs_callers_union"); w(union_n, "union_by_n_callers")
w(m96, "sbs96_matrix"); w(cosine, "sbs96_cosine"); w(sbs96_filter_n, "sbs96_snvs_by_filter")

saveRDS(list(files = files, raw = raw, per_caller = per_caller, per_caller_source = PER_CALLER_SOURCE,
             union_pass = union_pass, sets = sets, atom = atom, counts = counts, pairwise = pairwise,
             memb = memb, exclusive = exclusive, union_vs = union_vs, union_check = union_check, union_n = union_n,
             snv96 = snv96, m96 = m96, cosine = cosine,
             params = list(FASTA = FASTA)),
        file.path(OUT, paste0("callers", ".rds")))

fwrite(data.table(parameter = c("DIR_A", "OUT", "FASTA", "BCFTOOLS", "SAMTOOLS",
                                "PER_CALLER_SOURCE", paste0("TUMOUR_", names(TUMOURS))),
                  value = c(DIR_A, OUT, FASTA, BCFTOOLS, SAMTOOLS, PER_CALLER_SOURCE, TUMOURS)),
       file.path(OUT, "parameters.tsv"), sep = "\t")
fwrite(data.table(file = files$rel, md5 = unname(tools::md5sum(files$path))), file.path(OUT, "input_md5.tsv"), sep = "\t")
git_hash <- tryCatch(system2("git", c("-C", shQuote(REPO), "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) "unknown")
writeLines(c(R.version.string, bcf("--version")[1],
             if (file.exists(SAMTOOLS)) system2(SAMTOOLS, "--version", stdout = TRUE)[1] else "samtools not found",
             paste("data.table", packageVersion("data.table")), paste("ggplot2", packageVersion("ggplot2")),
             if (HAVE_QCVCF) paste("qcVCF", packageVersion("qcVCF")) else "qcVCF not installed",
             paste("FASTA", FASTA), paste("git", git_hash)), file.path(OUT, "versions.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

md_table <- function(x) {
  x <- as.data.table(x)
  c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
writeLines(c(
  "# mouse_wgs_snv_indel_callers summary", "",
  paste0("Run ", format(started, "%Y-%m-%d %H:%M"), "; ", "genome-wide",
         ". Per-caller sets from: ", PER_CALLER_SOURCE, ". PASS records only; MNVs split into SNVs."), "",
  "## Counts per caller and the union", "", md_table(counts), "",
  "## Pairwise shared (row caller = 100%)", "", md_table(pairwise[caller_1 != caller_2]), "",
  "## Variants by number of callers", "", md_table(by_n_callers), "",
  "## Union vs single callers", "", md_table(union_vs), "", md_table(union_check), "",
  if (nrow(cosine)) c("## SBS96 cosine similarity", "", md_table(dcast(cosine, tumour + set_1 ~ set_2, value.var = "cosine")), "") else
    "SBS96 skipped (no FASTA/samtools).",
  "## Notes", "",
  "- Agreement between callers shows robustness to the caller, not somatic origin.",
  "- FreeBayes PASS is a crude prefilter, not a somatic model; expect more germline/artefact calls.",
  "- Indels not left-aligned: FreeBayes indels may miss matches with Mutect2/Strelka for representation alone.",
  "- Doublet halves are left out of SBS96 (counted as one DBS event).",
  "- Counts, not mutations per Mb (no BAMs, no callable territory)."),
  file.path(OUT, "summary.md"))

log_file <- file.path(OUT, "RUN_LOG.md")
old_log <- if (file.exists(log_file)) readLines(log_file) else character()
writeLines(c(
  paste0("## ", format(Sys.time(), "%Y-%m-%d %H:%M"), " - mouse_wgs_snv_indel_callers"),
  paste0("- Run by:     ", Sys.info()[["user"]], "@", Sys.info()[["nodename"]], " (", R.version.string, ")"),
  paste0("- Script:     ", SCRIPT, "  git: ", git_hash),
  paste0("- Parameters: per-caller sets from ", PER_CALLER_SOURCE),
  paste0("- Inputs:     ", DIR_A, " (", nrow(files), " files; provider-split, not left-aligned); FASTA ", FASTA),
  paste0("- Outputs:    ", OUT),
  paste0("- Counts/QC:  records read ", nrow(raw), ", PASS per-caller rows ", nrow(per_caller), ", union PASS ", nrow(union_pass)),
  "", old_log), log_file)

message("Done in ", round(difftime(Sys.time(), started, units = "mins"), 1), " min. Outputs: ", OUT)
