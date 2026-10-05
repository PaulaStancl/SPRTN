# ---------------------------------------------------------------------------
# 02_sv_cnv.R - structural variants and copy number, RJALS tumour vs normal.
# SVs: Manta (sarek), ESVEE + LINX (oncoanalyser).  CNV: ASCAT (sarek), PURPLE.
# Run:  Rscript 02_sv_cnv.R   (or line by line in R - any of the folders above works)
# ---------------------------------------------------------------------------
# Finds 00_setup.R from this folder, from scripts/, or from the project root.
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
od <- file.path(OUT, "sv_cnv"); dir.create(od, recursive = TRUE, showWarnings = FALSE)

# ---- 1. input files ---------------------------------------------------------
ascat_dir <- file.path(SAREK, "variant_calling/ascat", PAIR)
files <- c(
  manta       = find_one(file.path(SAREK, "variant_calling/manta", PAIR), "\\.manta\\.somatic_sv\\.vcf\\.gz$", required = FALSE),
  esvee       = find_one(file.path(ONCO, "esvee"), "\\.esvee\\.somatic\\.vcf\\.gz$", required = FALSE),
  linx_svs    = find_one(file.path(ONCO, "linx"), "\\.linx\\.svs\\.tsv$", required = FALSE),
  ascat_seg   = find_one(ascat_dir, "\\.segments\\.txt$"),
  ascat_pp    = find_one(ascat_dir, "\\.purityploidy\\.txt$"),
  purple_cnv  = find_one(file.path(ONCO, "purple"), "\\.purple\\.cnv\\.somatic\\.tsv$"),
  purple_pur  = find_one(file.path(ONCO, "purple"), "\\.purple\\.purity\\.tsv$")
)
files <- files[!is.na(files)]
writeLines(paste(names(files), files, sep = "\t"), file.path(od, "inputs_used.tsv"))

# ---- 2. structural variants --------------------------------------------------
sv <- bind_rows(lapply(intersect(c("manta", "esvee"), names(files)), function(k)
  read_sv_vcf(files[[k]]) |> mutate(caller = k))) |>
  mutate(svtype = coalesce(svtype, "unknown")) |> filter(chr %in% STD_CHR)
sv_counts <- sv |> group_by(caller, svtype) |>
  summarise(total = n(), pass = sum(filter == "PASS"), .groups = "drop")
print(sv_counts); write_csv(sv_counts, file.path(od, "sv_counts.csv"))
save_plot(ggplot(filter(sv, filter == "PASS"), aes(svtype, fill = caller)) + geom_bar(position = "dodge") +
            labs(x = NULL, y = "PASS records", title = "SVs by type (BND = two records per event)"), "sv_types", od)

# TODO: match Manta and ESVEE breakpoints (StructuralVariantAnnotation) for a real overlap.
if ("linx_svs" %in% names(files)) {
  linx <- read_tsv(files[["linx_svs"]], show_col_types = FALSE)
  if ("type" %in% names(linx)) print(count(linx, type))     # LINX's classification of each SV
}

# ---- 3. copy number: ASCAT vs PURPLE -----------------------------------------
asc <- read_tsv(files[["ascat_seg"]], show_col_types = FALSE)
need_cols(asc, c("chr", "startpos", "endpos", "nMajor", "nMinor"), "ASCAT segments")
pur <- read_tsv(files[["purple_cnv"]], show_col_types = FALSE)
need_cols(pur, c("chromosome", "start", "end", "copyNumber", "minorAlleleCopyNumber"), "PURPLE cnv")

seg <- bind_rows(
  transmute(asc, source = "ASCAT",  chr = norm_chr(chr), start = startpos, end = endpos,
            cn = nMajor + nMinor, minor = nMinor),
  transmute(pur, source = "PURPLE", chr = norm_chr(chromosome), start = start, end = end,
            cn = copyNumber, minor = minorAlleleCopyNumber)
) |> filter(chr %in% STD_CHR) |>
  mutate(x0 = start + unname(CHR_OFFSET[chr]), x1 = end + unname(CHR_OFFSET[chr]))
write_csv(seg, file.path(od, "cn_segments.csv"))

save_plot(
  ggplot(seg, aes(x = x0, xend = x1, y = cn, yend = cn)) +
    geom_vline(xintercept = CHR_OFFSET, colour = "grey85", linewidth = 0.2) +
    geom_segment(linewidth = 1) + facet_grid(source ~ .) +
    scale_x_continuous(breaks = CHR_OFFSET + CHR_LEN / 2, labels = sub("chr", "", STD_CHR), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 8)) + labs(x = "chromosome", y = "total copy number", title = "Copy number"),
  "copy_number", od, w = 11, h = 5)

# Fraction of the autosomal genome per copy-number state, and LOH.
cn_frac <- seg |> filter(!chr %in% c("chrX", "chrY")) |>
  mutate(w = end - start + 1, loh = round(minor) == 0) |>
  group_by(source) |> summarise(LOH = sum(w[loh]) / sum(w), .groups = "drop")
state_frac <- seg |> filter(!chr %in% c("chrX", "chrY")) |>
  mutate(w = end - start + 1, state = pmin(round(cn), 6)) |>
  group_by(source, state) |> summarise(w = sum(w), .groups = "drop") |>
  group_by(source) |> mutate(fraction = w / sum(w)) |> ungroup() |> select(-w)
print(state_frac); print(cn_frac)
write_csv(state_frac, file.path(od, "cn_state_fraction.csv"))
write_csv(cn_frac, file.path(od, "loh_fraction.csv"))

# ---- 4. purity / ploidy, side by side ----------------------------------------
a  <- read_table(files[["ascat_pp"]], show_col_types = FALSE)
pp <- read_tsv(files[["purple_pur"]], show_col_types = FALSE)
purity <- tibble(source = c("ASCAT", "PURPLE"),
                 purity = c(a$AberrantCellFraction[1], pp$purity[1]),
                 ploidy = c(a$Ploidy[1], pp$ploidy[1]))
print(purity); write_csv(purity, file.path(od, "purity_ploidy.csv"))

# ---- 5. qcVCR -----------------------------------------------------------------
# TODO(Paula): call your qcVCR functions here and tell me what they take.
if (HAVE_QCVCR) {
  library(qcVCR)
  # qc <- qcVCR::<function>(files[["manta"]])
}
message("done: ", od)
