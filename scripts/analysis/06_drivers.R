# ---------------------------------------------------------------------------
# 06_drivers.R - driver candidates from oncoanalyser (PURPLE + LINX) next to tumourevo's
# driver annotation, per gene and per mutation.
#   Rscript 06_drivers.R          (r-mutation env; after 03 for the caller columns)
#
# What each source calls a driver:
#   PURPLE  driver catalogue: Hartwig driver gene panel; per gene a driver likelihood (dN/dS,
#           hotspot, in-frame, biallelic) for mutations, plus AMP / DEL from copy number.
#           Germline catalogue too (germline drivers in the panel).
#   LINX    driver catalogue (disruptions, e.g. HOM_DISRUPTION) and reported fusions.
#   tumourevo  sarek Mutect2 PASS calls with VEP 115; is_driver = gene in the IntOGen
#           compendium for HCC (CANCER_TYPE) AND VEP impact MODERATE or HIGH. That says the
#           gene is a known HCC driver gene, not that the variant itself is a known driver.
# PURPLE's somatic VCF adds, per mutation: PAVE's canonical effect and HGVS (INFO/IMPACT),
# SAGE's TIER (HOTSPOT = a known hotspot position) and REPORTED (PURPLE reports it as a driver).
#
# Output: <OUT>/drivers/
#   drivers_by_gene.csv       one row per gene found by any source: PURPLE (type, likelihood,
#                             method, biallelic, CN), LINX, tumourevo (variants), sources
#   driver_mutations.csv      one row per driver mutation (PURPLE-reported, hotspot, or
#                             tumourevo is_driver): PAVE and VEP annotation, VAF, which
#                             pipeline has it, which callers (03's mutation_callers.csv)
#   hcc_watchlist.csv         every PASS mutation in well-known HCC driver genes and the TERT
#                             promoter, in either pipeline, driver-labelled or not
#   linx_fusions_reported.csv LINX fusions with reported = true
#   drivers_overview.pdf      gene x source
#   drivers_somatic_table.csv / .md    slide table: somatic driver mutations, CN and SV drivers -
#                             Gene, Variant (c. / p.), Consequence, VAF, Callers, Driver likelihood,
#                             Evidence (pipeline), Annotation + Tier (empty: knowledge-base step)
#   drivers_germline_table.csv / .md   slide table: PURPLE germline drivers + SPRTN (not on
#                             Hartwig's germline panel; read counts from the SPRTN slice BAMs, phasing
#                             from 12_whatshap_sprtn.sh) - Zygosity, VAF normal -> tumour, tumour second
#                             hit (LOH / allelic imbalance), ClinVar, gene-disease link
# ---------------------------------------------------------------------------
source(Filter(file.exists, c("00_setup.R", "analysis/00_setup.R", "scripts/analysis/00_setup.R",
  "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/analysis/00_setup.R"))[1])
od <- file.path(OUT, "drivers"); dir.create(od, recursive = TRUE, showWarnings = FALSE)

# Well-known HCC driver genes (TCGA-LIHC 2017, Schulze 2015 / IntOGen) and the TERT promoter
# (hotspots C228T chr5:1,295,113 and C250T chr5:1,295,135, GRCh38; window around them)
HCC_GENES <- c("TERT", "CTNNB1", "TP53", "AXIN1", "ARID1A", "ARID2", "ALB", "APOB", "KEAP1", "NFE2L2",
               "RB1", "CDKN2A", "TSC1", "TSC2", "PTEN", "BAP1", "KMT2C", "KMT2D", "RPS6KA3", "NRAS",
               "MYC", "CCND1", "FGF19", "VEGFA", "MET")
TERT_PROMOTER <- list(chrom = "chr5", from = 1295000L, to = 1295400L)

# ---- 1. PURPLE + LINX ------------------------------------------------------------------------
pd <- file.path(ONCO, "purple"); ld <- file.path(ONCO, "linx")
read_cat <- function(f, src) if (is.na(f)) NULL else fread(f)[, source := src]
cat_all <- rbindlist(list(
  read_cat(find_one(pd, "\\.purple\\.driver\\.catalog\\.somatic\\.tsv$",  required = FALSE), "PURPLE somatic"),
  read_cat(find_one(pd, "\\.purple\\.driver\\.catalog\\.germline\\.tsv$", required = FALSE), "PURPLE germline"),
  read_cat(find_one(ld, "\\.linx\\.driver\\.catalog\\.tsv$",              required = FALSE), "LINX")), fill = TRUE)
if (nrow(cat_all)) {
  cat_all[, label := sprintf("%s %s (%s, likelihood %.2f%s)", tolower(source), driver, likelihoodMethod,
                             driverLikelihood, fifelse(as.character(biallelic) %chin% c("true", "TRUE"), ", biallelic", ""))]
  print(cat_all[, .(source, gene, driver, category, likelihoodMethod, driverLikelihood, biallelic, minCopyNumber, maxCopyNumber)])
}
fus_f <- find_one(ld, "\\.linx\\.fusion\\.tsv$", required = FALSE)
fus <- if (!is.na(fus_f)) fread(fus_f) else data.table()
if (nrow(fus) && "reported" %in% names(fus)) {
  fus <- fus[as.character(reported) %chin% c("true", "TRUE")]
  fwrite(fus, file.path(od, "linx_fusions_reported.csv"))
  message(nrow(fus), " reported LINX fusions")
}

# PURPLE's PASS mutations with PAVE's canonical impact (INFO/IMPACT: gene, transcript, effect,
# coding effect, splice region, HGVS c., HGVS p., other effects, worst coding effect, genes)
pv_f <- find_one(pd, "\\.purple\\.somatic\\.vcf\\.gz$")
pv <- read_vcf_table(pv_f)[FILTER == "PASS" & CHROM %chin% STD_CHR]
imp <- if ("info_IMPACT" %in% names(pv)) tstrsplit(as.character(pv$info_IMPACT), ",", fixed = TRUE, fill = NA_character_) else list()
pv[, `:=`(GENE = if (length(imp)) imp[[1]] else NA_character_,
          PAVE_EFFECT = if (length(imp) >= 3) imp[[3]] else NA_character_,
          PAVE_CODING = if (length(imp) >= 4) imp[[4]] else NA_character_,
          PAVE_HGVSC = if (length(imp) >= 6) imp[[6]] else NA_character_,
          PAVE_HGVSP = if (length(imp) >= 7) imp[[7]] else NA_character_,
          TIER = if ("info_TIER" %in% names(pv)) as.character(info_TIER) else NA_character_,
          REPORTED = if ("info_REPORTED" %in% names(pv)) as.character(info_REPORTED) %chin% c("TRUE", "true", "1") else FALSE)]
pv[, KEY := paste(CHROM, POS, REF, ALT1, sep = ":")]

# ---- 2. tumourevo ------------------------------------------------------------------------------
# driver_annotation/annotate_driver/<dataset>/<patient>/<tumour>/*_driver.rds: a list per sample,
# $mutations holds the Mutect2 PASS calls with VEP fields and is_driver / driver_label
tv_f <- list.files(file.path(TEVO, "driver_annotation"), "_driver\\.rds$", recursive = TRUE, full.names = TRUE)
tv <- data.table()
if (length(tv_f)) {
  obj <- readRDS(tv_f[1])
  hit <- Filter(function(x) is.list(x) && !is.null(x$mutations) && "is_driver" %in% names(x$mutations), obj)
  if (length(hit)) {
    tv <- as.data.table(hit[[1]]$mutations)
    tv <- tv[, !vapply(tv, is.list, logical(1)), with = FALSE]          # drop list columns (additional_info)
    tv[, `:=`(CHROM = as.character(chr), POS = as.integer(from), REF = ref, ALT1 = alt)]
    tv[, KEY := paste(CHROM, POS, REF, ALT1, sep = ":")]
    message("tumourevo: ", nrow(tv), " Mutect2 PASS mutations, ", tv[is_driver == TRUE, .N], " is_driver (",
            paste(unique(na.omit(tv$TUMOUR_TYPE)), collapse = ","), " driver genes)")
  }
} else message("no tumourevo *_driver.rds under ", file.path(TEVO, "driver_annotation"), " - tumourevo columns empty")
tcols <- function(...) intersect(c(...), names(tv))

# ---- 3. callers per mutation (03) ------------------------------------------------------------
mc_f <- file.path(OUT, "comparison", "overlap", "mutation_callers.csv")
mc <- if (file.exists(mc_f)) fread(mc_f, select = c("KEY", "CALLERS", "N_CALLERS")) else data.table(KEY = character(), CALLERS = character(), N_CALLERS = integer())
if (!file.exists(mc_f)) message("no ", mc_f, " - run 03 first for the caller columns")

# ---- 4. driver mutations -----------------------------------------------------------------------
drv_genes <- unique(cat_all[driver %chin% c("MUTATION", "GERMLINE_MUTATION") | grepl("MUTATION", driver), gene])
p_drv <- pv[REPORTED | TIER == "HOTSPOT" | GENE %chin% drv_genes,
            .(KEY, CHROM, POS, REF, ALT = ALT1, GENE, PAVE_EFFECT, PAVE_CODING, PAVE_HGVSC, PAVE_HGVSP,
              TIER, REPORTED, VAF_SAGE = round(VAF, 3), in_PURPLE = TRUE)]
t_drv <- if (nrow(tv)) tv[is_driver == TRUE, c("KEY", "CHROM", "POS", "REF", "ALT1", tcols("SYMBOL", "Consequence", "IMPACT", "HGVSc", "HGVSp", "driver_label", "VAF")), with = FALSE] else data.table(KEY = character())
if (nrow(t_drv)) { setnames(t_drv, c("ALT1", "SYMBOL", "Consequence", "IMPACT", "HGVSc", "HGVSp", "VAF"),
                            c("ALT", "VEP_SYMBOL", "VEP_CONSEQUENCE", "VEP_IMPACT", "VEP_HGVSC", "VEP_HGVSP", "VAF_MUTECT2"), skip_absent = TRUE)
                   t_drv[, in_tumourevo_driver := TRUE] }
dm <- merge(p_drv, t_drv, by = intersect(c("KEY", "CHROM", "POS", "REF", "ALT"), names(t_drv)), all = TRUE)
# a PURPLE driver mutation that tumourevo has but did not label (gene not in IntOGen HCC, or LOW impact)
if (nrow(tv)) dm[tv, on = "KEY", `:=`(in_Mutect2_PASS = TRUE, tumourevo_is_driver = i.is_driver)]
dm[, `:=`(in_PURPLE = !is.na(in_PURPLE) & in_PURPLE, in_tumourevo_driver = if ("in_tumourevo_driver" %in% names(dm)) !is.na(in_tumourevo_driver) else FALSE)]
dm[, GENE := fcoalesce(GENE, if ("VEP_SYMBOL" %in% names(dm)) VEP_SYMBOL else NA_character_)]
dm[, found_by := fcase(in_PURPLE & in_tumourevo_driver, "both",
                       in_PURPLE, "PURPLE only", default = "tumourevo only")]
dm[mc, on = "KEY", `:=`(CALLERS = i.CALLERS, N_CALLERS = i.N_CALLERS)]
dm[, REPORTED := !is.na(REPORTED) & REPORTED]
setorder(dm, -REPORTED, GENE, CHROM, POS, na.last = TRUE)
fwrite(dm, file.path(od, "driver_mutations.csv"))
print(dm[, .SD, .SDcols = intersect(c("GENE", "KEY", "PAVE_HGVSP", "VEP_HGVSP", "TIER", "REPORTED", "found_by", "CALLERS"), names(dm))])

# ---- 5. per gene -------------------------------------------------------------------------------
g_p <- if (nrow(cat_all)) cat_all[, .(PURPLE_LINX = paste(unique(label), collapse = "; "),
                                     likelihood = max(driverLikelihood), driver_type = paste(unique(driver), collapse = "+"),
                                     CN = paste0(round(min(minCopyNumber), 2), "-", round(max(maxCopyNumber), 2))), by = gene] else data.table(gene = character())
g_t <- if (nrow(tv)) tv[is_driver == TRUE, .(tumourevo = paste(unique(driver_label), collapse = "; "), tumourevo_n = .N),
                        by = .(gene = SYMBOL)] else data.table(gene = character())
g_f <- if (nrow(fus)) fus[, .(gene = unlist(strsplit(name, "_", fixed = TRUE))), by = name][
                            , .(LINX_fusion = paste(unique(name), collapse = "; ")), by = gene] else data.table(gene = character())
gl <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), list(g_p, g_t, g_f))
for (cl in c("PURPLE_LINX", "tumourevo", "LINX_fusion")) if (!cl %in% names(gl)) gl[, (cl) := NA_character_]
if (!"likelihood" %in% names(gl)) gl[, likelihood := NA_real_]
if (nrow(gl)) {
  gl[, sources := paste(c("PURPLE/LINX", "tumourevo", "LINX fusion")[c(!is.na(PURPLE_LINX), !is.na(tumourevo), !is.na(LINX_fusion))],
                        collapse = " + "), by = gene]
  gl[, HCC_gene := gene %chin% HCC_GENES]
  setorder(gl, -likelihood, gene, na.last = TRUE)
}
fwrite(gl, file.path(od, "drivers_by_gene.csv"))
print(gl)

# ---- 6. HCC watchlist: every PASS mutation in well-known HCC genes / TERT promoter --------------
in_tert <- function(chrom, pos) chrom == TERT_PROMOTER$chrom & pos >= TERT_PROMOTER$from & pos <= TERT_PROMOTER$to
wl <- rbindlist(list(
  pv[GENE %chin% HCC_GENES | in_tert(CHROM, POS),
     .(source = "PURPLE (SAGE)", GENE = fifelse(in_tert(CHROM, POS), "TERT promoter", GENE), KEY, effect = PAVE_EFFECT,
       hgvsp = PAVE_HGVSP, TIER, REPORTED, VAF = round(VAF, 3))],
  if (nrow(tv)) tv[SYMBOL %chin% HCC_GENES | in_tert(CHROM, POS),
     .(source = "tumourevo (Mutect2)", GENE = fifelse(in_tert(CHROM, POS), "TERT promoter", SYMBOL), KEY,
       effect = if ("Consequence" %in% names(tv)) Consequence else NA_character_,
       hgvsp = if ("HGVSp" %in% names(tv)) HGVSp else NA_character_, is_driver,
       VAF = if ("VAF" %in% names(tv)) round(VAF, 3) else NA_real_)]), fill = TRUE)
if (nrow(wl)) wl[mc, on = "KEY", CALLERS := i.CALLERS]
fwrite(wl[order(GENE, KEY)], file.path(od, "hcc_watchlist.csv"))
message("HCC watchlist: ", nrow(wl), " PASS records in ", uniqueN(wl$GENE), " genes")
cn_hcc <- if (nrow(cat_all)) cat_all[gene %chin% HCC_GENES & driver %chin% c("AMP", "PARTIAL_AMP", "DEL")] else data.table()
if (nrow(cn_hcc)) { message("HCC genes with copy-number drivers:"); print(cn_hcc[, .(gene, driver, minCopyNumber, maxCopyNumber)]) }

# ---- 7. overview plot: gene x source ----------------------------------------------------------------
pl <- rbindlist(list(
  if (nrow(cat_all)) cat_all[, .(gene, source, what = driver, likelihood = driverLikelihood)],
  if (nrow(tv)) tv[is_driver == TRUE, .(source = "tumourevo (IntOGen HCC)", what = "MUTATION", likelihood = NA_real_), by = .(gene = SYMBOL)],
  if (nrow(fus)) fus[, .(gene = unlist(strsplit(name, "_", fixed = TRUE))), by = name][
                      , .(gene, source = "LINX fusion", what = "FUSION", likelihood = NA_real_)]), fill = TRUE)
if (nrow(pl)) {
  pl <- unique(pl, by = c("gene", "source", "what"))
  pl[, gene := factor(gene, levels = rev(sort(unique(gene))))]
  save_plot(ggplot(pl, aes(source, gene, fill = what)) + geom_tile(colour = "white") +
              geom_text(aes(label = fifelse(is.na(likelihood), "", sprintf("%.2f", likelihood))), size = 2.6) +
              scale_x_discrete(position = "top") +
              labs(x = NULL, y = NULL, fill = "driver type",
                   title = "Driver candidates: oncoanalyser (PURPLE, LINX) vs tumourevo",
                   subtitle = "number = PURPLE/LINX driver likelihood\ntumourevo = IntOGen HCC driver gene with a MODERATE/HIGH VEP impact") +
              theme(axis.text.x = element_text(angle = 20, hjust = 0)),
            "drivers_overview", od, w = 8, h = 2.5 + 0.25 * uniqueN(pl$gene))
}

# ---- 8. slide tables: somatic and germline ------------------------------------------------------
strip_tx <- function(x) fifelse(is.na(x) | x == "", NA_character_, sub("^[^:]*:", "", x))   # ENST...:c.1A>G -> c.1A>G
caller_name <- c(muse = "MuSE", mutect2 = "Mutect2", sage = "SAGE", strelka = "Strelka2")
pretty_callers <- function(cl, is_snv) vapply(seq_along(cl), function(i) {
  if (is.na(cl[i])) return(NA_character_)
  k <- strsplit(cl[i], "+", fixed = TRUE)[[1]]
  sprintf("%s (%d of %d)", paste(fcoalesce(caller_name[k], k), collapse = ", "), length(k), if (is_snv[i]) 4L else 3L)
}, character(1))
fmt_lik <- function(x) fifelse(is.na(x), "—", sprintf("%.2f", x))
write_md <- function(d, f, note) writeLines(c(paste("|", paste(names(d), collapse = " | "), "|"),
  paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
  apply(d, 1, function(r) paste("|", paste(fifelse(is.na(r), "", r), collapse = " | "), "|")), "", note), f)
som_cat <- if (nrow(cat_all)) cat_all[source == "PURPLE somatic"] else data.table()
lik_mut <- if (nrow(som_cat)) som_cat[grepl("MUTATION", driver), .(lik = max(driverLikelihood)), by = gene] else data.table(gene = character(), lik = numeric())

# somatic: driver mutations (06 section 4) + TERT promoter + CN drivers + LINX disruptions / fusions
tert <- pv[in_tert(CHROM, POS), .(KEY, CHROM, POS, REF, ALT = ALT1, GENE = "TERT promoter", PAVE_EFFECT = "promoter (non-coding)",
                                   PAVE_HGVSC, PAVE_HGVSP, TIER, REPORTED, VAF_SAGE = round(VAF, 3), in_PURPLE = TRUE)]
mt <- rbindlist(list(dm, tert[!KEY %chin% dm$KEY]), fill = TRUE)
mt[mc, on = "KEY", `:=`(CALLERS = i.CALLERS)]
som <- if (nrow(mt)) mt[, .(
  Gene = GENE,
  `Variant (c. / p.)` = paste0(fcoalesce(strip_tx(PAVE_HGVSC), strip_tx(if ("VEP_HGVSC" %in% names(mt)) VEP_HGVSC else NA_character_), paste0(REF, ">", ALT)),
                               fifelse(!is.na(fcoalesce(strip_tx(PAVE_HGVSP), strip_tx(if ("VEP_HGVSP" %in% names(mt)) VEP_HGVSP else NA_character_))) &
                                         fcoalesce(strip_tx(PAVE_HGVSP), "") != "",
                                       paste0(" / ", fcoalesce(strip_tx(PAVE_HGVSP), strip_tx(if ("VEP_HGVSP" %in% names(mt)) VEP_HGVSP else NA_character_))), "")),
  Consequence = gsub("_", " ", fcoalesce(PAVE_EFFECT, if ("VEP_CONSEQUENCE" %in% names(mt)) VEP_CONSEQUENCE else NA_character_)),
  VAF = sprintf("%.2f", fcoalesce(VAF_SAGE, if ("VAF_MUTECT2" %in% names(mt)) as.numeric(VAF_MUTECT2) else NA_real_)),
  Callers = pretty_callers(CALLERS, nchar(REF) == 1 & nchar(ALT) == 1),
  `Driver likelihood` = fmt_lik(lik_mut$lik[match(GENE, lik_mut$gene)]),
  `Evidence (pipeline)` = trimws(paste0(fifelse(REPORTED %in% TRUE, "PURPLE reported; ", ""),
                                        fifelse(TIER %chin% "HOTSPOT", "known hotspot (SAGE); ", ""),
                                        fifelse(if ("in_tumourevo_driver" %in% names(mt)) in_tumourevo_driver %in% TRUE else FALSE, "IntOGen HCC gene (tumourevo)", "")), whitespace = "[ ;]"),
  rank = fifelse(REPORTED %in% TRUE, 1L, 2L))] else data.table()
cn_rows <- if (nrow(som_cat)) som_cat[driver %chin% c("AMP", "PARTIAL_AMP", "DEL", "HOM_DEL_DISRUPTION", "HOM_DUP_DISRUPTION"), .(
  Gene = gene, `Variant (c. / p.)` = fifelse(grepl("AMP", driver), sprintf("%s, CN %.1f", driver, maxCopyNumber), sprintf("%s, CN %.1f", driver, minCopyNumber)),
  Consequence = fifelse(grepl("AMP", driver), "amplification", "deletion / loss"), VAF = "—", Callers = "PURPLE (copy number)",
  `Driver likelihood` = fmt_lik(driverLikelihood), `Evidence (pipeline)` = sprintf("PURPLE driver (%s)", fifelse(category == "TSG", "tumour suppressor", "oncogene")), rank = 1L)] else data.table()
sv_rows <- rbindlist(list(
  if (nrow(cat_all)) cat_all[source == "LINX", .(Gene = gene, `Variant (c. / p.)` = driver, Consequence = "structural disruption", VAF = "—",
                                                Callers = "LINX (SV)", `Driver likelihood` = fmt_lik(driverLikelihood), `Evidence (pipeline)` = "LINX driver catalogue", rank = 1L)],
  if (nrow(fus)) fus[, .(Gene = name, `Variant (c. / p.)` = paste("fusion", name), Consequence = "gene fusion", VAF = "—", Callers = "LINX (SV)",
                         `Driver likelihood` = if ("likelihood" %in% names(fus)) as.character(likelihood) else "—",
                         `Evidence (pipeline)` = paste("LINX reported", if ("reportedType" %in% names(fus)) reportedType else ""), rank = 1L)]), fill = TRUE)
som <- rbindlist(list(som, cn_rows, sv_rows), fill = TRUE)
if (nrow(som)) {
  som <- unique(som)[order(rank, Gene)][, rank := NULL]
  som[, `:=`(Annotation = "", Tier = "")]                      # OncoKB / CIViC / ClinVar / COSMIC: knowledge-base step
}
fwrite(som, file.path(od, "drivers_somatic_table.csv"))
write_md(som, file.path(od, "drivers_somatic_table.md"),
         "Driver likelihood from PURPLE driver catalogue; Callers: SNVs out of 4 (Mutect2, Strelka2, MuSE, SAGE), indels out of 3; annotation sources: [OncoKB / CIViC / ClinVar / COSMIC]")
message("somatic table: ", nrow(som), " rows"); print(som)

# germline: PURPLE germline catalogue + its germline VCF (REPORTED or in a catalogue gene)
germ_cat <- if (nrow(cat_all)) cat_all[source == "PURPLE germline"] else data.table()
gv_f <- find_one(pd, "\\.purple\\.germline\\.vcf\\.gz$", required = FALSE)
germ <- data.table()
if (!is.na(gv_f)) {
  gv <- read_vcf_table(gv_f)[FILTER == "PASS" & CHROM %chin% STD_CHR]
  gi <- if ("info_IMPACT" %in% names(gv)) tstrsplit(as.character(gv$info_IMPACT), ",", fixed = TRUE, fill = NA_character_) else list()
  gv[, `:=`(GENE = if (length(gi)) gi[[1]] else NA_character_, EFF = if (length(gi) >= 3) gi[[3]] else NA_character_,
            HGVSC = if (length(gi) >= 6) gi[[6]] else NA_character_, HGVSP = if (length(gi) >= 7) gi[[7]] else NA_character_,
            REP = if ("info_REPORTED" %in% names(gv)) as.character(info_REPORTED) %chin% c("TRUE", "true", "1") else FALSE)]
  if (!"n_GT" %in% names(gv)) gv[, n_GT := NA_character_]
  info_chr <- function(k) if (k %in% names(gv)) as.character(gv[[k]]) else rep(NA_character_, nrow(gv))
  gv[, `:=`(CLN = info_chr("info_CLNSIG"), MACN = suppressWarnings(as.numeric(info_chr("info_PURPLE_MACN"))),
            BIAL = info_chr("info_BIALLELIC") %chin% c("TRUE", "true", "1"))]
  gsel <- gv[REP | GENE %chin% germ_cat$gene]
  if (nrow(gsel)) germ <- gsel[, .(
    Gene = GENE,
    `Variant (c. / p.)` = paste0(fcoalesce(HGVSC, paste0(REF, ">", ALT1)), fifelse(!is.na(HGVSP) & HGVSP != "", paste0(" / ", HGVSP), "")),
    Consequence = gsub("_", " ", EFF),
    Zygosity = fcase(n_GT %chin% c("1/1", "1|1"), "homozygous", n_GT %chin% c("0/1", "0|1", "1|0"), "heterozygous", default = as.character(n_GT)),
    `VAF normal → tumour` = sprintf("%.2f → %.2f", N_VAF, VAF),
    `Tumour 2nd hit` = fcase(BIAL | (GENE %chin% germ_cat[as.character(biallelic) %chin% c("true", "TRUE"), gene]), "biallelic (PURPLE)",
                             !is.na(MACN) & MACN < 0.5, sprintf("LOH (minor allele CN %.1f)", MACN),
                             VAF > N_VAF + 0.1, "allelic imbalance, variant allele gained", default = "none detected"),
    ClinVar = fcoalesce(gsub("_", " ", CLN), ""),
    `Driver likelihood` = fmt_lik(germ_cat$driverLikelihood[match(GENE, germ_cat$gene)]),
    `Gene–disease link` = "", Source = "PURPLE germline (Hartwig panel)")]
}

# SPRTN (not on Hartwig's germline panel): Y117C and c.718_718+3del, counted from the SPRTN slice
# BAMs (09) with the same pileup as 11_sprtn_checks.sh; ClinVar from sprtn_clinvar_variants.tsv;
# phasing verdicts from 12_whatshap_sprtn.sh; SPRTN copy number from PURPLE
WGS_SCRIPTS <- Filter(dir.exists, c("../wgs", "scripts/wgs", "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/wgs"))[1]
sl <- file.path(OUT, "igv_slices", "SPRTN_gene")
pile <- function(sample, pos) {
  bam <- file.path(sl, paste0(sample, "_sarek.slice.bam"))
  if (!file.exists(bam)) return(NA_character_)
  out <- tryCatch(system2("samtools", c("mpileup", "-Q", "20", "-q", "20", "-f", shQuote(FASTA), "-r", sprintf("chr1:%d-%d", pos, pos), shQuote(bam)),
                          stdout = TRUE, stderr = FALSE), error = function(e) character(0))
  if (!length(out)) return(NA_character_)
  strsplit(out[1], "\t", fixed = TRUE)[[1]]                     # chrom, pos, ref, depth, bases, quals
}
# as 11_sprtn_checks.sh: Y117C = G / (ref + G); deletion = reads with a deletion / depth
y117 <- function(sample) { f <- pile(sample, 231347825L); if (anyNA(f)) return(NA_real_)
  b <- toupper(gsub("\\^.|\\$", "", f[5])); r <- nchar(gsub("[^.,]", "", b)); g <- nchar(gsub("[^G]", "", b)); g / (r + g) }
del718 <- function(sample) { f <- pile(sample, 231351569L); if (anyNA(f)) return(NA_real_)
  n <- lengths(regmatches(f[5], gregexpr("-[0-9]+[ACGTNacgtn]+", f[5]))); n / as.numeric(f[4]) }
cv <- if (!is.na(WGS_SCRIPTS) && file.exists(file.path(WGS_SCRIPTS, "sprtn_clinvar_variants.tsv")))
  fread(file.path(WGS_SCRIPTS, "sprtn_clinvar_variants.tsv"), skip = "chrom") else data.table()
cv_of <- function(st) if (nrow(cv) && any(cv$start == st)) cv[start == st, sprintf("%s (%s)", classification[1], clinvar[1])] else ""
read_result <- function(f) if (file.exists(f)) { r <- grep("^# result:", readLines(f), value = TRUE); if (length(r)) sub("^# result: *", "", tail(r, 1)) else NA_character_ } else NA_character_
ph_dir <- file.path(OUT, "phasing", "SPRTN")
ph_imb <- read_result(file.path(ph_dir, "haplotype_imbalance.txt"))
ph_wh  <- read_result(file.path(ph_dir, "normal_tumour", "result.txt"))
in_trans <- !is.na(ph_imb) && grepl("IN TRANS", ph_imb, ignore.case = TRUE)
cng_f <- find_one(pd, "\\.purple\\.cnv\\.gene\\.tsv$", required = FALSE)
sprtn_cn <- if (!is.na(cng_f)) fread(cng_f)[gene == "SPRTN"][1] else data.table()
cn_txt <- if (nrow(sprtn_cn) && all(c("minCopyNumber", "minMinorAlleleCopyNumber") %in% names(sprtn_cn)))
  sprintf("CN %.1f, minor allele %.1f", sprtn_cn$minCopyNumber, sprtn_cn$minMinorAlleleCopyNumber) else "CN: see PURPLE"
vaf_txt <- function(f) { n <- f(NORMAL); t <- f(TUMOUR); if (is.na(n) || is.na(t)) "see 11_sprtn_checks.sh" else sprintf("%.2f → %.2f", n, t) }
zyg <- if (in_trans) "compound heterozygous (in trans)" else "heterozygous (phase: see 12)"
sprtn <- data.table(
  Gene = "SPRTN",
  `Variant (c. / p.)` = c("c.350A>G / p.Tyr117Cys", "c.718_718+3del"),
  Consequence = c("missense", "splice donor deletion"),
  Zygosity = zyg,
  `VAF normal → tumour` = c(vaf_txt(y117), vaf_txt(del718)),
  `Tumour 2nd hit` = c(paste0(cn_txt, if (in_trans) "; Y117C allele gained (allelic imbalance)" else "; phasing: see 12"),
                       paste0(cn_txt, if (in_trans) "; deletion allele on the less-amplified copy" else "")),
  ClinVar = c(cv_of(231347825L), cv_of(231351571L)),
  `Driver likelihood` = "— (not on Hartwig panel)",
  `Gene–disease link` = "Ruijs-Aalfs syndrome (progeroid features, early-onset HCC); biallelic SPRTN loss",
  Source = "own analysis (11, 12; not on Hartwig germline panel)")
germ <- rbindlist(list(sprtn, germ), fill = TRUE)
fwrite(germ, file.path(od, "drivers_germline_table.csv"))
write_md(germ, file.path(od, "drivers_germline_table.md"),
         paste0("Germline: PURPLE germline driver catalogue (Hartwig germline panel) + SPRTN (own analysis). VAF normal → tumour from read counts; ",
                "WhatsHap: ", fcoalesce(ph_wh, "not run"), ". Research findings - confirm in an accredited lab before any clinical use."))
message("germline table: ", nrow(germ), " rows"); print(germ)

message("done: ", od)
