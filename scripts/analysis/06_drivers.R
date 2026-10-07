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
message("done: ", od)
