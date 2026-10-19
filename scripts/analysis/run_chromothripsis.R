#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# run_chromothripsis.R - chromothripsis in one tumour/normal pair with ShatterSeek, called twice
# and cross-validated:
#   A  oncoanalyser: ESVEE SVs as PURPLE wrote them (purple.sv.vcf.gz) + PURPLE CN   (primary)
#   B  sarek:        Manta somatic SVs + ASCAT CN                                       (check)
# Then per chromosome called in either set: A vs B (tier, cluster, oscillations, boundaries),
# ESVEE-Manta breakpoint concordance, LOH in low vs high CN states, purity/ploidy, overlapping
# LINX clusters, breakends in centromeres / segmental duplications; plots and a short report.
#
#   Rscript run_chromothripsis.R --oncoanalyser_dir <.../oncoanalyser/RJALS/RJALS> \
#     --sarek_dir <.../sarek/RJALS> --tumor_id RJALS_Tm --normal_id RJALS_N --outdir <dir>
#
# Env: env_chromothripsis.yml (+ ShatterSeek from GitHub). Reads pipeline outputs only; outdir
# must lie outside both pipeline folders and outside cloud-synced folders (patient data).
# Step 1 (inputs, contig naming, columns): ../wgs/14_chromothripsis_inspect.sh
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(ggplot2)
  library(GenomicRanges); library(VariantAnnotation); library(StructuralVariantAnnotation)
  library(ShatterSeek); library(gridExtra); library(grid); library(circlize)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--oncoanalyser_dir", help = "oncoanalyser sample folder (contains purple/, linx/)"),
  make_option("--sarek_dir",        help = "sarek patient folder (contains variant_calling/)"),
  make_option("--tumor_id"), make_option("--normal_id"), make_option("--outdir"),
  make_option("--fragment_joins", default = "fdr",
              help = "fdr: BH q > 0.2 (Supplementary Note) [default]; raw: p > 0.05"),
  make_option("--match_bp", type = "integer", default = 500L, help = "breakpoint match window [500]"),
  make_option("--ref_dir", help = "UCSC hg38 centromeres.txt.gz + genomicSuperDups.txt.gz (downloaded if missing) [<outdir>/ref]"),
  make_option("--container_dir", default = Sys.getenv("NXF_SINGULARITY_CACHEDIR", "/common/WORK/pstancl/singularity_cache"),
              help = "Singularity image cache with oncoanalyser's LINX image"),
  make_option("--hmf_ref_dir", default = "/common/WORK/pstancl/references/hmf/oncoanalyser",
              help = "HMF reference data (searched for ensembl_data/ for the LINX command)")
)))
for (k in c("oncoanalyser_dir", "sarek_dir", "tumor_id", "normal_id", "outdir"))
  if (is.null(opt[[k]])) stop("missing --", k, call. = FALSE)
stopifnot(opt$fragment_joins %in% c("fdr", "raw"))
ONCO <- normalizePath(opt$oncoanalyser_dir, mustWork = TRUE)
SAREK <- normalizePath(opt$sarek_dir, mustWork = TRUE)
TUM <- opt$tumor_id; NOR <- opt$normal_id
OD <- normalizePath(opt$outdir, mustWork = FALSE)
if (grepl("OneDrive|CloudStorage|Dropbox|iCloud|Google Drive", OD, ignore.case = TRUE))
  stop("Refusing to write patient data to a cloud-synced folder: ", OD, call. = FALSE)
if (startsWith(OD, ONCO) || startsWith(OD, SAREK))
  stop("--outdir must be outside the pipeline output folders (they are not modified)", call. = FALSE)
PLOTS <- file.path(OD, "plots"); dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
REF <- if (is.null(opt$ref_dir)) file.path(OD, "ref") else opt$ref_dir
dir.create(REF, recursive = TRUE, showWarnings = FALSE)
TOL <- opt$match_bp
CHROMS <- c(1:22, "X")                     # ShatterSeek: chromosomes 1-22 and X, no "chr"
SET_LAB <- c(A = "A oncoanalyser (ESVEE/PURPLE)", B = "B sarek (Manta/ASCAT)")
log_lines <- character()
note <- function(...) { m <- paste0(...); message(m); log_lines <<- c(log_lines, m) }

# ---- 1. inputs (found, not assumed) -------------------------------------------------------
find_one <- function(dir, pattern, label, required = TRUE) {
  f <- sort(list.files(dir, pattern, recursive = TRUE, full.names = TRUE))
  if (!length(f)) {
    if (required) stop(label, " not found under ", dir, " (", pattern, ")", call. = FALSE)
    note("MISSING ", label, " (", pattern, ")"); return(NA_character_)
  }
  if (length(f) > 1) note(label, ": ", length(f), " matches, using ", f[1])
  f[1]
}
esc <- function(x) gsub(".", "\\.", x, fixed = TRUE)   # sample IDs: letters, digits, _ and .
pair <- paste0(esc(TUM), "_vs_", esc(NOR))
F <- list(
  purple_sv  = find_one(ONCO, paste0("^", esc(TUM), "\\.purple\\.sv\\.vcf\\.gz$"), "PURPLE SV VCF"),
  purple_cn  = find_one(ONCO, paste0("^", esc(TUM), "\\.purple\\.cnv\\.somatic\\.tsv$"), "PURPLE CN"),
  purple_pur = find_one(ONCO, paste0("^", esc(TUM), "\\.purple\\.purity\\.tsv$"), "PURPLE purity"),
  linx_vis   = find_one(ONCO, paste0("^", esc(TUM), "\\.linx\\.vis_sv_data\\.tsv$"), "LINX vis_sv_data"),
  linx_clu   = find_one(ONCO, paste0("^", esc(TUM), "\\.linx\\.clusters\\.tsv$"), "LINX clusters"),
  manta      = find_one(SAREK, paste0("^", pair, "\\.manta\\.(somatic_sv|somaticSV)\\.vcf\\.gz$"), "Manta somatic SV VCF"),
  ascat_seg  = find_one(SAREK, paste0("^", pair, "\\.segments\\.txt$"), "ASCAT segments"),
  ascat_pp   = find_one(SAREK, paste0("^", pair, "\\.purityploidy\\.txt$"), "ASCAT purity/ploidy")
)
fwrite(data.table(input = names(F), path = unlist(F)), file.path(OD, "inputs_used.tsv"), sep = "\t")

need_cols <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) stop(label, ": missing column(s) ", paste(miss, collapse = ", "), call. = FALSE)
}
norm_chr <- function(x) sub("^chr", "", as.character(x))

# both VCFs must be GRCh38 with the same contig naming; chr1 length tells hg38 from hg19
for (k in c("purple_sv", "manta")) {
  si <- seqinfo(scanVcfHeader(F[[k]]))
  c1 <- intersect(c("chr1", "1"), seqnames(si))
  len1 <- if (length(c1)) seqlengths(si)[c1[1]] else NA
  note(k, ": contigs ", if (any(startsWith(seqnames(si), "chr"))) "chr-prefixed" else "without chr",
       ", chr1 length ", len1, if (identical(unname(len1), 248956422L)) " (GRCh38)" else " (NOT GRCh38?)")
}

# ---- 2. SVs -> one record per SV, ShatterSeek types from breakend strands ------------------
# StructuralVariantAnnotation strands match ShatterSeek's convention (a deletion is +/-):
# +/- DEL, -/+ DUP, +/+ h2hINV, -/- t2tINV, different chromosomes TRA.
sv_type <- function(c1, s1, c2, s2)
  fifelse(c1 != c2, "TRA", fcase(s1 == "+" & s2 == "-", "DEL", s1 == "-" & s2 == "+", "DUP",
                                 s1 == "+" & s2 == "+", "h2hINV", default = "t2tINV"))
# lower breakend first (intrachromosomal pos1 < pos2, which ShatterSeek needs), keep 1-22, X
tidy_svs <- function(d, label) {
  d[, `:=`(chr1 = norm_chr(chr1), chr2 = norm_chr(chr2))]
  n0 <- nrow(d)
  d <- d[chr1 %chin% CHROMS & chr2 %chin% CHROMS]
  swap <- d[, match(chr1, CHROMS) > match(chr2, CHROMS) | (chr1 == chr2 & pos1 > pos2)]
  d[swap, `:=`(chr1 = chr2, pos1 = pos2, strand1 = strand2, chr2 = chr1, pos2 = pos1, strand2 = strand1)]
  d[, type := sv_type(chr1, strand1, chr2, strand2)]
  note(label, ": ", n0, " SVs, ", nrow(d), " on 1-22/X: ",
       paste(d[, .N, by = type][order(type), paste0(type, " ", N)], collapse = ", "))
  d[]
}
read_sv_vcf <- function(path, label) {
  vcf <- readVcf(path, "hg38")
  flt <- rowRanges(vcf)$FILTER
  vcf <- vcf[flt %in% c("PASS", ".")]
  gr <- breakpointRanges(vcf)                      # breakpoints only: single breakends are dropped
  bp <- data.table(id = names(gr), partner = gr$partner, chr = as.character(seqnames(gr)),
                   pos = as.integer((start(gr) + end(gr)) / 2), strand = as.character(strand(gr)))
  d <- bp[id < partner][bp, on = .(partner = id), nomatch = NULL,
                         .(id, chr1 = chr, pos1 = pos, strand1 = strand,
                           chr2 = i.chr, pos2 = i.pos, strand2 = i.strand)]
  note(label, ": ", length(flt), " records, ", sum(flt %in% c("PASS", ".")), " PASS -> ", nrow(d), " breakpoint SVs")
  tidy_svs(d, label)
}
sv <- list(A = read_sv_vcf(F$purple_sv, "A PURPLE SV VCF (ESVEE)"),
           B = read_sv_vcf(F$manta, "B Manta"))

# LINX vis_sv_data: orientation 1 -> "+", -1 -> "-"; single / inferred breakends left out
vis <- fread(F$linx_vis)
need_cols(vis, c("ClusterId", "SvId", "Type", "ResolvedType", "ChrStart", "ChrEnd", "PosStart",
                 "PosEnd", "OrientStart", "OrientEnd"), "LINX vis_sv_data")
vis[, `:=`(ChrStart = as.character(ChrStart), ChrEnd = as.character(ChrEnd))]
vis_bp <- vis[!Type %chin% c("SGL", "INF") & !ChrEnd %chin% c("-1", "0")]
linx_sv <- tidy_svs(vis_bp[, .(id = as.character(SvId), chr1 = ChrStart, pos1 = PosStart,
                               strand1 = fifelse(OrientStart == 1, "+", "-"), chr2 = ChrEnd, pos2 = PosEnd,
                               strand2 = fifelse(OrientEnd == 1, "+", "-"))], "LINX vis_sv_data")
# does LINX see the same SVs as the PURPLE VCF? breakends per chromosome
bends <- function(d) rbind(d[, .(chr = chr1)], d[, .(chr = chr2)])[, .N, by = chr]
sv_check <- merge(bends(sv$A)[, .(chr, purple_vcf = N)], bends(linx_sv)[, .(chr, linx = N)], by = "chr", all = TRUE)
for (k in c("purple_vcf", "linx")) set(sv_check, which(is.na(sv_check[[k]])), k, 0L)
sv_check[, `:=`(agree = purple_vcf == linx, ord = match(chr, CHROMS))]
setorder(sv_check, ord)[, ord := NULL]
fwrite(sv_check, file.path(OD, "check_linx_vs_purple_vcf_breakends.tsv"), sep = "\t")
note("LINX vs PURPLE SV VCF breakends per chromosome: ", sv_check[agree == TRUE, .N], "/", nrow(sv_check),
     " chromosomes agree", if (any(!sv_check$agree))
       paste0(" (differ: ", paste(sv_check[agree == FALSE, paste0(chr, " ", purple_vcf, "/", linx)], collapse = ", "), ")"))

# ---- 3. copy number ------------------------------------------------------------------------
pcn <- fread(F$purple_cn)
need_cols(pcn, c("chromosome", "start", "end", "copyNumber", "minorAlleleCopyNumber"), "PURPLE CN")
asc <- fread(F$ascat_seg)
need_cols(asc, c("chr", "startpos", "endpos", "nMajor", "nMinor"), "ASCAT segments")
# raw segments (kept for LOH and plots): total = rounded CN clamped at 0, LOH from the minor allele
seg <- list(
  A = pcn[, .(chr = norm_chr(chromosome), start, end, cn_raw = copyNumber, minor = minorAlleleCopyNumber,
              total_cn = pmax(round(copyNumber), 0), loh = minorAlleleCopyNumber < 0.5)],
  B = asc[, .(chr = norm_chr(chr), start = startpos, end = endpos, cn_raw = nMajor + nMinor, minor = nMinor,
              total_cn = nMajor + nMinor, loh = nMinor == 0)])
seg <- lapply(seg, function(x) setorder(x[chr %chin% CHROMS], chr, start))
# segments without a CN (ASCAT can leave nMajor/nMinor empty) are dropped, as in the ShatterSeek README
for (k in names(seg)) if (any(is.na(seg[[k]]$total_cn))) {
  note(k, " CN: dropped ", sum(is.na(seg[[k]]$total_cn)), " segment(s) without a copy number")
  seg[[k]] <- seg[[k]][!is.na(total_cn)]
}
# ShatterSeek wants adjacent segments with different CN: merge runs of equal total CN
merge_cn <- function(x) x[, .(start = min(start), end = max(end), total_cn = total_cn[1]),
                          by = .(chr, run = rleid(chr, total_cn))][, run := NULL][]
cn <- lapply(seg, merge_cn)
for (k in names(cn)) note(k, " CN: ", nrow(seg[[k]]), " segments -> ", nrow(cn[[k]]), " after merging equal total CN")

# purity / ploidy from both callers
pp_p <- fread(F$purple_pur); need_cols(pp_p, c("purity", "ploidy"), "PURPLE purity")
pp_a <- fread(F$ascat_pp);   need_cols(pp_a, c("AberrantCellFraction", "Ploidy"), "ASCAT purity/ploidy")
pp <- data.table(caller = c("PURPLE", "ASCAT"), purity = c(pp_p$purity[1], pp_a$AberrantCellFraction[1]),
                 ploidy = c(pp_p$ploidy[1], pp_a$Ploidy[1]))
# my thresholds for "substantial": 0.10 purity, 0.3 ploidy (about a tenth of a copy per 3n genome)
pp_flag <- abs(diff(pp$purity)) > 0.10 || abs(diff(pp$ploidy)) > 0.3
note("purity/ploidy PURPLE ", pp$purity[1], "/", round(pp$ploidy[1], 2), ", ASCAT ", pp$purity[2], "/",
     round(pp$ploidy[2], 2), if (pp_flag) " - DISAGREE" else " - agree")

# ---- 4. ShatterSeek + Cortes-Ciriano et al. 2020 criteria -----------------------------------
run_ss <- function(s, c) shatterseek(
  SV.sample = SVs(chrom1 = s$chr1, pos1 = as.numeric(s$pos1), chrom2 = s$chr2, pos2 = as.numeric(s$pos2),
                  SVtype = s$type, strand1 = s$strand1, strand2 = s$strand2),
  seg.sample = CNVsegs(chrom = c$chr, start = c$start, end = c$end, total_cn = c$total_cn),
  min.Size = 1, genome = "hg38")
ss <- list(A = run_ss(sv$A, cn$A), B = run_ss(sv$B, cn$B))

# Criteria - ShatterSeek tutorial (inst/tutorial/tutorial.pdf, section 3 "Recommended cut-off
# values") = Cortes-Ciriano et al., Nat Genet 2020, 52:331-341, Supplementary Note p. 22:
#   High: >= 6 interleaved intrachromosomal SVs, >= 7 adjacent segments oscillating between 2 CN
#         states, the fragment joins test, and either the chromosomal enrichment or the
#         exponential distribution of breakpoints test.
#   High: >= 3 interleaved intrachromosomal SVs and >= 4 interchromosomal SVs, >= 7 adjacent
#         segments oscillating between 2 CN states and the fragment joins test.
#   Low:  >= 6 interleaved intrachromosomal SVs, 4, 5 or 6 adjacent segments oscillating between
#         2 CN states, the fragment joins test, and either the chromosomal enrichment or the
#         exponential distribution of breakpoints test.
# Tests (Supplementary Note p. 22): "The P values from these three tests were corrected using the
# FDR method. The level of significance was set to q <= 0.2." Here BH across the chromosomes of
# this sample (the paper corrected across its cohort). The fragment joins test is passed when the
# SV types do NOT differ from equal proportions: q > 0.2 (--fragment_joins raw: p > 0.05).
# Columns: clusterSize = interleaved intrachromosomal SVs (sum of the four number_* by type, see
# call_tiers), number_TRA = interchromosomal SVs in the
# cluster, max_number_oscillating_CN_segments_2_states, pval_exp_cluster = exponential test on the
# cluster's breakpoints, chr_breakpoint_enrichment = enrichment test.
bh <- function(p) { q <- rep(NA_real_, length(p)); ok <- !is.na(p); q[ok] <- p.adjust(p[ok], "BH"); q }
yes <- function(x) !is.na(x) & x
call_tiers <- function(o, set) {
  s <- as.data.table(o@chromSummary)
  # statistical_criteria() replaces chromSummary and drops clusterSize; its number_DEL / DUP /
  # h2hINV / t2tINV count the intrachromosomal SVs of the interleaved cluster (connComp) by type
  need_cols(s, c("number_DEL", "number_DUP", "number_h2hINV", "number_t2tINV"), "ShatterSeek chromSummary")
  s[, clusterSize := number_DEL + number_DUP + number_h2hINV + number_t2tINV]
  need_cols(s, c("chrom", "start", "end", "clusterSize", "number_TRA", "pval_fragment_joins",
                 "chr_breakpoint_enrichment", "pval_exp_cluster", "max_number_oscillating_CN_segments_2_states"),
            "ShatterSeek chromSummary")
  s[, `:=`(q_fragment_joins = bh(pval_fragment_joins), q_chr_enrichment = bh(chr_breakpoint_enrichment),
           q_exp_cluster = bh(pval_exp_cluster))]
  s[, fragment_joins_ok := yes(if (opt$fragment_joins == "fdr") q_fragment_joins > 0.2 else pval_fragment_joins > 0.05)]
  s[, enrich_or_exp_ok := yes(q_chr_enrichment <= 0.2) | yes(q_exp_cluster <= 0.2)]
  osc <- s$max_number_oscillating_CN_segments_2_states
  s[, `:=`(high_1 = yes(clusterSize >= 6) & yes(osc >= 7) & fragment_joins_ok & enrich_or_exp_ok,
           high_2 = yes(clusterSize >= 3) & yes(number_TRA >= 4) & yes(osc >= 7) & fragment_joins_ok,
           low    = yes(clusterSize >= 6) & yes(osc %in% 4:6) & fragment_joins_ok & enrich_or_exp_ok)]
  s[, call := fcase(high_1 | high_2, "high", low, "low", default = "none")]
  s[, chrom := as.character(chrom)]
  cbind(set = set, input = SET_LAB[[set]], s)
}
summ <- rbindlist(lapply(names(ss), function(k) call_tiers(ss[[k]], k)), fill = TRUE)
summ[, ord := match(chrom, CHROMS)]; setorder(summ, set, ord); summ[, ord := NULL]
fwrite(summ, file.path(OD, "shatterseek_summary.tsv"), sep = "\t")
called <- summ[call != "none", unique(chrom)]
called <- called[order(match(called, CHROMS))]
note("called (either set): ", if (length(called)) paste(summ[call != "none", paste0(set, ":chr", chrom, " ", call)], collapse = ", ") else "none")
# nothing called: still look at the largest clusters (>= 6 interleaved SVs) as candidates
cand <- if (length(called)) called else summ[yes(clusterSize >= 6), unique(chrom)]

# ---- 5. cross-validation and annotation per chromosome ---------------------------------------
# centromeres + segmental duplications (UCSC hg38) for the artefact flag
ucsc <- function(f) {
  p <- file.path(REF, f)
  if (!file.exists(p)) tryCatch(download.file(paste0("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/", f),
                                              p, quiet = TRUE, mode = "wb"),
                                error = function(e) note("could not download ", f, " - artefact flag incomplete"))
  if (file.exists(p)) fread(p, header = FALSE, select = 2:4, col.names = c("chr", "start", "end"))[
    , chr := norm_chr(chr)][chr %chin% CHROMS] else NULL
}
artef <- list(centromere = ucsc("centromeres.txt.gz"), segdup = ucsc("genomicSuperDups.txt.gz"))
in_bed <- function(cc, pos, bed) {
  if (is.null(bed)) return(rep(NA, length(cc)))
  overlapsAny(GRanges(cc, IRanges(pos, width = 1)), GRanges(bed$chr, IRanges(bed$start + 1L, bed$end)))
}
for (k in names(sv)) sv[[k]][, `:=`(
  centromere = in_bed(chr1, pos1, artef$centromere) | in_bed(chr2, pos2, artef$centromere),
  segdup     = in_bed(chr1, pos1, artef$segdup)     | in_bed(chr2, pos2, artef$segdup))]

# one SV matches another when both breakends lie within TOL bp with the same strands
near <- function(c1, p1, s1, c2, p2, s2, b)
  b$chr1 == c1 & abs(b$pos1 - p1) <= TOL & b$strand1 == s1 & b$chr2 == c2 & abs(b$pos2 - p2) <= TOL & b$strand2 == s2
match_both <- function(a, b) if (!nrow(a) || !nrow(b)) rep(FALSE, nrow(a)) else
  vapply(seq_len(nrow(a)), function(i) a[i, any(near(chr1, pos1, strand1, chr2, pos2, strand2, b) |
                                              near(chr2, pos2, strand2, chr1, pos1, strand1, b))], logical(1))
match_one <- function(a, b) if (!nrow(a) || !nrow(b)) rep(FALSE, nrow(a)) else {
  bb <- rbind(b[, .(c = chr1, p = pos1, s = strand1)], b[, .(c = chr2, p = pos2, s = strand2)])
  hit <- function(c0, p0, s0) any(bb$c == c0 & abs(bb$p - p0) <= TOL & bb$s == s0)
  vapply(seq_len(nrow(a)), function(i) a[i, hit(chr1, pos1, strand1) || hit(chr2, pos2, strand2)], logical(1))
}
# NA (not 0) when the UCSC tracks could not be read
nflag <- function(x) if (!length(x) || all(is.na(x))) NA_integer_ else sum(x, na.rm = TRUE)
in_region <- function(d, cc, s0, e0) d[(chr1 == cc & pos1 >= s0 & pos1 <= e0) | (chr2 == cc & pos2 >= s0 & pos2 <= e0)]

# LOH: the two most common total-CN states (by bp) in the region; chromothripsis typically has
# LOH in the low state and retained heterozygosity in the high state (Korbel & Campbell 2013)
loh_region <- function(x, cc, s0, e0) {
  x <- x[chr == cc & end >= s0 & start <= e0][, w := pmin(end, e0) - pmax(start, s0) + 1]
  st <- x[, .(bp = sum(w)), by = total_cn][order(-bp)]
  if (nrow(st) < 2) return(list(states = if (nrow(st)) as.character(st$total_cn[1]) else NA_character_,
                                low_state_loh_frac = NA_real_, high_state_het_frac = NA_real_))
  lo <- min(st$total_cn[1:2]); hi <- max(st$total_cn[1:2])
  list(states = paste0(lo, "/", hi),
       low_state_loh_frac  = round(x[total_cn == lo, sum(w[yes(loh)]) / sum(w)], 2),
       high_state_het_frac = round(x[total_cn == hi, sum(w[yes(!loh)]) / sum(w)], 2))
}

clu <- fread(F$linx_clu)
need_cols(clu, c("clusterId", "resolvedType", "clusterCount"), "LINX clusters")
vis[, `:=`(c1 = norm_chr(ChrStart), c2 = norm_chr(ChrEnd))]
linx_clusters <- function(cc, s0, e0) {
  ids <- vis[(c1 == cc & PosStart >= s0 & PosStart <= e0) | (c2 == cc & PosEnd >= s0 & PosEnd <= e0), unique(ClusterId)]
  if (!length(ids)) return(data.table())
  out <- clu[clusterId %in% ids, .(clusterId, resolvedType, clusterCount)]
  out[, chromosomes := vapply(clusterId, function(i) {
    ch <- vis[ClusterId == i, unique(c(c1, c2))]; ch <- ch[ch %chin% c(CHROMS, "Y")]
    paste(ch[order(match(ch, c(CHROMS, "Y")))], collapse = ",") }, "")]
  out[, other_chromosomes := vapply(strsplit(chromosomes, ","), function(ch) paste(setdiff(ch, cc), collapse = ","), "")]
  out[order(-clusterCount)]
}

reg <- function(k, cc) summ[set == k & chrom == cc]
conc <- rbindlist(lapply(cand, function(cc) {
  a <- reg("A", cc); b <- reg("B", cc)
  has <- function(r) nrow(r) && !is.na(r$start)
  s0 <- min(c(if (has(a)) a$start, if (has(b)) b$start)); e0 <- max(c(if (has(a)) a$end, if (has(b)) b$end))
  ov <- if (has(a) && has(b)) max(0, min(a$end, b$end) - max(a$start, b$start)) else NA_real_
  ra <- in_region(sv$A, cc, s0, e0); rb <- in_region(sv$B, cc, s0, e0)
  la <- loh_region(seg$A, cc, s0, e0); lb <- loh_region(seg$B, cc, s0, e0)
  lx <- linx_clusters(cc, s0, e0)
  data.table(
    chrom = cc, region_start = s0, region_end = e0,
    call_A = a$call, call_B = b$call,
    clusterSize_A = a$clusterSize, clusterSize_B = b$clusterSize,
    number_TRA_A = a$number_TRA, number_TRA_B = b$number_TRA,
    osc2_A = a$max_number_oscillating_CN_segments_2_states, osc2_B = b$max_number_oscillating_CN_segments_2_states,
    start_A = a$start, end_A = a$end, start_B = b$start, end_B = b$end,
    overlap_bp = ov, jaccard = if (is.na(ov)) NA_real_ else round(ov / (max(a$end, b$end) - min(a$start, b$start)), 3),
    sv_A_in_region = nrow(ra), sv_B_in_region = nrow(rb),
    A_matched_both_ends = sum(match_both(ra, rb)), A_matched_one_end = sum(match_one(ra, rb)),
    B_matched_both_ends = sum(match_both(rb, ra)),
    A_centromere = nflag(ra$centromere), A_segdup = nflag(ra$segdup),
    B_centromere = nflag(rb$centromere), B_segdup = nflag(rb$segdup),
    cn_states_A = la$states, low_state_loh_A = la$low_state_loh_frac, high_state_het_A = la$high_state_het_frac,
    cn_states_B = lb$states, low_state_loh_B = lb$low_state_loh_frac, high_state_het_B = lb$high_state_het_frac,
    purity_purple = pp$purity[1], ploidy_purple = pp$ploidy[1], purity_ascat = pp$purity[2], ploidy_ascat = pp$ploidy[2],
    purity_ploidy_disagree = pp_flag,
    linx_clusters = if (nrow(lx)) paste(lx$clusterId, collapse = ",") else "",
    linx_resolvedType = if (nrow(lx)) paste(lx$resolvedType, collapse = ",") else "",
    linx_clusterCount = if (nrow(lx)) paste(lx$clusterCount, collapse = ",") else "",
    linx_chromosomes = if (nrow(lx)) paste(lx$chromosomes, collapse = ";") else "",
    linx_other_chromosomes = if (nrow(lx)) paste(unique(unlist(strsplit(lx$other_chromosomes, ","))), collapse = ",") else "")
}), fill = TRUE)
fwrite(conc, file.path(OD, "concordance.tsv"), sep = "\t")
# every SV in the candidate regions with its match and artefact flags, for a look by eye
sv_detail <- rbindlist(lapply(cand, function(cc) {
  r <- conc[chrom == cc]
  ra <- in_region(sv$A, cc, r$region_start, r$region_end); rb <- in_region(sv$B, cc, r$region_start, r$region_end)
  rbind(cbind(set = "A", ra, matched_both_ends = match_both(ra, rb), matched_one_end = match_one(ra, rb)),
        cbind(set = "B", rb, matched_both_ends = match_both(rb, ra), matched_one_end = match_one(rb, ra)))[, region_chrom := cc]
}), fill = TRUE)
if (nrow(sv_detail)) fwrite(sv_detail, file.path(OD, "concordance_svs.tsv"), sep = "\t")
# breakends of SVs that only one caller has in a candidate region, as IGV loci (+/-300 bp) for
#   REGIONS=<this file> SOFT_CLIPS=true OUT_IGV=<...>/igv/chromothripsis ../wgs/05_igv_snapshots.sh
only <- if (nrow(sv_detail)) sv_detail[matched_both_ends == FALSE] else sv_detail
if (nrow(only)) {
  locus <- function(cc, p) sprintf("chr%s:%d-%d", cc, pmax(1L, as.integer(p) - 300L), as.integer(p) + 300L)
  only[, tag := sprintf("%s_only_chr%s_%s_chr%s_%d_chr%s_%d", set, region_chrom, type, chr1, as.integer(pos1), chr2, as.integer(pos2))]
  igv <- rbind(only[, .(name = paste0(tag, "_end1"), locus = locus(chr1, pos1))],
               only[, .(name = paste0(tag, "_end2"), locus = locus(chr2, pos2))])
  writeLines(c("# SVs in candidate regions found by one caller only (run_chromothripsis.R); name <TAB> locus",
               paste(igv$name, igv$locus, sep = "\t")), file.path(OD, "igv_regions_unmatched.tsv"))
  note("IGV loci for ", nrow(only), " single-caller SVs: igv_regions_unmatched.tsv (",
       paste(only[, .N, by = .(set, region_chrom)][, paste0(set, " chr", region_chrom, " ", N)], collapse = ", "), ")")
}

# ---- 6. plots ---------------------------------------------------------------------------------
save_both <- function(p, name, w, h) {
  ggsave(file.path(PLOTS, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(PLOTS, paste0(name, ".png")), p, width = w, height = h, dpi = 150, bg = "white")
}
# ShatterSeek's own plot (SVs, CN, region), set A next to set B. plot_chromothripsis() sets its
# x-axis from ALL intrachromosomal SVs on the chromosome (plot_function.R, min/max of pos1/pos2),
# so one SV far from the cluster zooms a panel out to the whole chromosome. The main plot keeps
# only SVs with both breakends in the shared window (union of the A and B regions), so A and B
# share an x-axis; *_full.png is ShatterSeek's unmodified view.
# plot_chromothripsis()'s table (4th grob) labels a row "Oscillating CN (2 and 3 states)" but shows
# only the 3-state count (plot_function.R cols_sel); the criteria use the 2-state count. Relabel
# rows (and optionally replace their value) in the tableGrob: rows are text grobs in column 1,
# values in column 2.
relabel_table <- function(tg, map) {
  is_txt <- vapply(tg$grobs, inherits, logical(1), "text")
  for (i in which(is_txt)) {
    m <- map[[as.character(tg$grobs[[i]]$label)]]
    if (is.null(m) || tg$layout$l[i] != 1) next
    tg$grobs[[i]]$label <- m[[1]]
    if (!is.null(m[[2]])) {
      j <- which(is_txt & tg$layout$t == tg$layout$t[i] & tg$layout$l == 2)
      if (length(j)) tg$grobs[[j[1]]]$label <- m[[2]]
    }
  }
  tg
}
# replace a text label anywhere inside a grob tree (ShatterSeek's CN panel is labelled "CN")
relabel_grob <- function(g, from, to) {
  if (inherits(g, "text") && identical(as.character(g$label), from)) g$label <- to
  if (!is.null(g$grobs)) g$grobs <- lapply(g$grobs, relabel_grob, from, to)
  if (!is.null(g$children)) g$children <- do.call(gList, lapply(g$children, relabel_grob, from, to))
  g
}
ss_panel <- function(k, cc, win = NULL) tryCatch({
  o <- ss[[k]]
  if (!is.null(win)) {
    s0 <- o@detail$SV
    o@detail$SV <- s0[s0$chrom1 != cc | (s0$pos1 >= win[1] & s0$pos2 <= win[2]), ]
  }
  p <- plot_chromothripsis(ShatterSeek_output = o, chr = paste0("chr", cc),
                           sample_name = paste(TUM, "-", SET_LAB[[k]]), genome = "hg38")
  r <- summ[set == k & chrom == cc]
  p[[4]] <- relabel_table(p[[4]], list(
    "Oscillating CN (2 and 3 states)" = list("Oscillating CN segments, 2 / 3 states",
      paste(r$max_number_oscillating_CN_segments_2_states, "/", r$max_number_oscillating_CN_segments_3_states)),
    "Links with other chrs" = list("Linked SV cluster on other chr (>= 2 transl.)", NULL)))
  # minor-allele CN panel under ShatterSeek's absolute CN, as in Cortes-Ciriano et al. 2020 Fig. 1
  # (ShatterSeek's own BAF= panel needs read counts). Same x range as the package's panels: min
  # pos1 to max pos2 of the intrachromosomal SVs plotted (plot_function.R)
  sc <- o@detail$SV[o@detail$SV$chrom1 == cc, ]
  xr <- if (nrow(sc)) c(min(sc$pos1), max(sc$pos2)) else c(0, 10)
  m <- seg[[k]][chr == cc & end >= xr[1] & start <= xr[2] & !is.na(minor),
                .(start = pmax(start, xr[1]), end = pmin(end, xr[2]), minor = pmax(round(minor), 0))]
  gm <- ggplot(m) + geom_segment(aes(x = start, xend = end, y = minor, yend = minor), linewidth = 2) +
    scale_x_continuous(expand = c(0.01, 0.01), limits = xr, labels = function(x) paste(x / 1e6, "MB")) +
    scale_y_continuous(minor_breaks = NULL, breaks = sort(unique(m$minor)),
                       limits = range(c(0, m$minor)) + c(-0.35, 0.35)) +
    labs(x = NULL, y = "Minor CN") + theme_bw() +    # ShatterSeek's common_ggplot2 sizes (plot_function.R)
    theme(panel.grid = element_blank(), axis.text = element_text(size = 7), axis.title.y = element_text(size = 7),
          axis.line = element_line(linewidth = 0.5, colour = "black"), plot.background = element_blank(),
          plot.margin = unit(c(0, 0.5, 0.2, 0), "cm"))
  g5 <- ggplotGrob(gm)
  if (length(g5$widths) == length(p[[2]]$widths)) g5$widths <- p[[2]]$widths
  p[[3]] <- relabel_grob(p[[3]], "CN", "Absolute CN")
  ttl <- sprintf("%s - %s: chr. %s: %s-%s", TUM, SET_LAB[[k]], cc,
                 format(r$start, big.mark = ","), format(r$end, big.mark = ","))
  arrangeGrob(p[[1]], p[[2]], p[[3]], g5, p[[4]], nrow = 5, ncol = 1, heights = c(0.2, 0.4, 0.35, 0.25, 0.4),
              top = textGrob(ttl, gp = gpar(fontsize = 12)))
}, error = function(e) textGrob(paste0(SET_LAB[[k]], ": no ShatterSeek plot for chr", cc, "\n(", conditionMessage(e), ")")))
for (cc in cand) {
  win <- conc[chrom == cc, c(region_start, region_end)]
  for (k in names(ss)) {
    s0 <- ss[[k]]@detail$SV; out <- sum(s0$chrom1 == cc & !(s0$pos1 >= win[1] & s0$pos2 <= win[2]))
    if (out) note(k, " chr", cc, ": ", out, " intrachromosomal SV(s) outside the plot window ", win[1], "-", win[2], " (shown in _full)")
  }
  save_both(arrangeGrob(ss_panel("A", cc, win), ss_panel("B", cc, win), ncol = 2), paste0("shatterseek_chr", cc), 18, 12)
  save_both(arrangeGrob(ss_panel("A", cc), ss_panel("B", cc), ncol = 2), paste0("shatterseek_chr", cc, "_full"), 18, 12)
}

# total CN and minor-allele CN along the chromosome, both callers, candidate region shaded
for (cc in cand) {
  d <- rbind(seg$A[chr == cc, .(caller = "PURPLE", start, end, total = cn_raw, minor)],
             seg$B[chr == cc, .(caller = "ASCAT", start, end, total = cn_raw, minor)])
  d <- melt(d, measure.vars = c("total", "minor"), variable.name = "cn_type", value.name = "cn")
  r <- conc[chrom == cc]
  p <- ggplot(d) +
    annotate("rect", xmin = r$region_start / 1e6, xmax = r$region_end / 1e6, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "firebrick") +
    geom_segment(aes(x = start / 1e6, xend = end / 1e6, y = cn, yend = cn, colour = cn_type), linewidth = 1.4) +
    facet_wrap(~ caller, ncol = 1) +
    scale_colour_manual(values = c(total = "black", minor = "#1f78b4"), labels = c("total CN", "minor allele CN")) +
    labs(x = paste0("chr", cc, " (Mb)"), y = "copy number", colour = NULL,
         title = paste0(TUM, " chr", cc, ": total and minor-allele copy number"),
         subtitle = "shaded: candidate region (union of A and B); minor CN 0 = LOH") +
    theme_bw() + theme(legend.position = "bottom")
  save_both(p, paste0("cn_loh_chr", cc), 11, 6)
}

# genome-wide circos: called regions (A, B), PURPLE CN, ASCAT CN, set A SV links by type
lens <- seqlengths(seqinfo(scanVcfHeader(F$purple_sv)))
names(lens) <- norm_chr(names(lens)); lens <- lens[CHROMS]
SV_COL <- c(DEL = "#E69F00", DUP = "#0072B2", h2hINV = "#333333", t2tINV = "#009E73", TRA = "#CC79A7")
TIER_COL <- c(high = "firebrick", low = "orange")
circos_draw <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.degree = c(rep(1, length(CHROMS) - 1), 6), cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.004, 0.004))
  circos.genomicInitialize(data.frame(chr = CHROMS, start = 0, end = as.numeric(lens)), plotType = c("axis", "labels"),
                           axis.labels.cex = 0.35, labels.cex = 0.8)
  for (k in c("A", "B")) {
    r <- summ[set == k & call != "none"]
    circos.track(ylim = c(0, 1), track.height = 0.035, bg.border = "grey85", panel.fun = function(x, y) {
      rr <- r[chrom == CELL_META$sector.index]
      if (nrow(rr)) circos.rect(rr$start, 0, rr$end, 1, col = TIER_COL[rr$call], border = NA)
    })
  }
  for (k in c("A", "B")) {
    s <- seg[[k]]; cap <- 8
    circos.track(ylim = c(0, cap), track.height = 0.11, bg.border = "grey85", panel.fun = function(x, y) {
      for (h in c(2, 4, 6)) circos.lines(CELL_META$xlim, c(h, h), col = "grey90")
      ss_ <- s[chr == CELL_META$sector.index]
      if (nrow(ss_)) circos.segments(ss_$start, pmin(ss_$cn_raw, cap), ss_$end, pmin(ss_$cn_raw, cap),
                                     col = if (k == "A") "black" else "#1f78b4", lwd = 1.5)
    })
  }
  for (i in seq_len(nrow(sv$A))) with(sv$A[i], circos.link(chr1, pos1, chr2, pos2, col = adjustcolor(SV_COL[type], 0.7),
                                                         lwd = 0.8, h.ratio = if (chr1 == chr2) 0.3 else 0.6))
  legend("topleft", bty = "n", cex = 0.7, title = "set A SV links",
         legend = names(SV_COL), col = SV_COL, lwd = 2)
  legend("topright", bty = "n", cex = 0.7, legend = c(
    "tracks, outside in:", "called regions A, B (red high / orange low)", "PURPLE total CN (0-8)", "ASCAT total CN (0-8)"))
  title(paste(TUM, "- chromothripsis overview"), cex.main = 0.9)
  circos.clear()
}
pdf(file.path(PLOTS, "circos_genome.pdf"), 9, 9); circos_draw(); invisible(dev.off())
png(file.path(PLOTS, "circos_genome.png"), 9, 9, units = "in", res = 150); circos_draw(); invisible(dev.off())

# LINX plots of the overlapping clusters; if missing, the SvVisualiser command (oncoanalyser 3.0.0
# ran LINX 2.3.1, image hmftools-linx:2.3.1--hdfd78af_0, modules/local/linx/visualiser/main.nf)
linx_png <- list.files(ONCO, "\\.png$", recursive = TRUE, full.names = TRUE)
linx_png <- linx_png[grepl("linx", linx_png, ignore.case = TRUE)]
ids <- unique(as.integer(unlist(strsplit(conc$linx_clusters, ","))))
ids <- ids[!is.na(ids)]
plot_hits <- rbindlist(lapply(ids, function(i) {
  f <- linx_png[grepl(paste0("cluster-0*", i, "[^0-9]"), basename(linx_png))]
  data.table(clusterId = i, plot = if (length(f)) f else NA_character_)
}))
missing_ids <- plot_hits[is.na(plot), unique(clusterId)]
linx_cmd <- character()
if (length(missing_ids)) {
  sif <- Sys.glob(file.path(opt$container_dir, "*hmftools-linx-2.3.1*"))[1]
  if (is.na(sif)) sif <- "docker://quay.io/biocontainers/hmftools-linx:2.3.1--hdfd78af_0"
  ens <- list.files(opt$hmf_ref_dir, "^ensembl_gene_data\\.csv$", recursive = TRUE, full.names = TRUE)[1]
  ens <- if (is.na(ens)) "<ensembl_data dir: find <HMF ref dir> -name ensembl_gene_data.csv>" else dirname(ens)
  linx_cmd <- sprintf(paste0(
    "singularity exec -B /common %s bash -c 'linx -Xmx8g com.hartwig.hmftools.linx.visualiser.SvVisualiser ",
    "-sample %s -vis_file_dir %s -purple_dir %s -amber_dir %s -cobalt_dir %s -ref_genome_version 38 ",
    "-ensembl_data_dir %s -circos $(which circos) -clusterId %s -plot_out %s/ -data_out %s/'"),
    sif, TUM, dirname(F$linx_vis), dirname(F$purple_cn), file.path(ONCO, "amber"), file.path(ONCO, "cobalt"),
    ens, paste(missing_ids, collapse = ","), file.path(PLOTS, "linx"), file.path(OD, "linx_data"))
}

# ---- 7. report ----------------------------------------------------------------------------------
fmt <- function(x) if (is.null(x) || !length(x) || is.na(x)) "-" else as.character(x)
verdict <- if (!length(called)) "No chromosome meets the high- or low-confidence criteria in either set." else
  paste(vapply(called, function(cc) {
    r <- conc[chrom == cc]
    both <- r$call_A != "none" && r$call_B != "none"
    sprintf("- **chr%s**: set A %s, set B %s - %s. ESVEE-Manta: %d of %d set A SVs in the region matched at both ends (+/-%d bp, same strands). LOH low state %s / heterozygous high state %s (PURPLE). LINX cluster(s) %s (%s)%s.",
            cc, r$call_A, r$call_B, if (both) "supported by both pipelines" else "supported by one pipeline only",
            r$A_matched_both_ends, r$sv_A_in_region, TOL, fmt(r$low_state_loh_A), fmt(r$high_state_het_A),
            fmt(r$linx_clusters), fmt(r$linx_resolvedType),
            if (nzchar(r$linx_other_chromosomes)) paste0(", also on chr ", r$linx_other_chromosomes, " (possible multi-chromosome event)") else "")
  }, ""), collapse = "\n")
show <- function(x) paste(c(paste(names(x), collapse = " | "), paste(rep("---", ncol(x)), collapse = " | "),
                            apply(x, 1, function(r) paste(trimws(format(r)), collapse = " | "))), collapse = "\n")
kc <- c("chrom", "call_A", "call_B", "clusterSize_A", "clusterSize_B", "osc2_A", "osc2_B", "jaccard",
        "sv_A_in_region", "A_matched_both_ends", "A_centromere", "A_segdup", "low_state_loh_A", "high_state_het_A", "linx_clusters")
writeLines(c(
  paste0("# Chromothripsis - ", TUM, " vs ", NOR), "",
  paste0("Run ", format(Sys.time(), "%Y-%m-%d %H:%M"), ". ShatterSeek ", as.character(packageVersion("ShatterSeek")),
         ", genome hg38. Set A: ", SET_LAB[["A"]], " (primary); set B: ", SET_LAB[["B"]], " (orthogonal check)."), "",
  "## Verdict", "", verdict, "",
  "## Criteria", "",
  "Cortes-Ciriano et al. 2020 (ShatterSeek tutorial section 3; Supplementary Note p. 22): high = >= 6 interleaved intrachromosomal SVs, >= 7 adjacent 2-state CN oscillations, fragment joins test, and chromosomal enrichment or exponential breakpoint test; or >= 3 interleaved SVs + >= 4 interchromosomal SVs, >= 7 oscillations and fragment joins test. Low = as the first, with 4-6 oscillations.",
  paste0("Tests BH-corrected across this sample's chromosomes, significance q <= 0.2; fragment joins passed when ",
         if (opt$fragment_joins == "fdr") "q > 0.2" else "raw p > 0.05", "."), "",
  "## Purity / ploidy", "", show(pp), "",
  if (pp_flag) "PURPLE and ASCAT disagree (> 0.10 purity or > 0.3 ploidy)." else "PURPLE and ASCAT agree (within 0.10 purity and 0.3 ploidy).", "",
  "## Candidate chromosomes", "",
  if (nrow(conc)) show(conc[, ..kc]) else "none (no cluster of >= 6 interleaved SVs in either set)", "",
  "Full tables: `shatterseek_summary.tsv` (all chromosomes x set), `concordance.tsv`, `concordance_svs.tsv`.", "",
  "## LINX plots", "",
  if (nrow(plot_hits)) paste0("- cluster ", plot_hits$clusterId, ": ", fifelse(is.na(plot_hits$plot), "no plot found", plot_hits$plot)) else "no overlapping LINX cluster",
  if (length(linx_cmd)) c("", "Generate the missing ones:", "", "```", linx_cmd, "```") else character(), "",
  "## Log", "", paste0("- ", log_lines)
), file.path(OD, "report.md"))
note("done: ", OD)
