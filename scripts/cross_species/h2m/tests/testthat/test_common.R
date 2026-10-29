# Unit tests on small synthetic inputs for the pure functions in scripts/common.R.
# Run: Rscript tests/run_tests.R   (from scripts/cross_species/h2m)

# ---- ortholog mapping -------------------------------------------------------------------
orth <- data.table(
  gene_id_h     = c("H1", "H2", "H2", "H3", "H4"),
  symbol_h      = c("A", "B", "B", "C", "D"),
  gene_id_m     = c("M1", "M2a", "M2b", NA, "M4"),
  homology_type = c("ortholog_one2one", "ortholog_one2many", "ortholog_one2many", NA, "ortholog_many2many"))

test_that("classify_orthologs: one2one / one2many / no_ortholog / not_in_release", {
  cl <- classify_orthologs(orth, c("H1", "H9"))
  expect_equal(cl["H1", homology], "ortholog_one2one")
  expect_equal(cl["H1", gene_id_m], "M1")
  expect_equal(cl["H2", homology], "ortholog_one2many")
  expect_equal(cl["H2", n_mouse], 2L)
  expect_true(is.na(cl["H2", gene_id_m]))            # 1:many never gets a single mouse gene
  expect_equal(cl["H2", mouse_targets], "M2a;M2b")
  expect_equal(cl["H3", homology], "no_ortholog")
  expect_equal(cl["H9", homology], "not_in_release")
  expect_equal(cl["H4", homology], "ortholog_many2many")
})

test_that("map_orthologs: one row per target, no-ortholog genes kept", {
  m <- map_orthologs(c("H1", "H2", "H3", "H9", NA), orth)
  expect_equal(m[gene_id_h == "H2", .N], 2L)
  expect_equal(m[gene_id_h == "H3", homology_type], "no_ortholog")
  expect_equal(m[gene_id_h == "H9", homology_type], "not_in_release")
  expect_false(anyNA(m$gene_id_h))
})

test_that("ortholog_release_diff reports changed genes only", {
  o2 <- copy(orth)[gene_id_h == "H1", gene_id_m := "M1x"]
  o2 <- rbind(o2, data.table(gene_id_h = "H5", symbol_h = "E", gene_id_m = "M5", homology_type = "ortholog_one2one"))
  d <- ortholog_release_diff(orth, o2, 110, 115)
  expect_setequal(d$gene_id_h, c("H1", "H5"))
  expect_equal(d[gene_id_h == "H1", change], "mouse_target_changed")
  expect_equal(d[gene_id_h == "H5", change], "new_gene")
})

# ---- CN state calling -----------------------------------------------------------------------
t <- list(deep_del_abs = 0.5, loss_rel = 0.75, gain_rel = 1.25, amp_rel = 3, loh_minor_max = 0.5)

test_that("cn_state is ploidy-relative", {
  expect_equal(cn_state(c(0.2, 1.5, 2, 3.9, 9.5, NA), 2, t),
               c("DEEP_DEL", "NEUTRAL", "NEUTRAL", "GAIN", "AMP", "not_assessable"))
  # ploidy 3.15 (this tumour): 2 copies is a loss, 4 a gain
  expect_equal(cn_state(c(2, 3, 4), 3.15, t), c("LOSS", "NEUTRAL", "GAIN"))
})

test_that("gene_cn_state: loss from min, gain from max, partial, LOH", {
  s <- gene_cn_state(min_cn = c(1, 3, 2, 0.1, 3), max_cn = c(1, 10, 3, 3, 3),
                     min_minor = c(0, 1, 0.2, 0, NA), ploidy = 3, t = t)
  expect_equal(s$state, c("LOSS", "AMP", "LOSS", "DEEP_DEL", "NEUTRAL"))
  expect_equal(s$direction, c("loss", "gain", "loss", "loss", "neutral"))
  expect_equal(s$partial, c(FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(s$loh, c(TRUE, FALSE, TRUE, FALSE, FALSE))    # no LOH call on DEEP_DEL or without minor CN
})

# ---- synteny blocks and interval merging ----------------------------------------------------------
test_that("synteny_blocks break on human chromosome and rank jumps", {
  d <- data.table(gene_id_h = paste0("h", 1:6),
                  h_chrom = c("chr1", "chr1", "chr1", "chr2", "chr1", "chr1"),
                  h_start = c(100, 200, 300, 100, 5000, 400),
                  m_chrom = "chr4", m_start = c(10, 20, 30, 40, 50, 60))
  b <- synteny_blocks(d, max_rank_jump = 1)
  expect_equal(b[1], b[2]); expect_equal(b[2], b[3])
  expect_false(b[4] == b[3])                     # other human chromosome
  expect_false(b[5] == b[4])
  expect_equal(b[6], b[5])                       # human ranks h5 = 5, h6 = 4: |4 - 5| = 1, same block
  expect_false(synteny_blocks(d, max_rank_jump = 0)[6] == synteny_blocks(d, max_rank_jump = 0)[5])
  # human ranks in mouse order h1 h2 h3 h5 h6 = 1 2 4 3 5: jumps of 2 at h3 and at h6
  d2 <- copy(d)[gene_id_h == "h5", h_start := 250][h_chrom == "chr1"]
  b2 <- synteny_blocks(d2, 1)
  expect_equal(as.integer(factor(b2, levels = unique(b2))), c(1L, 1L, 2L, 2L, 3L))
})

genes <- data.table(chrom = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr2"),
                    start = c(100, 300, 500, 700, 20e6, 100), end = c(200, 400, 600, 800, 20e6 + 100, 200),
                    gene = paste0("g", 1:6), human_gene = paste0("G", 1:6), segment = c("s1", "s1", "s1", "s2", "s2", "s3"),
                    direction = c("gain", "gain", "neutral", "gain", "gain", "gain"),
                    block = c("B1", "B1", "B1", "B1", "B1", "B2"))

test_that("merge_runs merges same-direction neighbours and breaks correctly", {
  iv <- merge_runs(genes, "direction", c("gain", "loss"), max_gap = 5e6)
  expect_equal(nrow(iv), 4L)                     # g1-g2 | (neutral g3) | g4 | g5 (gap) | g6 (chr2)
  expect_equal(iv$genes[1], "g1,g2")
  expect_equal(iv$start[1], 100); expect_equal(iv$end[1], 400)
  expect_equal(iv$n_genes, c(2L, 1L, 1L, 1L))
  g2 <- copy(genes)[gene == "g3", direction := "gain"]
  iv2 <- merge_runs(g2, "direction", c("gain"), max_gap = 5e6)
  expect_equal(iv2$genes[1], "g1,g2,g3,g4")
  expect_equal(iv2$segments[1], "s1,s2")
  g3 <- copy(g2)[gene == "g4", block := "B9"]
  expect_equal(merge_runs(g3, "direction", "gain", 5e6, block_col = "block")$genes[1], "g1,g2,g3")
})

test_that("merge_overlapping, overlap_bp, total_bp", {
  a <- data.table(chrom = "chr1", start = c(0, 50, 200), end = c(100, 150, 300))
  m <- merge_overlapping(a)
  expect_equal(nrow(m), 2L); expect_equal(m$end[1], 150)
  expect_equal(nrow(merge_overlapping(a, gap = 60)), 1L)
  b <- data.table(chrom = c("chr1", "chr2"), start = c(100, 0), end = c(250, 10))
  expect_equal(overlap_bp(a, b), 100)            # 100-150 and 200-250
  expect_equal(total_bp(a), 250)
})

# ---- fusion feasibility --------------------------------------------------------------------------
test_that("sv_class from gene positions and strands", {
  expect_equal(sv_class("chr1", 100, "+", "chr1", 900, "+"), "DEL")
  expect_equal(sv_class("chr1", 900, "+", "chr1", 100, "+"), "DUP")
  expect_equal(sv_class("chr1", 900, "-", "chr1", 100, "-"), "DEL")
  expect_equal(sv_class("chr1", 100, "+", "chr1", 900, "-"), "INV")
  expect_equal(sv_class("chr1", 100, "+", "chr2", 900, "+"), "TRA")
  expect_true(is.na(sv_class(NA, 1, "+", "chr1", 1, "+")))
})

test_that("fusion_feasibility is a flag", {
  expect_equal(fusion_feasibility(0, 1, "DEL", NA), "not_possible_no_ortholog")
  expect_equal(fusion_feasibility(1, 1, "DEL", "DEL"), "feasible_same_sv_class")
  expect_equal(fusion_feasibility(1, 1, "DEL", "TRA"), "feasible_needs_TRA")
  expect_equal(fusion_feasibility(1, 1, "DEL", NA), "not_assessable")
})

# ---- SNVs ------------------------------------------------------------------------------------------
test_that("vcf_to_maf follows bioh2m's convention", {
  m <- vcf_to_maf(c(100L, 100L, 100L, 100L), c("C", "A", "ACG", "AT"), c("T", "AGG", "A", "AC"))
  expect_equal(m$type_h, c("SNP", "INS", "DEL", "DNP"))
  expect_equal(m$start_h, c(100L, 100L, 101L, 100L))
  expect_equal(m$end_h, c(100L, 101L, 102L, 101L))
  expect_equal(m$ref_seq_h, c("C", "-", "CG", "AT"))
  expect_equal(m$alt_seq_h, c("T", "GG", "-", "AC"))
})

test_that("short_protein / same_protein_change", {
  expect_equal(short_protein(c("p.Arg172His", "p.R172H", "(p.Gly12Asp)", NA)), c("R172H", "R172H", "G12D", NA))
  expect_equal(same_protein_change(c("p.Arg172His", NA, "p.V637E"), c("R172H", "R172H", "V600E")), c(TRUE, FALSE, FALSE))
})

test_that("pave_fields reads the header order", {
  h <- '##INFO=<ID=IMPACT,Number=10,Type=String,Description="Variant Impact [Gene, Transcript, CanonicalEffect, CanonicalCodingEffect, SpliceRegion, HgvsCodingImpact, HgvsProteinImpact, OtherReportableEffects, WorstCodingEffect, GenesAffected]">'
  expect_equal(pave_fields(h), PAVE_DEFAULT)
  expect_equal(pave_fields("no brackets"), PAVE_DEFAULT)
})

# ---- spectrum rescaling --------------------------------------------------------------------------
test_that("SBS96 / TRI32 definitions and labels", {
  expect_length(SBS96, 96); expect_length(unique(SBS96), 96); expect_length(TRI32, 32)
  expect_equal(sbs96_label(c("ACG", "CGT"), c("C", "G"), c("T", "A")), c("A[C>T]G", "A[C>T]G"))
})

test_that("collapse_trinuc folds purine-centred contexts", {
  k64 <- as.vector(outer(outer(BASES, BASES, paste0), BASES, paste0))
  c64 <- setNames(rep(1, 64), k64)
  t32 <- collapse_trinuc(c64)
  expect_equal(sum(t32), 64); expect_true(all(t32 == 2))
})

test_that("rescale_spectrum: identity with equal genomes, reweights otherwise, keeps total", {
  tri <- setNames(rep(1, 32), TRI32)
  cnt <- setNames(rep(0, 96), SBS96); cnt["A[C>T]G"] <- 10; cnt["T[T>C]T"] <- 10
  s <- rescale_spectrum(cnt, tri, tri)
  expect_equal(s$expected_mouse_count, s$human_count)
  tri_m <- tri; tri_m["ACG"] <- 0.5                         # CpG half as frequent in "mouse"
  s2 <- rescale_spectrum(cnt, tri, tri_m)
  expect_equal(sum(s2$expected_mouse_count), 20)
  expect_lt(s2[MutationType == "A[C>T]G", expected_mouse_count], 10)
  expect_gt(s2[MutationType == "T[T>C]T", expected_mouse_count], 10)
})

# ---- direction / role -------------------------------------------------------------------------
test_that("event_direction, role_consistency, combine_directions", {
  expect_equal(event_direction("LOSS", NA), "inactivating")
  expect_equal(event_direction("AMP", NA), "activating")
  expect_equal(event_direction("MISSENSE", "oncogene"), "activating")
  expect_equal(event_direction("MISSENSE", "TSG"), "inactivating")
  expect_equal(event_direction("MISSENSE", "oncogene, TSG"), "unknown")
  expect_equal(role_consistency("inactivating", "TSG"), "consistent")
  expect_equal(role_consistency("inactivating", "oncogene"), "discordant")
  expect_equal(role_consistency("activating", NA), "")
  expect_equal(combine_directions(c("activating", "unknown")), "activating")
  expect_equal(combine_directions(c("activating", "inactivating")), "mixed")
})

test_that("config substitution and chain names", {
  expect_equal(subst("{a}/x/{b}", list(a = "{b}/y", b = "z")), "z/y/x/z")
  expect_equal(ucsc_chain_name("hg38", "mm10"), "hg38ToMm10.over.chain.gz")
  expect_equal(norm_chrom(c("1", "chrX", "MT")), c("chr1", "chrX", "chrM"))
})

test_that("orthologs_from_compara merges the human- and mouse-folder files", {
  hdr <- "gene_stable_id\tprotein_stable_id\tspecies\tidentity\thomology_type\thomology_gene_stable_id\thomology_protein_stable_id\thomology_species\thomology_identity\tdn\tds\tgoc_score\twga_coverage\tis_high_confidence\thomology_id"
  fh <- tempfile(fileext = ".tsv.gz"); fm <- tempfile(fileext = ".tsv.gz")
  writeLines(c(hdr, "H1\tP\thomo_sapiens\t90\tortholog_one2one\tM1\tQ\tmus_musculus\t80\tNULL\tNULL\tNULL\tNULL\t1\t1",
               "H1\tP\thomo_sapiens\t50\tortholog_one2one\tX1\tQ\tmus_spretus\t50\tNULL\tNULL\tNULL\tNULL\t1\t2",
               "H3\tP\thomo_sapiens\t40\tother_paralog\tH4\tQ\thomo_sapiens\t40\tNULL\tNULL\tNULL\tNULL\tNULL\t3"), gzfile(fh))
  writeLines(c(hdr, "M2\tQ\tmus_musculus\t70\tortholog_one2one\tH2\tP\thomo_sapiens\t60\tNULL\tNULL\tNULL\tNULL\t0\t4"), gzfile(fm))
  hg <- data.table(gene_id = c("H1", "H2", "H3"), gene_name = c("A", "B", "C"), chrom = "chr1", start = 1L, end = 2L)
  mg <- data.table(gene_id = c("M1", "M2"), gene_name = c("a", "b"), chrom = "chr4", start = 1L, end = 2L)
  o <- orthologs_from_compara(c(homo_sapiens = fh, mus_musculus = fm), hg, mg, 110)
  expect_equal(o[gene_id_h == "H1", gene_id_m], "M1")                  # other Mus species ignored
  expect_equal(o[gene_id_h == "H2", .(gene_id_m, perc_id_h_to_m, perc_id_m_to_h)],
               data.table(gene_id_m = "M2", perc_id_h_to_m = 60, perc_id_m_to_h = 70))   # mouse-side row swapped
  expect_true(is.na(o[gene_id_h == "H3", gene_id_m]))                  # paralogs ignored, gene kept
  expect_equal(o[gene_id_h == "H2", symbol_m], "b")
})

test_that("parse_hgvsp, map_aligned_pos, cds_protein", {
  p <- parse_hgvsp(c("p.Arg250Trp", "p.(Gly12Asp)", "p.Lys120fs", "p.?", NA, "p.R175H"))
  expect_equal(p$aa, c("R", "G", "K", NA, NA, "R")); expect_equal(p$pos, c(250L, 12L, 120L, NA, NA, 175L))
  # human MKRLE vs mouse MK-LE: R (pos 3) faces a gap; L (pos 4) is mouse pos 3
  expect_equal(map_aligned_pos("MKRLE", "MK-LE", 3)$mouse_aa, "-")
  expect_equal(map_aligned_pos("MKRLE", "MK-LE", 4)[c("mouse_pos", "mouse_aa")], list(mouse_pos = 3L, mouse_aa = "L"))
  expect_equal(map_aligned_pos("MKRLE", "MKKLE", 3)$mouse_aa, "K")
  # + strand: ATG AAA | CGT TAA -> MKR*; - strand pieces given as + strand sequence
  expect_equal(cds_protein(data.table(start = c(1, 20), strand = "+", frame = "0"), c("ATGAAA", "CGTTAA")), "MKR*")
  rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
  expect_equal(cds_protein(data.table(start = c(1, 20), strand = "-", frame = "0"), c(rc("CGTTAA"), rc("ATGAAA"))), "MKR*")
})
