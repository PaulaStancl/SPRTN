# Mouse WGS (Sprtn knock-in HCC): status, findings and how to continue

Last updated 2026-11-25 (section 6: everything agreed and done on 2026-11-25). Plan: Claude Docs "Mouse WGS → human comparison: analysis plan"
(https://claude.ai/artifact/Wm6n1M8AwjP9VZsCRqocbK). This README records what has been done since.

**Most important finding so far: the setting-A tumour/normal pairs for mice 1 and 2 are swapped.**
HCC_1 belongs with KI2 and HCC_2 with KI1 (section 3). Until setting A is re-called with the correct
pairs, all setting-A results for HCC_1 and HCC_2 are unreliable as somatic calls.

---

## 1. Data

Everything is on ssi under `/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/data/mouse/wgs/`.
There are no BAMs; all analyses use the VCF/BCF files. Reference GRCm38, contigs `1..19, X, Y`.

| Setting | Folder | Content |
|---|---|---|
| A: tumour vs matched KI normal | `tumor_vs_normal/snv_indel/` | Strelka2, Mutect2 and FreeBayes per tumour, plus `merged/*.union.ann.bcf` (Patrik's merge, INFO/CALLER) |
| B: tumour vs joint WT normals | `tumor_vs_joint_wt/` | Mutect2 joint calls; the only files with the 3 WT normals |
| C: tumour-only + WT PoN | `tumor_only_pon/` | Mutect2 |
| SVs | `tumor_vs_normal/sv/{calls,merged}/` | Delly, SvABA, Manta; OctopuSV merge |

Samples (short names used in all scripts):

| Short | Full sample name | Labelled pairing |
|---|---|---|
| HCC_1 | id_37264_sfb_dikic_2025_05_10_KI_HCC_1 | KI1 |
| HCC_2 | id_13708_sfb_dikic_2025_05_11_KI_HCC_2 | KI2 |
| HCC_3 | id_11440_sfb_dikic_2025_05_12_KI_HCC_3 | KI3 |
| KI1 | id_37264_sfb_dikic_2025_05_05_KI1 | |
| KI2 | id_13708_sfb_dikic_2025_05_06_KI2 | |
| KI3 | id_11440_sfb_dikic_2025_05_07_KI3 | |

The `id_` prefixes group HCC_1 with KI1, HCC_2 with KI2 and HCC_3 with KI3, matching the labels, not the
genotypes. If these are per-mouse submission IDs, the swap happened before sequencing (tissue or labels).

PASS rules applied by the provider (Patrik):
- **Strelka2:** sarek defaults.
- **Mutect2:** FilterMutectCalls with sarek's orientation priors.
- **FreeBayes:** crude prefilter `INFO/MQM >= 40 && GT[normal]="0/0" && GT[tumor]="0/1" && FMT/DP[normal] >= 30 && FMT/DP[tumor] >= 30 && QUAL/INFO/AO >= 5 && SAF > 0 && SAR > 0`.
  FreeBayes writes PASS when the expression passes and `.` when it fails.

Files were split with `bcftools norm -m -any` (no `-f`, so indels are **not left-aligned**). The annotation
is **VEP CSQ**, not SnpEff ANN.

---

## 2. Scripts

All are tutorial-style, numbered blocks, run line by line in RStudio or with `Rscript` / `bash` on ssi.

- **Env:** `r-mutation` (`/common/WORK/pstancl/envs/r-mutation`). It contains bcftools, but only version
  **1.4.1** (2017).
- **Outputs:** `results/mouse_wgs/<script>/`, each with a `RUN_LOG.md`.

| Script | What it does | Status |
|---|---|---|
| `00_input_inspection.R` | Task 0: provenance and md5s, headers, normals used per setting, record counts, SV inventory, manifest, issues | written, result not reviewed |
| `mouse_wgs_snv_indel_callers.R` | Callers within setting A. Blocks 3–8: counts, pairwise sharing, 1/2/3-caller categories, merged union, SBS96 (GRCm38 FASTA from h2m resources). R1–R7: raw single-caller PASS → intersect vs merged file, merged FILTER rule check (R5c), Venn and bar plot of my intersect vs merged | run; R1–R7 reviewed; revised 2026-11-25 (section 6.1), re-run needed |
| `mouse_wgs_snv_indel_settings.R` | Settings compared per tumour, 9 sets (PASS in each file's own FILTER): A_freebayes (file PASS), A_freebayes_pass_expr (PASS & prefilter expression), A_freebayes_expr (expression only, FILTER not used), A_mutect2_strelka (PASS then intersect), A_ge2_callers (PASS in >=2 of FreeBayes PASS&expr, Mutect2, Strelka), A_strelka, A_mutect2, B_mutect2_joint_wt, C_mutect2_pon. Counts (all vs PASS), pairwise shared, exclusive categories, proportion shared, SBS96 (ggplot + qcVCF), cosine between settings and between tumours | 6-set version run 2026-11-25; 8-set version (FreeBayes x3) not run yet; see 2.2 |
| `mouse_wgs_freebayes_pass_check.R` | Checks that every FreeBayes PASS record (setting A, own files) satisfies the provider's prefilter (MQM, GT normal 0/0 and tumour 0/1, DP >= 30 both, QUAL/AO >= 5, SAF/SAR > 0), and the reverse; failing conditions per record | run 2026-11-25 (result 3.6); split-piece table not reviewed yet |
| `driver_genes_human/mouse_wgs_driver_genes_human.R` | Harmful (VEP IMPACT HIGH/MODERATE) PASS mutations (SNVs, MNVs split, and indels; any of the 9 settings sets) in mouse orthologs (Ensembl Compara, h2m table) of human cancer genes; separate analyses: IntOGen HCC (Compendium, TUMOUR_TYPE == HCC, 106 genes), IntOGen all tumour types, COSMIC CGC (optional, file not on ssi yet); patient_driver = gene is_driver in the RJALS tumourevo driver.rds (AXIN1, KMT2A); one TRUE/FALSE column per set; gene x set heatmaps. Needs settings.rds and the IntOGen + CGC file paths | new 2026-11-25, not run yet |
| `driver_genes_human/mouse_wgs_driver_lollipop.R` | trackViewer lollipop plots of the patient's driver genes (RJALS tumourevo is_driver: AXIN1, KMT2A): human protein with the patient's coding mutations + mouse ortholog with PASS coding mutations (any impact), one track per settings set; Pfam domains and protein length from Ensembl REST; table of coding variants x 9 sets. Runs in its own env r-lollipop (trackViewer cannot be solved into r-mutation; create command in the script header) | new 2026-11-25, not run yet |
| `mouse_wgs_example_sites.R` | Diagnostics for chr19:30419577 T>C and chr10:3184122 G>T from the single-caller files: per-caller metrics, FreeBayes criteria, per-caller read counts in all samples, plots 01–08 | run, reviewed |
| `mouse_wgs_sample_identity.R` | Sample identity from FreeBayes genotypes: IBS0 rates per shared confident sites (WT from Mutect2 joint), LD-thinned rate, IBS0 blocks and chromosomes, depth and GQ/MQM at IBS0 sites, optional MGP strains, bcftools gtcheck | v1 run (counts in section 3); **current version not run yet** |
| `gtcheck_freebayes.sh` | bcftools gtcheck: biallelic autosomal SNVs, QUAL/DP/GQ filter, merge `-m none`, all-sites and common-sites runs, per-1000 table | run; common-sites rows not reviewed yet |
| `mouse_wgs_somalier.R` | Identity: normal-homozygous discordance (primary), IBS0, somalier and gtcheck on one masked 6-sample genotype file; panels A/B/C × thinning, bootstrap, jackknife; see 2.1 | revised 2026-11-25 (f5e1226); `--test` (swap + null) passes on ssi incl. somalier n check and gtcheck 1.24; real data run 2026-11-25: swap confirmed, HCC_3-KI3 confirmed (section 3.3) |
| `mouse_wgs_identity_informativeness.R` | Reads the outputs of `mouse_wgs_somalier.R` (masked 6-sample BCF, primary sites): genotype patterns across the KI normals, informative sites, opposite homozygotes per pair of normals, private sites, tumour-dependent retention, identity on a tumour-independent site set (all normals non-ref), sites per chromosome, and `bam_check_sites_for_patrik.csv` | run on real data 2026-11-25 (section 6.3) |
| `mouse_wgs_identity_strain_check.R` | Are the identity sites and the tumour vs labelled-normal differences germline? Panel A sites annotated with the Mouse Genomes Project (GRCm38 release: rsID/ID, strains carrying the ALT, C57BL/6NJ and 129 strains), WT normals from the setting-B Mutect2 files, and a neighbouring-SNP check (inherited blocks vs isolated sites, ±5 Mb); annotates `bam_check_sites_for_patrik.csv` | run 2026-11-25 without MGP (0 records over ftp); re-run with the local MGP copy pending (section 6.4) |
| `mouse_wgs_identity_discordant_sites.R` | Lists the germline SNVs where a tumour differs from a KI normal, from all sites of the masked BCF (red triangles of the slide "Allele fractions at all germline SNPs" and the ticks of "Where each tumour differs from its labelled normal", ~520 per labelled pair; genome_slide column = blue inherited / red unexplained): genotypes, read counts and VAF of all six samples, type (opposite homozygotes / het vs hom), whether in Patrik's BAM list; counts_3x3.tsv to match the slide | written 2026-11-26, not run |
| `for_patrick/for_patrick.R` | Three CSVs for Patrik: all identity sites wide (GT, REF/ALT reads per sample, discordance vs labelled normal); 20 discordant (opposite homozygotes in both labelled pairs) + 10 concordant (all 0/1) sites for IGV; IntOGen driver genes with mouse orthologs. Run line by line | written 2026-11-27, not run |
| `mouse_wgs_snv_indel.R` | Task 1: settings A/B/C (counts, VAF, FILTER, overlaps, annotation) | written, not run; block 9 reads ANN and must switch to **CSQ** |

Run commands (from the repo on ssi):
```
cd /common/WORK/pstancl/projects/MariaBoskovic/SPRTN
git pull
Rscript scripts/mouse_wgs/<script>.R          # or: bash scripts/mouse_wgs/gtcheck_freebayes.sh
```

### 2.1 `mouse_wgs_somalier.R` (rewritten 2026-11-24)

**Revised 2026-11-25 after review (CHANGE 22–39 in the script).**
- **Why:** the ascertainment means IBS0 can detect a swap but can't confirm correct labels. Each FreeBayes
  file only holds sites where its pair is non-ref, so a correctly paired animal is never 0/0 at the
  retained sites.
- **New primary metric:** discordance at sites where the **normal is homozygous** (Conpair-style).
- **New statistics:**
  - paired sign tests, using sites where at least one candidate normal is homozygous; the spec's
    both-homozygous version is in the `*_bothhom` columns;
  - 10 Mb block bootstrap, leave-one-chromosome-out jackknife and discordant-segment counts;
  - max-depth filter and optional `MASK_BED`.
- **FreeBayes PASS sites are kept by default** (`EXCLUDE_FB_PASS = FALSE`); excluding them is a
  sensitivity run.
- **Other sensitivity runs:** the files-1+2 ascertainment and MIN_DP 15/20.
- **Test:** `--test` now runs a swap scenario and a null scenario (all labels correct). Both pass on the
  Mac (R path only).
- **On ssi:** run `--test` first (somalier n check, gtcheck with bcftools 1.24), then the real data. The
  design notes below describe the 2026-11-24 version; the script header has the full changelog.

Why it was rewritten:
- The first real run crashed at block 5: somalier v0.3.5 `pairs.tsv` has no `hom_concordance` column. The
  site list had built fine (814,298 candidates, 62,946 germline-like, 12,339 after 20 kb bins).
- Paula's review spec also asked for stricter site selection, masking, panels, three methods and full logging.
  The spec is in the conversation of 2026-11-24. Each code change is tagged `[CHANGE n]`, and the script
  header has a CHANGELOG (1–21).

Design:
- **Sites:**
  - biallelic SNVs on autosomes 1–19, QUAL ≥ 30, from the 3 FreeBayes BCFs;
  - REF+ALT ≥ 30 in **all six** samples. This is the ascertainment fix: each BCF only holds sites variant in
    its own pair, so the old ≥ 4-of-6 rule compared labelled pairs at more sites than cross pairs;
  - genotype bins on VAF: ≤ 0.05 hom-ref, 0.35–0.65 het, ≥ 0.95 hom-alt, anything else ambiguous;
  - every KI normal confidently genotyped, at least one of them non-ref;
  - the site filters use only the KI normals. **Correction 2026-11-25:** the site list still depends on the
    tumours. A FreeBayes file only holds a site if its tumour or its normal is non-ref, so where a labelled
    normal is 0/0 the site is in that file only because the tumour is non-ref ("tumour-dependent
    retention"). It cannot create the swap result: the zeros for the genetic matches and the 15-23%
    discordance of the labelled pairs are genotypes at the sites, not which sites were kept. Checked on a
    tumour-independent site set (all three normals non-ref) in `mouse_wgs_identity_informativeness.R`;
  - tumour genotypes that are ambiguous are masked (`./.`), never set to hom-ref.
- **One masked file:** `masked/masked_6samples.bcf` (GT:AD, INFO/AF from the KI normals, no GL/PL).
  somalier, gtcheck and the custom IBS0 all read this file.
- **Panels:**
  - A: all sites;
  - B: KI1–3 not all identical;
  - C: opposite homozygotes in ≥ 2 KI normals.

  Each panel is thinned by greedy minimum distance, keeping the site with the best minimum depth: 20 kb
  (primary), 100 kb, 500 kb and 1 Mb. A panel × thinning with fewer than 50 sites is marked underpowered.
- **Methods:**
  - **custom IBS0** (primary): IBS0 / jointly callable sites, with exact 95% CIs;
  - **somalier** `relate -s sites -d 30 --min-ab 0.35`, with `SOMALIER_AB_HOM_CUTOFF=0.05`;
  - **bcftools gtcheck** `-G 1 -T targets`.
- **Matching:**
  - per tumour, the lowest IBS0 rate wins and IBS2 breaks ties;
  - a call is ambiguous if the best and second-best CIs overlap, the joint sites are < 50, or the panel is
    underpowered;
  - all 3! one-to-one assignments are scored;
  - the labelled pairing and the hypothesis HCC_1–KI2, HCC_2–KI1, HCC_3–KI3 are both tested.
- **Outputs:**
  - tables: `verdict.tsv`, `robustness.tsv`, `assignment_scores.tsv`, `pairwise_metrics.tsv`,
    `filter_counts.tsv`, `panel_sizes.tsv`;
  - plots: heatmaps for 6 metrics, `ibs0_per_tumour`, `ibs0_vs_ibs2`, `robustness`,
    `chromosomal_distribution`;
  - `RUN_LOG.md`, which holds the changelog, assumptions, limitations and DEVIATION lines.

Tool facts, checked against the source code:
- **somalier v0.3.5:**
  - The binary is at `/common/WORK/pstancl/tools/somalier/somalier`, installed by the osteosarcoma
    `00_somalier.sh`. It is not in a micromamba env.
  - The version is written to stderr on every call (`--version` is not an option).
  - `extract` needs an indexed VCF/BCF.
  - Genotypes are called from AD: hom if AB < `SOMALIER_AB_HOM_CUTOFF` (default 0.01) or > 1 minus it; het
    if `min-ab` ≤ AB ≤ 1 − `min-ab`; unknown otherwise or below `-d`.
  - relatedness = 2·(shared_hets − 2·ibs0) / het_ab (KING-robust style, no allele frequencies).
  - `concordance` is somalier's adjusted homozygous concordance.
  - n = sites where both samples have a known genotype.
- **bcftools 1.4.1 gtcheck** (no `-g`): cross-checks all samples and writes
  `ERR  error_rate  n_sites  sample_i  sample_j`. With `-G 1` it compares GTs only.

Run on ssi, the synthetic test first:
```
cd /common/WORK/pstancl/projects/MariaBoskovic/SPRTN
git pull origin main
Rscript scripts/mouse_wgs/mouse_wgs_somalier.R --test 2>&1 | tee somalier_test.log
grep -E "TEST_MODE|DEVIATION|Error" somalier_test.log
Rscript scripts/mouse_wgs/mouse_wgs_somalier.R 2>&1 | tee mouse_wgs_somalier.log
grep -E "\[filter\]|DEVIATION|Error|Done" mouse_wgs_somalier.log
column -t results/mouse_wgs/mouse_wgs_somalier/verdict.tsv
```
- **Test mode:** writes 3 synthetic pairs (HCC_1/HCC_2 swapped) to
  `results/mouse_wgs/mouse_wgs_somalier_TEST/`. It must end with "TEST_MODE: all checks passed", including the
  six-`.somalier`-files check, which is skipped on the Mac.
- **Things to watch in the real run:**
  - how many sites survive "all six samples" (expect far fewer than 62,946);
  - whether panels B/C at 500 kb or 1 Mb are underpowered;
  - whether somalier and gtcheck agree with the custom IBS0.
- **Open choices, defaults in block 0:**
  - `MIN_JOINT_SITES = 50`;
  - an assignment is "clearly separated" when its gap to the next is > 2·√(summed binomial variances);
  - `KNOCKIN_LOCUS = NA`. Set the Sprtn knock-in position to annotate IBS0 near it.

---

### 2.2 Settings comparison and FreeBayes PASS check (2026-11-25)

**`mouse_wgs_snv_indel_settings.R`** compares, per tumour, 9 sets. PASS is always each file's own FILTER;
the merged file is not used. MNVs split into SNVs; doublet halves left out of SBS96.

| Set | Definition |
|---|---|
| `A_freebayes` | tumor_vs_normal, FreeBayes FILTER PASS as in the file |
| `A_freebayes_pass_expr` | FreeBayes PASS **and** the prefilter expression (3.6) on the split record |
| `A_freebayes_expr` | FreeBayes, the prefilter expression only; the FILTER column is not used (per variant, MNVs split) |
| `A_mutect2_strelka` | PASS first in Mutect2's and in Strelka's own file, then intersected ("all records" = found by both) |
| `A_ge2_callers` | PASS in >= 2 of: `A_freebayes_pass_expr` (FreeBayes PASS AND expression), `A_mutect2` (PASS), `A_strelka` (PASS); "all records" = found in >= 2 of the 3 files |
| `A_strelka` | Strelka PASS |
| `A_mutect2` | Mutect2 PASS |
| `B_mutect2_joint_wt` | `tumor_vs_joint_wt/`, Mutect2 PASS |
| `C_mutect2_pon` | `tumor_only_pon/`, Mutect2 PASS |

Outputs (`results/mouse_wgs/mouse_wgs_snv_indel_settings/`): counts all vs PASS (`counts_per_setting*`,
qcVCF `plot_mutation_counts`), pairwise shared heatmap (`pairwise_shared`, labels = % of row, 1 decimal, as
qcVCF), qcVCF pairwise heatmaps over all tumour × set (`qcVCF_pairwise_shared_SNV/INDEL`), exclusive
categories, proportion shared, SBS96 per tumour (ggplot + qcVCF `plot96_matrix` all / PASS / not PASS),
cosine between sets and between tumours, `sbs96_matrix.tsv`, `pass_variants_per_setting.tsv.gz`.
First run 2026-11-25 (6 sets, 4.2 min); the 8-set version (commit 443d3cc) is not run yet.

**qcVCF pairwise fix.** `qcVCF::plot_pairwise_shared_mutations()` merged the table with itself without
`allow.cartesian`. data.table refuses a join with more rows than `nrow(x) + nrow(i)`; a variant in k columns
gives k² rows. SNVs failed (428,433 rows > 316,290 = 2 × 158,145), INDELs passed only because most
FreeBayes indels are in one column (171,222 < 212,616). Not a property of the data: fixed in qcVCF commit
7cb2088 (`allow.cartesian = TRUE`). On ssi: `cd /common/WORK/pstancl/tools/qcVCF && git pull`, then
`remotes::install_local("qcVCF", force = TRUE)` in r-mutation. The settings script keeps a fallback
(`pairwise_shared_sparse()`, same numbers and look) if the installed qcVCF is the old one.

**`mouse_wgs_freebayes_pass_check.R`** evaluates the FreeBayes prefilter expression on every record of the
setting-A FreeBayes files and compares it with FILTER in both directions; writes the agreement table,
failing conditions, the disagreeing records and the split-piece table. Result in 3.6.

**Union check (chat code, not a script).** R code given in chat to read one tumour's `*.union.ann.bcf`, set the
tumour/normal sample columns by hand, evaluate the expression and cross-tabulate it with the merged FILTER
for records whose INFO/CALLER includes FreeBayes. Caveat: the merged QUAL is the bcftools-merge QUAL (not
necessarily FreeBayes'), so QUAL/AO can differ from the FreeBayes file. Result not reviewed yet.

## 3. Findings

### 3.1 The merged (union) file's FILTER is the bcftools-merge union, not a consensus

The rule was checked on all 6.8 M merged records with no exception (R5c): 6,707,589 FALSE/FALSE and
112,549 TRUE/TRUE.
- A merged record is PASS if at least 1 caller is PASS and no caller that found the variant has a named
  filter. FreeBayes `.` never blocks PASS.
- INFO/CALLER lists every caller that **called** the variant, passed or not.
- So "merged PASS + ≥ 2 callers in INFO/CALLER" is **not** "≥ 2 callers PASS".

Evidence:
- **Raw route vs merged route,** PASS ≥ 2 callers, chr19 test run: my intersect (each caller's own PASS,
  then overlap) gave HCC_1 173, HCC_2 220, HCC_3 86. The merged file gave 175, 272, 84.
- **Merged non-PASS despite 2 callers PASS in their own files** (genome-wide):
  - Strelka2 PASS + Mutect2 filtered + FreeBayes PASS: 2,241 / 537 / 230.
  - Strelka2 filtered + Mutect2 PASS + FreeBayes PASS: 58 / 25 / 6.
- **chr19:30419577 T>C, HCC_1:** Mutect2 PASS, FreeBayes PASS, Strelka2 LowEVS → merged LowEVS.
- **chr10:3184122 G>T, HCC_1:** Strelka2 PASS, FreeBayes `.` → merged PASS with CALLER "strelka,freebayes".

### 3.2 The two example sites are not good somatic calls

- **chr19:30419577 T>C:**
  - VAF is similar in tumour and matched normal: HCC_1 0.10–0.25 vs KI1 0.09–0.12; HCC_3 0.19 vs KI3 0.10.
  - Mutect2 PASS is borderline: TLOD 6.54 (cutoff 6.3), only 20/15 reads, ALT base quality 20, NLOD 4.07.
  - Strelka2 SomaticEVS 0.1.
  - Likely an artefact or low-level background. The "≥ 2 callers PASS" rule keeps it, so Task 2 (alt
    support in the matched normal) must remove it.
- **chr10:3184122 G>T:**
  - The VAFs look germline: HCC_1 1.00, KI1 0.00, HCC_2 0.00, KI2 0.99, HCC_3 0.65, KI3 0.43.
  - As labelled this is impossible within one animal, but it is consistent once KI1 and KI2 are swapped.
    This is the site that led to section 3.3.

### 3.3 Sample swap: HCC_1 = KI2 animal, HCC_2 = KI1 animal (putative until the records confirm it)

**High-confidence opposite homozygotes (IBS0):** DP ≥ 30 in both samples, ref VAF ≤ 0.05, alt VAF ≥ 0.95,
FreeBayes genotypes.

| Pair | Concordance | IBS0 |
|---|---|---|
| HCC_1 – KI2 | 0.988 | 60 |
| HCC_2 – KI1 | 0.982 | 73 |
| HCC_3 – KI3 (labelled) | 0.951 | 190 |
| HCC_1 – KI1 (labelled) | 0.734 | 1,638 |
| HCC_2 – KI2 (labelled) | 0.737 | 1,402 |
| KI1 – KI2 | 0.753 | 1,182 |
| HCC_1 – HCC_2 | 0.756 | 1,242 |

The HCC_3–KI3 IBS0 sites are spread over all autosomes, roughly proportional to chromosome size.

**bcftools gtcheck** (`gtcheck_freebayes.sh`, bcftools 1.4.1; sites kept when both samples of a pair have
DP ≥ 30 and GQ ≥ 30). Discordance per 1000 sites, all-sites run:

| Pair | per 1000 |
|---|---|
| HCC_2 – KI1 | 40.2 |
| HCC_1 – KI2 | 43.0 |
| KI3 – KI2 | 46.2 |
| HCC_3 – HCC_1 | 50.1 |
| HCC_3 – KI3 (labelled) | 67.7 |
| HCC_1 – KI1 (labelled) | 118.8 |
| HCC_2 – KI2 (labelled) | 130.4 |
| KI1 – KI2 | 140.9 |
| HCC_1 – HCC_2 | 145.2 |

- In the all-sites run, labelled pairs share one file and are compared at about twice as many sites
  (220–274k vs ~100–127k), which inflates their discordance.
- The **common-sites** run puts every pair on the same 82,949 sites. Its rows are in `gtcheck_table.tsv` but
  have not been reviewed yet. They decide mouse 3.
- The first gtcheck run (QUAL ≥ 30 only) gave the same order with weaker separation: 85.6 / 91.3 vs
  125.1 / 134.5.

**Interpretation:**
- Mice 1 and 2 are each consistent with one animal (A = HCC_1 + KI2, B = HCC_2 + KI1). The labelled pairs
  are as different as the two normals or the two tumours are from each other.
- Mouse 3 is not resolved: it is genetically close to animal A, and HCC_3–KI3 is not as clean as the A and B
  pairs.
- Terminology: VAF gives an inferred genotype (diploid assumption), not copy number. The ~60–73 IBS0 in true
  pairs is an unexplained baseline (misclassification, mapping, somatic + LOH, allele-specific coverage).
- These checks use the same FreeBayes calls, so they are two methods but not independent evidence. Lab
  records, or a pileup if BAMs ever become available, would be independent.

**somalier / custom IBS0 rewrite** (`mouse_wgs_somalier.R`, real data, 2026-11-25). The run used one
masked 6-sample file: sites with REF+ALT ≥ 30 in all six samples; KI normals confidently genotyped;
ambiguous tumour genotypes masked. The table is the primary result: panel A, thinned at 20 kb, about 7,200
jointly callable sites per tumour–normal pair.

| Tumour | Labelled normal: IBS0 rate | Best normal | Best IBS0 (95% CI upper) | Second best: IBS0 (CI lower) | Ambiguous | somalier best |
|---|---|---|---|---|---|---|
| HCC_1 | KI1: 0.0714 (n = 7,228) | **KI2** | 0 (0.0005) | KI3: 0.0167 (0.0139) | no | KI2 |
| HCC_2 | KI2: 0.0710 (n = 7,254) | **KI1** | 0 (0.0005) | KI3: 0.0306 (0.0268) | no | KI1 |
| HCC_3 | KI3: 0 (n = 7,190) | **KI3** | 0 (0.0005) | KI2: 0.0168 (0.0140) | no | KI3 |

- **The swap is confirmed with zero opposite homozygotes.** HCC_1 matches KI2 and HCC_2 matches KI1, each at 0 IBS0 over about 7,200 sites. The labelled
  pairs have about 7% IBS0.
- **Mouse 3 is now resolved: HCC_3 and KI3 are the same animal.** This replaces the open point above. They
  have 0 IBS0, and the next-best normal is clearly separated.
- **somalier agrees on all three tumours.**
- **gtcheck was not run.**
  - r-mutation now has bcftools 1.24. Its gtcheck has no `-G` option (it uses `-u GT`) and prints `DCv2`
    lines, whereas the script handles only the 1.4.1 syntax. The script logs this as a DEVIATION and skips
    gtcheck.
  - Decided 2026-11-25 not to add 1.24 support. gtcheck would run on the same masked genotypes, so it adds
    no independent evidence. The custom IBS0 and somalier already agree, and the result is clear-cut.
- **Still to review:**
  - `robustness.tsv` (panels B/C and the coarser thinning);
  - `assignment_scores.tsv`;
  - the plots in `results/mouse_wgs/mouse_wgs_somalier/`.
- All three methods use the same FreeBayes calls. Lab records remain the independent confirmation.

**Revised identity analysis** (`mouse_wgs_somalier.R` f5e1226, real data, 2026-11-25). This is the
final result.
- **Primary metric:** discordance at normal-homozygous sites.
- **Sites:** panel A, thinned at 20 kb (2,465 sites).
- **Defaults:** PASS sites kept, max-depth filter 2.5 × median.

| Tumour | Labelled normal: discordance | Genetic match: discordance | Next best | Paired p | Evidence vs the labelled / next normal |
|---|---|---|---|---|---|
| HCC_1 | KI1: 22.7% | **KI2: 0 / 1,481** | KI3: 8.3% | 3e-36 | 194 segments on 18 chromosomes (vs KI1) |
| HCC_2 | KI2: 15.3% | **KI1: 0 / 1,630** | KI3: 9.3% | 9e-41 | 129 segments on 19 chromosomes (vs KI2) |
| HCC_3 | KI3: 0 / 1,425 | **KI3** | KI2: 11.0% | 7e-49 | 101 segments on 18 chromosomes (vs KI2) |

- **No tumour is ambiguous.** IBS0, somalier and gtcheck (bcftools 1.24) all give the same best normal.
- **HCC_3–KI3 is now positively confirmed.** IBS0 could not confirm a correct label because of the per-pair
  ascertainment; the normal-homozygous discordance can.
- **Robustness:**
  - block bootstrap (255 × 10 Mb blocks): D = 0.38 (95% CI 0.33–0.44), D > 0 in 100% of replicates, and
    support 1.0 for the assignment HCC_1–KI2, HCC_2–KI1, HCC_3–KI3;
  - jackknife: no chromosome changes the verdict or an ambiguity flag (weakest paired p ≤ 6e-27);
  - files-1+2 ascertainment: p = 1e-308 and 1e-227;
  - MIN_DP 15 and 20: same verdict;
  - PASS exclusion: same verdict (it removes only 118 of 3,188 sites).
- **The evidence is genome-wide** (18–19 of the 19 autosomes), not one inherited block.
- **Open:** the max-depth filter removes many sites (panel A dropped from about 7,200 to 3,188 sites).
  Check the per-sample medians in `filter_counts.tsv`. If one sample's median is low, run a sensitivity
  check with a larger `MAX_DP_FACTOR`. The verdict does not depend on it.
- **Genetics can't say which label is wrong** (tumour or normal) or where the swap happened; that needs
  the lab records.

### 3.4 The mice are not genetically identical

Over 1,000 confident opposite homozygotes between animals means the KI colony is not isogenic. Possible
causes:
- residual 129 / ES-cell background or incomplete backcrossing;
- B6J/B6N substrain mixing;
- segregating heterozygosity in the colony.

This is real germline variation, and it inflates any comparison against another animal's normal, i.e.
settings B and C and any mis-paired A. The collaborator questions in the plan (ES-cell origin, backcross
generations, substrain) are needed.

### 3.5 FreeBayes is weak as a somatic caller

It is a germline caller with a crude prefilter:
- tumour GT must be 0/1, so it drops VAF ≈ 1 (LOH, clonal chrX in these male mice) and low-VAF calls;
- normal GT 0/0 still allows alt reads in the normal (e.g. 7/59).

In a 2-of-3 rule it casts the deciding vote wherever Strelka2 and Mutect2 disagree.

### 3.6 FreeBayes FILTER does not follow the provider's prefilter expression (2026-11-25)

Checked with `mouse_wgs_freebayes_pass_check.R` on the setting-A FreeBayes files
(`calls/*.freebayes.ann.norm.bcf`, fully split: no multiallelic ALT left). The expression

`INFO/MQM >= 40 && GT[normal]="0/0" && GT[tumor]="0/1" && FMT/DP[normal] >= 30 && FMT/DP[tumor] >= 30 && QUAL/INFO/AO >= 5 && SAF > 0 && SAR > 0`

was evaluated in R on every record (tumour/normal column from the sample names; a missing value fails) and
compared with FILTER:

| Tumour | records | FILTER PASS | expression TRUE | both | PASS, expression FALSE | expression TRUE, not PASS |
|---|---|---|---|---|---|---|
| HCC_1 | 1,272,624 | 55,611 | 54,890 | 24,579 | 31,032 | 30,311 |
| HCC_2 | 1,138,614 | 34,433 | 35,866 | 11,164 | 23,269 | 24,702 |
| HCC_3 | 1,158,627 | 34,672 | 33,674 | 10,076 | 24,596 | 23,598 |

- Only 29–44% of PASS records satisfy the expression, and about as many records satisfy it but are `.`.
- Among PASS records that fail it, the tumour genotype `0/1` fails most (HCC_3: 21,536 of 24,596), then
  SAR > 0, SAF > 0, QUAL/AO, MQM. Depth and the normal genotype did not show up as failing.
- PASS records that fail often come in groups with one QUAL at the same or nearby positions (e.g.
  10:3202200 T>TA and T>G; 10:3450091 CAG>C and C>CAG; 10:3239159 and 10:3239164), with tumour GT `0/0`.
- **Most likely explanation:** the filter was applied to the original FreeBayes records (multiallelic /
  complex) before they were split and normalised (`.norm`). The pieces kept the original FILTER, but GT, AO,
  SAF and SAR now describe the piece. Not yet confirmed: the script's `split_piece` table (same QUAL within
  50 bp) is still to be reviewed, and records that pass the expression but are `.` are not fully explained.
- **The expression as written cannot run as-is:** bcftools takes sample indices (`GT[0]`, `FMT/DP[1]`), not
  names, and the samples are not called `normal` / `tumor`. Ask Patrik for the exact command and the order of
  filter vs `norm`. A per-variant version, on the split file:
  `bcftools filter -s FB_FAIL -i 'INFO/MQM>=40 && GT[N]="ref" && GT[T]="het" && FMT/DP[N]>=30 && FMT/DP[T]>=30 && QUAL/INFO/AO>=5 && INFO/SAF>0 && INFO/SAR>0'`
  (N, T = normal and tumour index from `bcftools query -l`; "ref"/"het" also match phased or `1/0` genotypes).
- Consequence: the file's FreeBayes PASS is not the expression; it is **not** enough to add the expression
  on top of PASS (that drops the records that pass the expression but are `.`). Hence three FreeBayes sets in
  the settings comparison (2.2). Mutect2 and Strelka PASS are unaffected (their own filters).

---

## 4. Open decisions (ask Paula before changing)

- **Consensus definition.** Proposed:
  - main set = PASS in ≥ 2 callers' own files, or Strelka2 ∩ Mutect2 with FreeBayes only as an annotation;
  - sensitivity tier = 1 somatic caller PASS plus the other filtered only for low support (`weak_evidence`,
    `LowEVS`, `LowDepth`), never for `germline` / `normal_artifact` / `panel_of_normals`.
  - Not agreed yet.
- **Callers script.** Done 2026-11-25: every per-caller set in `mouse_wgs_snv_indel_callers.R` is PASS in the
  caller's own file; the merged file is only the separate "union (merged)" set and never filters a caller.
- **MNV comparison unit.** Decided 2026-11-25: MNVs split for caller comparison; doublets count once (DBS)
  for signatures (`scripts/analysis/README.md`, Open questions 1).
- **FreeBayes PASS definition.** Decided 2026-11-25: keep the file's PASS as `A_freebayes` and compare it with
  `A_freebayes_pass_expr` and `A_freebayes_expr` (settings script). Open: whether to add the same to
  `mouse_wgs_snv_indel_callers.R` as a switch (`FREEBAYES_PASS = "filter" | "pass_expr" | "expr"`, one run
  per value, separate output folders; changes FreeBayes, >=2 callers, all 3; not Mutect2, Strelka,
  Mutect2 & Strelka or the merged union). Not done; Paula to decide.
- **High-confidence set** stays Mutect2 & Strelka (PASS in both own files); unaffected by the FreeBayes issue.

## 5. Next steps, in order

0. **Settings / FreeBayes (2026-11-25):** reinstall qcVCF on ssi (2.2); run the 8-set
   `mouse_wgs_snv_indel_settings.R`; rerun `mouse_wgs_freebayes_pass_check.R` and review the `split_piece`
   table; run the union chat code for one tumour; ask Patrik for the exact FreeBayes filter command and
   whether it ran before or after `bcftools norm`; then decide on the callers-script switch (section 4).

1. Review the **common-sites gtcheck rows**: `grep -E '^(site_set|common_sites)' results/mouse_wgs/gtcheck/gtcheck_table.tsv | column -t`.
2. **Upgrade bcftools** in r-mutation (`"bcftools>=1.17"`), then re-run `gtcheck_freebayes.sh` to get
   per-sample genotype masking.
3. Run the current **`mouse_wgs_sample_identity.R`**. Review `verdict.tsv`, `ibs0_rate_per_tumour.png`,
   `ibs0_genome_distribution.png` / `ibs0_blocks.tsv` and `quality_ibs0_vs_other.png`.
4. ~~Run `mouse_wgs_somalier.R`~~ done 2026-11-25 (result in 3.3). Still to review: `robustness.tsv`,
   `assignment_scores.tsv`, the plots, and the max-depth filter (3.3).
5. Download MGP REL-1505 to `/common/WORK/pstancl/references/mgp/` and re-run
   **`mouse_wgs_identity_strain_check.R`** (commands in 6.4). Review `summary_discordant_sites.tsv`,
   `discordant_sites_genome.png` and `isolated_discordant_sites.tsv`; then send
   `bam_check_sites_for_patrik_annotated.csv` to Patrik (he fills in the `BAM_reads_ref_alt_*` columns).
   Re-run `mouse_wgs_snv_indel_callers.R` for the revised plots (6.1); run `mouse_wgs_snv_indel_settings.R`
   and `mouse_wgs_freebayes_pass_check.R`.
6. **Ask Patrik and the lab** where the KI1/KI2 swap happened and what the `id_` numbers mean; ask about
   HCC_3/KI3.
7. **Re-call setting A** with HCC_1 vs KI2 and HCC_2 vs KI1 (Strelka2, Mutect2, FreeBayes, merge). Compare
   the old and corrected burden, VAF, SBS96 and germline overlap; the false positives should drop.
8. Then continue the plan: Task 1 (`mouse_wgs_snv_indel.R`, after switching annotation to CSQ), Task 2
   (cross-sample evidence from FORMAT fields), Tasks 3–6.

Not possible without BAMs: pileups, callable regions and mutations per Mb, SV read support, copy number,
new caller runs for strategies D (GM normal vs WT) and E (WT vs WT).

---

## 6. Agreed and done on 2026-11-25

### 6.1 Callers in setting A (`mouse_wgs_snv_indel_callers.R`)

Agreed rules (do not change without Paula):
- **Per-caller PASS = FILTER == PASS in that caller's own file**, by its own internal filters. The merged
  file (`union.ann.bcf`) is a separate set, labelled **"union (merged)"**, and never filters, defines or
  falls back for a caller's variants. The script stops if the single-caller files are missing.
  Why: taking caller sets from the merged FILTER dropped variants (HCC_1 Strelka: 28,817 own-file PASS vs
  21,439 plotted before the fix; Mutect2 6,838 PASS).
- **MNVs are split into SNVs for every caller comparison** (counts, pairwise, intersects, Venn). Signatures
  count a doublet once (DBS); doublet halves are left out of SBS96. Why: Strelka2 writes a doublet as two SNVs,
  Mutect2/FreeBayes as one MNV, and FreeBayes also joins non-adjacent SNVs into one haplotype record
  (HCC_1 FreeBayes PASS: 55,611 records = 29,234 SNV + 1,864 MNV + 24,513 indel; the MNVs hold 4,386 SNVs).
  Also recorded in `scripts/analysis/README.md`, Open questions 1.
- **Sets** (all from own-file PASS): each caller; `>=2 callers` (PASS in ≥ 2 own files, any pair);
  **`mutect2+strelka`** (PASS in Mutect2's and in Strelka's own file, then present in both; FreeBayes ignored,
  because Paula distrusts FreeBayes; the high-confidence set); `all 3 callers`; and `union (merged)` (merged
  FILTER PASS, includes 1-caller calls).
- **MNV split everywhere** now also covers the R7 intersects, the Venn and `intersect_mutect2_strelka_PASS.tsv`
  (6963be1). Before, the Venn did not split: HCC_1 Mutect2 & Strelka was 3,291 in the Venn vs 3,351 in the
  pairwise plot; the difference was doublet halves.
- **Why the pairwise % dropped** versus older plots (e.g. HCC_1 Mutect2 → Strelka 94.6 % → 78.4 %): the old
  plots predate 391d39a and took per-caller sets from merged PASS + INFO/CALLER, which inflates sharing. The
  current numbers are correct; `parameters.tsv` `PER_CALLER_SOURCE` shows which source a run used.

Changes to the outputs:
- `intersect_mutect2_strelka_PASS.tsv` and `_counts.tsv`; the set also appears in the raw-vs-merged bar
  plot, the SBS96 sets and both cosine heatmaps.
- **Venn** with three sets: my intersect (≥ 2 callers PASS), Mutect2 & Strelka (PASS in both), merged file
  (PASS, ≥ 2 callers) (`venn_snv_my_intersect_vs_mutect2_strelka_vs_merged.png/.tsv`).
- **counts_per_caller** plot: all records vs PASS as rows (`counts_per_caller_all_vs_PASS.tsv`).
- **qcVCF pairwise shared:** one plot per class (SNV, INDEL) across all tumour × caller sets
  (`HCC_1 mutect2`, `HCC_2 strelka`, ...) instead of one plot per tumour.
- **SBS96:** qcVCF plots stacked in one column, each set split into all SNVs / PASS / not PASS
  (`sbs96_snvs_by_filter.tsv`); the ggplot `sbs96_<tumour>` shows PASS only, one row per set. The SBS96
  matrix and cosine similarities stay PASS-only. Contexts come from the GRCm38 FASTA (samtools faidx), folded
  to the pyrimidine strand, raw counts, no rescaling.
- **Between-tumour cosine** (`sbs96_cosine_between_tumours.tsv/.png`, also in `summary.md` and the rds as
  `cosine_samples`): HCC_1/2/3 per set; colour scale from the lowest off-diagonal value (rounded down to
  0.05) to 1, diagonal grey.
- Commits: 391d39a 67b47d7 65de2c2 4159641 72191e9 657c29f a81c9d4 536528f 6dc6393 ecb9d69 ba8fed4 b795527
  beb1a5e 3331a41 6963be1. Re-run needed for all of it.

FreeBayes checks done in the R console:
- The provider's prefilter is the rule above (section 1). The check that every FreeBayes PASS satisfies it,
  and the reverse, is now `mouse_wgs_freebayes_pass_check.R` (not run yet).
- FreeBayes sample names differ from the other callers; read them with `bcftools query -l` (normal first).
- The BCFs are not indexed: query sites with `-i 'CHROM=="10" && POS==3184122'`, not `-r`.

### 6.2 Sample identity: how the evidence is read

- **Primary evidence:** opposite homozygotes (IBS0: one sample 0/0, the other 1/1). One animal cannot have
  them except by error, LOH or a somatic event, so they separate animals best.
- **Complementary evidence:** heterozygous vs homozygous differences (e.g. chr1:5609923). They also arise from
  LOH, copy number or allelic imbalance in a tumour, so they support but do not decide.
- **Genotype bins on VAF:** ≤ 0.05 hom-ref, 0.35–0.65 het, ≥ 0.95 hom-alt, anything else ambiguous.
- **The ascertainment critique** (sites chosen through each pair's FreeBayes file) does not overturn the
  result: the 0 discordance for the genetic matches and 15–23 % for the labelled pairs are genotypes at the
  sites, not which sites were kept. Checked on a tumour-independent site set (all three KI normals non-ref);
  the README and `mouse_wgs_somalier.R` header wording was corrected.

### 6.3 Informativeness and BAM check (`mouse_wgs_identity_informativeness.R`, run)

- Outputs in `results/mouse_wgs/mouse_wgs_identity_informativeness/`: genotype patterns across KI1–3,
  informativeness summary, opposite homozygotes per pair of normals, private sites, tumour-dependent
  retention, discordance on the tumour-independent site set, sites per chromosome.
- **`bam_check_sites_for_patrik.csv`:** 5 opposite-homozygous sites per pair of normals (`N_PER_TEST`), one per
  chromosome, deepest first, plus the example sites chr10:3184122 and chr1:5609923. Patrik fills in the
  `BAM_reads_ref_alt_*` and comment columns. Five per pair is enough for an independent read-level check;
  raise `N_PER_TEST` for more.

### 6.4 Are the differences germline? (`mouse_wgs_identity_strain_check.R`)

Question: could the tumour-vs-labelled-normal differences be real somatic mutations instead of another
animal's germline? Three checks per panel A site:
- **Mouse Genomes Project** (MGP REL-1505, GRCm38, dbSNP142 rsIDs, 36 strains): known polymorphism with the
  same ALT, which strains carry it, C57BL/6NJ and the 129 strains (donor blocks of the knock-in line).
  EBI has no REL-2005; REL-2112 is GRCm39 and must not be used.
- **WT normals** (WT1–3) from the setting-B Mutect2 files.
- **Neighbouring SNPs** (± 5 Mb, ≥ 3 neighbours): inherited differences come in blocks where the tumour
  matches another mouse's normal; a real somatic mutation is isolated. Isolated sites are split into
  "matches another mouse's normal" and candidates whose tumour genotype is in no KI normal
  (LOH-compatible vs tumour gained an allele).
- Not in MGP does not mean somatic (private germline variants exist).

Status:
- The first run (20.5 min) read MGP over ftp and got **0 MGP records**; the MGP columns are NA and the
  summaries show n/a. A single-site ftp query works (1:3000287 → rs587313017), cause unknown.
- Agreed: download MGP locally. The script uses
  `/common/WORK/pstancl/references/mgp/mgp.v5.merged.snps_all.dbSNP142.vcf.gz` (+ `.tbi`) when present,
  else ftp. Delete the cached `mgp_panelA_positions*.vcf.gz` in the output folder before re-running.

```
cd /common/WORK/pstancl/projects/MariaBoskovic/SPRTN
git pull
micromamba activate /common/WORK/pstancl/envs/r-mutation
rm -f results/mouse_wgs/mouse_wgs_identity_strain_check/mgp_panelA_positions*.vcf.gz*
mkdir -p /common/WORK/pstancl/references/mgp
nohup bash -c '
  cd /common/WORK/pstancl/references/mgp &&
  wget -c -q https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz{,.tbi,.md5} &&
  md5sum -c mgp.v5.merged.snps_all.dbSNP142.vcf.gz.md5 &&
  cd /common/WORK/pstancl/projects/MariaBoskovic/SPRTN &&
  Rscript scripts/mouse_wgs/mouse_wgs_identity_strain_check.R
' > results/mouse_wgs/strain_check.log 2>&1 &
```
Then check `grep -n MGP results/mouse_wgs/strain_check.log` and `summary_discordant_sites.tsv`.

### 6.5 Still open from today

- Mutect2 low-ROQ artefacts: proposed ROQ < 30 filter and extra tumourevo runs on filtered calls; not agreed.
- Consensus definition for the main somatic set (section 4).
- Optional: recolour the IBS0-vs-IBS2 plot so IBS0 and IBS2 are labelled directly (offered, not confirmed).
