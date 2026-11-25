# Mouse WGS (Sprtn knock-in HCC): status, findings and how to continue

Last updated 2026-11-24. Plan: Claude Docs "Mouse WGS → human comparison: analysis plan"
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
| `mouse_wgs_snv_indel_callers.R` | Callers within setting A. Blocks 3–8: counts, pairwise sharing, 1/2/3-caller categories, merged union, SBS96 (GRCm38 FASTA from h2m resources). R1–R7: raw single-caller PASS → intersect vs merged file, merged FILTER rule check (R5c), Venn and bar plot of my intersect vs merged | run; R1–R7 reviewed |
| `mouse_wgs_example_sites.R` | Diagnostics for chr19:30419577 T>C and chr10:3184122 G>T from the single-caller files: per-caller metrics, FreeBayes criteria, per-caller read counts in all samples, plots 01–08 | run, reviewed |
| `mouse_wgs_sample_identity.R` | Sample identity from FreeBayes genotypes: IBS0 rates per shared confident sites (WT from Mutect2 joint), LD-thinned rate, IBS0 blocks and chromosomes, depth and GQ/MQM at IBS0 sites, optional MGP strains, bcftools gtcheck | v1 run (counts in section 3); **current version not run yet** |
| `gtcheck_freebayes.sh` | bcftools gtcheck: biallelic autosomal SNVs, QUAL/DP/GQ filter, merge `-m none`, all-sites and common-sites runs, per-1000 table | run; common-sites rows not reviewed yet |
| `mouse_wgs_somalier.R` | Identity with three methods (somalier, bcftools gtcheck, custom IBS0) on one masked 6-sample genotype file; panels A/B/C x thinning; see 2.1 | rewritten 2026-11-24 (94c7d76); `--test` passes on ssi 2026-11-25 (all checks incl. somalier, 9a19573); **real data not run yet** |
| `mouse_wgs_snv_indel.R` | Task 1: settings A/B/C (counts, VAF, FILTER, overlaps, annotation) | written, not run; block 9 reads ANN and must switch to **CSQ** |

Run commands (from the repo on ssi):
```
cd /common/WORK/pstancl/projects/MariaBoskovic/SPRTN
git pull
Rscript scripts/mouse_wgs/<script>.R          # or: bash scripts/mouse_wgs/gtcheck_freebayes.sh
```

### 2.1 `mouse_wgs_somalier.R` (rewritten 2026-11-24)

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
  - tumours are never used to select sites. Their ambiguous genotypes are masked (`./.`), never set to hom-ref.
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
- **somalier agrees on all three tumours.** gtcheck was skipped in this run because of a bug in how the
  script checked the gtcheck options, fixed in be67ed8; it is to be re-run.
- **Still to review:**
  - `robustness.tsv` (panels B/C and the coarser thinning);
  - `assignment_scores.tsv`;
  - the plots in `results/mouse_wgs/mouse_wgs_somalier/`.
- All three methods use the same FreeBayes calls. Lab records remain the independent confirmation.

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

---

## 4. Open decisions (ask Paula before changing)

- **Consensus definition.** Proposed:
  - main set = PASS in ≥ 2 callers' own files, or Strelka2 ∩ Mutect2 with FreeBayes only as an annotation;
  - sensitivity tier = 1 somatic caller PASS plus the other filtered only for low support (`weak_evidence`,
    `LowEVS`, `LowDepth`), never for `germline` / `normal_artifact` / `panel_of_normals`.
  - Not agreed yet.
- **Callers script.** Blocks 4–8 of `mouse_wgs_snv_indel_callers.R` still take caller sets and the union
  from the merged file. Switching them to the raw route was proposed but not agreed.
- **MNV comparison unit.** Split MNVs vs join doublets is still open (memory, 2026-10-07).

## 5. Next steps, in order

1. Review the **common-sites gtcheck rows**: `grep -E '^(site_set|common_sites)' results/mouse_wgs/gtcheck/gtcheck_table.tsv | column -t`.
2. **Upgrade bcftools** in r-mutation (`"bcftools>=1.17"`), then re-run `gtcheck_freebayes.sh` to get
   per-sample genotype masking.
3. Run the current **`mouse_wgs_sample_identity.R`**. Review `verdict.tsv`, `ibs0_rate_per_tumour.png`,
   `ibs0_genome_distribution.png` / `ibs0_blocks.tsv` and `quality_ibs0_vs_other.png`.
4. Run **`mouse_wgs_somalier.R`**: `--test` first, then the real data (commands in 2.1). somalier is already
   installed at `/common/WORK/pstancl/tools/somalier/somalier` (v0.3.5). Review `verdict.tsv`, `robustness.png`
   and the heatmaps, then add the result to section 3.3.
5. **Ask Patrik and the lab** where the KI1/KI2 swap happened and what the `id_` numbers mean; ask about
   HCC_3/KI3.
6. **Re-call setting A** with HCC_1 vs KI2 and HCC_2 vs KI1 (Strelka2, Mutect2, FreeBayes, merge). Compare
   the old and corrected burden, VAF, SBS96 and germline overlap; the false positives should drop.
7. Then continue the plan: Task 1 (`mouse_wgs_snv_indel.R`, after switching annotation to CSQ), Task 2
   (cross-sample evidence from FORMAT fields), Tasks 3–6.

Not possible without BAMs: pileups, callable regions and mutations per Mb, SV read support, copy number,
new caller runs for strategies D (GM normal vs WT) and E (WT vs WT).
