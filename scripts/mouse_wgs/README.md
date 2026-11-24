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
| `mouse_wgs_somalier.R` | somalier with a mouse germline site list built from the FreeBayes calls | written; **somalier not installed** |
| `mouse_wgs_snv_indel.R` | Task 1: settings A/B/C (counts, VAF, FILTER, overlaps, annotation) | written, not run; block 9 reads ANN and must switch to **CSQ** |

Run commands (from the repo on ssi):
```
cd /common/WORK/pstancl/projects/MariaBoskovic/SPRTN
git pull
Rscript scripts/mouse_wgs/<script>.R          # or: bash scripts/mouse_wgs/gtcheck_freebayes.sh
```

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
4. **Install somalier** (`micromamba create -y -p /common/WORK/pstancl/envs/somalier -c conda-forge -c bioconda somalier`),
   then run `mouse_wgs_somalier.R`.
5. **Ask Patrik and the lab** where the KI1/KI2 swap happened and what the `id_` numbers mean; ask about
   HCC_3/KI3.
6. **Re-call setting A** with HCC_1 vs KI2 and HCC_2 vs KI1 (Strelka2, Mutect2, FreeBayes, merge). Compare
   the old and corrected burden, VAF, SBS96 and germline overlap; the false positives should drop.
7. Then continue the plan: Task 1 (`mouse_wgs_snv_indel.R`, after switching annotation to CSQ), Task 2
   (cross-sample evidence from FORMAT fields), Tasks 3–6.

Not possible without BAMs: pileups, callable regions and mutations per Mb, SV read support, copy number,
new caller runs for strategies D (GM normal vs WT) and E (WT vs WT).
