# Mouse WGS (SPRTN KI HCC) - scripts, findings, status

Author: Paula Breyer. Last updated: 2026-11-24.

Read **Status and next steps** first when the work resumes.

## Data

Location on ssi: `/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/data/mouse/wgs/`. Patrik processed
the data; his notes are in `patrik_mouse_wgs_processed.txt` and the md5 lists are in `md5sums*.txt`.
There are no BAMs, so there is no pileup, no callable territory (counts only, not mutations per Mb),
no SV read support and no copy number.

| Mouse | Tumour | KI normal (non-adjacent liver), as labelled |
|---|---|---|
| HCC_1 | sfb_dikic_2025_05_10_KI_HCC_1 | sfb_dikic_2025_05_05_KI1 |
| HCC_2 | sfb_dikic_2025_05_11_KI_HCC_2 | sfb_dikic_2025_05_06_KI2 |
| HCC_3 | sfb_dikic_2025_05_12_KI_HCC_3 | sfb_dikic_2025_05_07_KI3 |

There are also three WT normal livers. They appear only in settings B and C.

| Setting | Folder | Normal | Callers |
|---|---|---|---|
| A | `tumor_vs_normal/snv_indel/` | matched KI normal | Strelka2, Mutect2, FreeBayes. Single-caller files plus `merged/*.union.ann.bcf` (INFO/CALLER) |
| B | `tumor_vs_joint_wt/` | the 3 WT normals, joint | Mutect2 (exploratory) |
| C | `tumor_only_pon/` | none; WT panel of normals | Mutect2 |
| SV | `tumor_vs_normal/sv/calls/`, `sv/merged/` | matched KI normal | Delly, SvABA, Manta; OctopuSV merge |

- **Normalisation:** small variants were split (`bcftools norm -m -any`) but not left-aligned, so the
  same indel can be written differently by different callers.
- **Reference:** GRCm38 FASTA (Ensembl 102), from the cross-species work:
  `/common/WORK/pstancl/references/cross_species/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa`.
- **Annotation:** VEP CSQ, not SnpEff ANN.
- **Formats:** most files are BCF, which only bcftools (or R's VariantAnnotation) can read. `zcat` does
  not work on them.

## Scripts (`scripts/mouse_wgs/`)

Output goes to `results/mouse_wgs/<script name>/`, and each script writes a RUN_LOG.md. Env:
`r-mutation`, which has bcftools; if RStudio doesn't see the env's PATH, the scripts fall back to
`/common/WORK/pstancl/envs/r-mutation/bin`.

| Script | What it does | State |
|---|---|---|
| `00_input_inspection.R` | Task 0: provenance, md5, headers, normals used per setting, records/FILTER, SV inventory, manifest, issues | runs; fixed 2026-11-24 (setting labels A/B/C, headers without `#CHROM`); summary not yet reviewed |
| `mouse_wgs_snv_indel_callers.R` | Setting A callers: R1-R7 compare the single-caller files with the merged file; blocks 4-8 do counts, pairwise sharing, exclusive categories, union vs callers, SBS96 | R1-R7 run and reviewed; **blocks 4-8 need fixing** (see below) |
| `mouse_wgs_snv_indel.R` | Task 1: settings A/B/C - counts, Ts/Tv, FILTER reasons, VAF/depth, overlap between settings, annotation | written; **block 9 reads SnpEff ANN but the files have VEP CSQ** |
| `mouse_wgs_example_sites.R` | Diagnostics for chr19:30419577 T>C and chr10:3184122 G>T from the single-caller files (all nine samples) | run; chr10 led to the sample swap |
| `mouse_wgs_sample_identity.R` | Tumour/KI normal identity: inferred genotypes, concordance, IBS0, LD-thinned rates, optional MGP strains, bcftools gtcheck | run; found the swap |
| `mouse_wgs_somalier.R` | Identity with somalier, using a mouse germline site list built from the FreeBayes calls | **not yet run**: somalier binary missing (see below) |
| `gtcheck_freebayes.sh` | Identity with bcftools gtcheck, FreeBayes SNVs, DP/GQ masking, all-sites and common-sites runs | written 2026-11-24; run in progress |

Run order: `00_input_inspection.R`, then `mouse_wgs_snv_indel_callers.R`, then `mouse_wgs_snv_indel.R`. The identity scripts
are independent of that chain.

## Findings so far

### 1. The merged setting A file can't give per-caller PASS sets (confirmed 2026-11-24)

Checked on all 6.8 M merged records with R5c in `mouse_wgs_snv_indel_callers.R`; there were no exceptions.

- **Merged FILTER:** follows the bcftools merge default. A record is PASS if at least 1 caller is PASS
  and no caller that found it has a named filter. FreeBayes writes `.` when its prefilter fails, so it
  never blocks PASS.
- **INFO/CALLER:** lists every caller that called the variant, passed or not. So "merged PASS with ≥2
  callers" is not the same as "≥2 callers PASS".
- **Consequence:** caller sets, intersects and the union must come from the single-caller files' own
  PASS. Use the merged file only to look up which callers called a variant.
- **Examples:**
  - HCC_1 10:3184122 is merged PASS with `strelka,freebayes`, but FreeBayes is `.`, so only Strelka passed it.
  - HCC_1 19:30419577 is FreeBayes PASS but merged LowEVS.

### 2. Sample swap in setting A: HCC_1 and HCC_2 were called against the wrong normal

Found with `mouse_wgs_sample_identity.R`, using FreeBayes genotypes (DP ≥ 30; ref VAF ≤ 0.05, alt ≥ 0.95,
het 0.35-0.65).

| Pair | Concordance | Opposite homozygotes | Reading |
|---|---|---|---|
| HCC_1 - KI2 | 0.988 | 60 | same animal |
| HCC_2 - KI1 | 0.982 | 73 | same animal |
| HCC_1 - KI1 (labelled) | 0.734 | 1638 | different mice |
| HCC_2 - KI2 (labelled) | 0.737 | 1402 | different mice |
| HCC_3 - KI3 (labelled) | 0.951 | 190 | not confirmed |

- HCC_3 - KI3 has 190 opposite homozygotes, close to the counts between mouse 3 and the HCC_1/KI2
  animal (201-232). Mouse 3 and that animal look closely related.
- Consequence: setting A calls for HCC_1 and HCC_2 include germline differences between two mice,
  which explains HCC_1's inflated counts. Treat setting A results for HCC_1 and HCC_2 as invalid until
  they are re-called with the right pairs. Settings B and C don't use the KI normals, so they aren't
  affected.
- somalier and gtcheck are being run as second methods. They use the same FreeBayes calls, so they
  are not independent evidence.

## Status and next steps

**Waiting on others**
1. Paula to confirm the source of the swap with Patrik and the lab: is the label or the sample wrong,
   and is HCC_3 - KI3 correct?
2. Setting A re-called with the corrected pairs (HCC_1 vs KI2, HCC_2 vs KI1) by Patrik. Until then, do
   not interpret the HCC_1 and HCC_2 caller comparisons.

**Identity checks (finish these first)**

3. **somalier:** the script expects the binary at `/common/WORK/pstancl/tools/somalier/somalier`
   (v0.3.5, installed by the osteosarcoma `00_somalier.sh`). Check it with `--version`. If it's gone,
   reinstall with the line in the script header, then run `mouse_wgs_somalier.R`.
4. **gtcheck:** collect the result of `gtcheck_freebayes.sh` (`results/mouse_wgs/gtcheck/gtcheck_table.tsv`).
   Use the `common_sites` rows.
5. **Record the outcome:** write the result of 3 and 4 into "Findings" section 2 above.

**Fixes before the analysis continues**

6. In `mouse_wgs_snv_indel_callers.R`, block 3b builds `per_caller` from the merged file. Change it to
   use the single-caller files' PASS, with multiallelics split as in R2. Then:
   - base R7 on `atom` / `memb`, so MNVs are split as in blocks 4-6. R7 now classifies by
     `nchar(REF) == nchar(ALT)` without splitting.
   - rewrite the block 7 `union_check` comment: once the per-caller sets come from the single-caller
     files, 0 / 0 is no longer the expected result.
   - R1-R6 stay as the "check the merged file" section.

   Proposed on 2026-11-24; waiting for Paula's OK.
7. `mouse_wgs_snv_indel.R` block 9: switch from SnpEff ANN to VEP CSQ, with the field order taken from
   the `##INFO=<ID=CSQ` header.
8. Run `00_input_inspection.R` end to end and read its `summary.md` and `issues.tsv`.

**Analysis (after 1-8)**

9. Caller comparison and SBS96 per tumour (`mouse_wgs_snv_indel_callers.R`). Use only HCC_3 until
   setting A is re-called.
10. Settings A/B/C comparison (`mouse_wgs_snv_indel.R`). B and C are valid for all three tumours.
11. SVs: inventory is in `00_input_inspection.R`; no analysis script yet.

**Open decisions (ask before changing)**
- Comparison unit: split MNVs into SNVs (current) or join doublets. Undecided since 2026-10-07.
- Mutect2 low-ROQ artefacts / ROQ filter: carried over from the human analysis, not applied here.
