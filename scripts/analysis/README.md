# Mutation analysis (RJALS)

Scripts that analyse the pipeline outputs under `results/wgs/` (sarek, oncoanalyser, tumourevo).
Results go to `results/wgs/analysis/`. Caller-specific filters and fields are in [CALLERS.md](CALLERS.md).

## Run order and environments

Envs are under `/common/WORK/pstancl/envs/` (recipes: `env_*.yml` in this folder).
Run `../wgs/06_normalize_vcfs.sh` once before, so indels are compared after `bcftools norm`.

| step | script | env | what it does |
|---|---|---|---|
| 1 | `01_sarek.R` | r-mutation | sarek: SNV/indel QC (raw + PASS), Mutect2/Strelka2/MuSE overlap, 96-context, SV, CN, signature inputs |
| 1 | `02_oncoanalyser.R` | r-mutation | oncoanalyser: SAGE/PURPLE QC, 96-context, ESVEE/LINX SV, PURPLE CN, signature inputs |
| 1 | `03_compare_callers.R` | r-mutation | all four SNV callers: overlap, VAF agreement, set QC (Wilcoxon), 96-context per caller and consensus, signature inputs per set |
| 1 | `04_summary_table.R` | r-mutation | sarek vs oncoanalyser results table (`summary/summary_table.csv/.md`) |
| 2 | `01b_sarek_signatures.py [--pipeline sarek\|oncoanalyser\|comparison]` | sigprofiler | SigProfilerMatrixGenerator + SigProfilerAssignment `cosmic_fit` (SBS96, DBS78, ID83), one fit per set |
| 3 | `01c_sarek_signatures_organ.R <pipeline>` | signature-tools (`Rscript --no-environ`) | FitMS, liver common + rare signatures (signature.tools.lib; needs NNLM) |
| 4 | `05_signature_summary.R <pipeline>` | r-mutation | attributions of all sets side by side (SigProfiler + FitMS), fit quality |
| 4 | `05b_hartwig_sigs_check.R` | r-mutation | SAGE set: Hartwig SIGS vs our SBS96 counts, and SIGS / SigProfiler / FitMS signatures |

qcVCF (private repo PaulaStancl/qcVCF) must be installed in r-mutation, otherwise the qcVCF
plots (96-context, pairwise shared) are skipped with a message. Install from a clone:
`git clone https://github.com/PaulaStancl/qcVCF.git` then `Rscript -e 'remotes::install_local("qcVCF")'`.

## Method decisions

**What counts as a call.** FILTER = PASS for Mutect2, Strelka2 and SAGE/PURPLE. MuSE: PASS
and Tier1-5, as its authors recommend for WGS (CALLERS.md; `SPRTN_MUSE_TIERS` restricts tiers).

**MNVs are split for comparing callers.** SAGE and Mutect2 report adjacent substitutions on the
same reads as one MNV record; Strelka2 and MuSE write them as separate SNVs. For overlap,
concordance, caller totals, set QC and the summary table, each MNV is split into its SNVs
(`atomize_mnv()`), so e.g. SAGE's 2,124 SNV + 17 MNV records give 2,158 SNVs to match.

**Doublets count once in mutational-signature analyses (since 2026-10-07).** An SNV with
another SNV at the neighbouring position in the same caller or set is half of a doublet base
substitution (DBS): a split MNV, or the two SNVs Strelka2/MuSE write for the same event
(`in_doublet()` in `00_setup.R`). Doublets are left out of SBS96 (96-context plots, substitution
spectrum, SigProfiler and FitMS SBS96) and counted in DBS78 only. This follows COSMIC/PCAWG
(Alexandrov et al. 2020, Nature), ICAMS (`SplitSBSVCF` removes both SBSs of a DBS) and Hartwig
SIGS (single-base SNVs only). SigProfilerMatrixGenerator (1.3.6) would count a doublet in both
DBS78 and SBS96, so `01b` builds SBS96 from `signatures/input/vcf_sbs/` (doublet halves
removed) and DBS78/ID83 from `signatures/input/vcf/` (all calls).
Methods wording: *"SNVs adjacent to another SNV in the same call set were treated as doublet
base substitutions (DBS78) and excluded from the SBS96 catalogue."*

**Consensus sets.** sarek SNVs: PASS in >= 2 of Mutect2, Strelka2, MuSE; sarek indels:
Mutect2 AND Strelka2 (MuSE calls no indels). Shared with oncoanalyser: sarek consensus AND
SAGE. Matched on CHROM:POS:REF:ALT after splitting MNVs and normalising indels.

**Signatures.** One fit per set (no de novo extraction - one tumour). SigProfilerAssignment
`cosmic_fit` with COSMIC reference signatures; FitMS with liver-specific common signatures
(tier T1), rare tier T2, KLD, bootstrap (`SPRTN_NBOOT`, default 200). oncoanalyser's own SIGS
fits the 30 COSMIC v2 signatures by least squares; `05b` matches it to the others on the
signature number. Plot colours in `05` follow SigProfiler's `plotActivity` palette (COSMIC
artefact signatures in grey).

## Changes

**2026-10-07**
- MuSE (extra sarek run `RJALS_vc`) added to QC, overlaps and the summary table; PASS + Tier1-5 count as calls.
- SV and CN outputs moved to `sarek|oncoanalyser/{sv,cnv}/`.
- Plots: `facet_grid(space = "free")` where panels have different categories; bars with `position_dodge(preserve = "single")`.
- 96-context per caller: rows by number of callers; `snv_96context_per_caller` shows each caller, then all callers, then >= 2.
- Concordance plots and `snv_indel_caller_totals.csv`: each caller's total.
- 03 set QC: VAF, ALT reads, depth, normal VAF/depth and caller score for this-caller-only / >= 2 / all callers, Wilcoxon tests (BH-adjusted) in `set_qc_tests.csv`; each caller's top 1% left out of the plots only.
- 04: sarek SNV consensus = >= 2 of Mutect2, Strelka2, MuSE.
- Signatures per comparison set (muse, mutect2, strelka, sage, two_plus, all_callers); `05_signature_summary.R` added; SigProfiler colours.
- `01c` stops early with install instructions when NNLM is missing (FitMS needs it).
- `05b_hartwig_sigs_check.R` added (SAGE set vs Hartwig SIGS).
- Doublets counted once (see Method decisions): `vcf_sbs/` signature inputs, 01b runs the matrix generator twice (`matrix_generator/` for SBS96, `matrix_generator_all/` for DBS78 + ID83). Expected effect: SAGE SBS96 2,158 -> 2,124, identical to SIGS. Rerun 01-03, 01b, 01c, 05, 05b.
