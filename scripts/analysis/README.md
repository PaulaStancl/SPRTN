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
| 5 | `06_drivers.R` | r-mutation | driver candidates: PURPLE + LINX vs tumourevo (IntOGen HCC), per gene and per mutation, callers per driver mutation, HCC watchlist incl. TERT promoter; slide tables `drivers_somatic_table` and `drivers_germline_table` (.csv/.md, germline incl. SPRTN) |

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
Not final - see Open questions 1.

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

**Drivers (06).** PURPLE's driver catalogue (Hartwig gene panel, driver likelihood, AMP/DEL)
and LINX (disruptions, reported fusions) are compared with tumourevo's annotation, which only
flags Mutect2 PASS mutations in IntOGen HCC driver genes with a MODERATE/HIGH VEP impact (the
gene is a known driver gene, not the variant a known driver). A curated HCC watchlist
(CTNNB1, TP53, AXIN1, ARID1A/2, ... and the TERT promoter, which is non-coding and never flagged
by tumourevo) lists every PASS mutation there. Knowledge-base annotation (OncoKB, CIViC,
ClinVar, hotspots) is a later step.

**Signatures.** One fit per set (no de novo extraction - one tumour). SigProfilerAssignment
`cosmic_fit` with COSMIC reference signatures; FitMS with liver-specific common signatures
(tier T1), rare tier T2, KLD, bootstrap (`SPRTN_NBOOT`, default 200). oncoanalyser's own SIGS
fits the 30 COSMIC v2 signatures by least squares; `05b` matches it to the others on the
signature number. Plot colours in `05` follow SigProfiler's `plotActivity` palette (COSMIC
artefact signatures in grey).

**FitMS: two plots, one fit.** 01c fits the liver organ signatures once per set; the cosine and
% unassigned belong to that fit. `05` shows it twice: `attribution_FitMS_liver_organ_signatures_unassigned`
(organ signatures + unassigned, shares of all SNVs - for QC) and
`attribution_FitMS_liver_as_reference_signatures` (the same exposures converted to reference
signatures with `convertExposuresFromOrganToRefSigs`, unassigned left out, shares of the assigned
SNVs - for interpretation next to SigProfiler's COSMIC fit). RefSig names match COSMIC processes
but the profiles are not identical. Report RefSig plus the % unassigned (`fit_quality.csv`).

**Germline drivers.** PURPLE's germline catalogue covers Hartwig's germline panel only (cancer
predisposition genes, ClinVar pathogenic / loss of function). Read each row by: (1) is it truly
pathogenic (ClinVar classification and stars); (2) is the gene relevant to this tumour; (3) did
the tumour hit the other allele (`biallelic`, LOH, or tumour VAF above the normal's); (4) does the
tumour show the phenotype (CHORD HRD / SBS3 for HR genes, MSI for MMR). SPRTN is not on the panel;
it is the project's own germline finding (compound heterozygous Y117C + c.718_718+3del, in
trans by tumour allelic imbalance; Ruijs-Aalfs syndrome, early-onset HCC) and is added to
`06`'s `drivers_germline_table` from our analysis (`../wgs/11_sprtn_checks.sh`, `12_whatshap_sprtn.sh`).
Germline findings are research results: confirm in an accredited lab and with genetic
counselling before any clinical use.

**Mutect2 read-orientation artefacts (ROQ).** `ROQ` (Mutect2 INFO) is the Phred-scaled quality
that the ALT allele is *not* a read-orientation artefact: damage on one DNA strand before
sequencing (8-oxoG oxidation -> C>A/G>T; cytosine deamination -> C>T/G>A) puts the ALT on reads
of one pair orientation only (F1R2 vs F2R1). sarek runs `LearnReadOrientationModel`, which learns
per-trinucleotide-context artefact priors from this sample's F1R2 counts; `FilterMutectCalls`
combines all its error models and filters at the probability threshold that maximises the
estimated F score - there is no fixed ROQ cut-off, so low-ROQ calls can PASS (Benjamin et al.
2019). Analogues in the other callers are weaker: Strelka2 has strand bias (`SNVSB`) inside its
`SomaticEVS` score, SAGE a strand-bias filter and base-quality recalibration per context, MuSE none.
In RJALS (`01` section 5b; `qc_vcf/pass/mutect2_roq_qc.pdf`, `mutect2_roq_spectrum.*`,
`mutect2_roq_cutoff.*`, `snv_96context_mutect2_by_ROQ.*`):
- 2,100 of 5,393 Mutect2 PASS SNVs have ROQ < 20: median VAF 0.06 (vs 0.30), 4 ALT reads (vs 19),
  80% with VAF < 0.1, C>T 59% (vs 28%) and T>C 23% (vs 17%), C>A unchanged (no oxidation
  signal), 26% confirmed by another caller (vs 73%). The same holds in sarek's original Mutect2
  VCF (not normalised, MNVs not split).
- Per ROQ bin the C>T excess persists up to 20-30 (72% at 0-5, 36% at 20-25) and is gone from 30
  (~20%, clean group >= 60: 30%); 30-60 have normal C>T but rising confirmation and VAF - real
  low-VAF calls (a call with few ALT reads cannot reach a high ROQ).
- Conclusion: mostly deamination-like orientation artefacts passing `FilterMutectCalls`; the
  consensus sets remove most of them (three quarters have no second caller). Affected:
  Mutect2-only counts, the `mutect2` signature set, and tumourevo (Mutect2 PASS only). See Open
  questions 2.

Literature: Costello et al. 2013, Nucleic Acids Res 41:e67 (8-oxoG C>A/G>T artefacts from
oxidation during shearing; orientation-based filtering); Chen et al. 2017, Science 355:752-756
(DNA damage is a pervasive cause of sequencing errors in public data incl. TCGA and 1000 Genomes,
dominating low-frequency variant calls; detected by strand/read imbalance); Do & Dobrovic 2015,
Clin Chem 61:64-71 (cytosine deamination -> C>T/G>A artefacts, FFPE and beyond); Benjamin et al.
2019, bioRxiv 10.1101/861054 (Mutect2 / FilterMutectCalls: F-score threshold, orientation-bias
model "crucial for good performance on FFPE samples"); GATK documentation for
LearnReadOrientationModel and FilterMutectCalls.

## Open questions

**1. Compare callers in single bases (MNVs split) or in events (doublets joined)?** *Open since 2026-10-07.*

The callers write the same doublet differently: SAGE and Mutect2 as one MNV record
(`162683852 GC>CT`), Strelka2 and MuSE as two SNVs (`162683852 G>C`, `162683853 C>T`). The
overlap matches on CHROM:POS:REF:ALT, so the formats have to be made the same first.

| | A: split MNVs into SNVs (current) | B: join adjacent SNVs into doublets |
|---|---|---|
| unit of overlap, set QC, 04 table | single bases | mutational events |
| SAGE count | 2,158 SNVs | 2,124 SNVs + 17 DBS |
| a caller that found only one half of a doublet | still matches that half | does not match |
| same unit as the signature analysis (doublets = 1 DBS) | no | yes |
| one number per caller everywhere | no (2,158 here, 2,124 in signatures) | yes |

To decide: which unit the presentation and paper should use for "SNVs per caller" and the
overlaps. If B: join Strelka2/MuSE adjacent SNVs into MNV records before matching, report DBS as
their own row in the 04 table, and rerun 01-04. Until decided, the scripts use A.

**2. Filter Mutect2's low-ROQ calls?** *Open since 2026-10-07; nothing applied yet.*

Proposed rule (from the ROQ bins, see Method decisions): remove Mutect2 PASS SNVs with
ROQ < 30 unless another caller confirms them. 30 removes the bins still carrying the C>T excess
(20-30) and keeps the 30-60 calls, which look real; "unless confirmed" keeps ~126 confirmed calls
from 20-30 that a plain cut would lose. 40 or 60 would be too harsh.
Once agreed: an *additional* tumourevo run (next to the existing one) on (a) the ROQ-filtered
Mutect2 calls and (b) the high-confidence consensus calls - PyClone-VI / MOBSTER got ~2,200
low-VAF artefacts (~40% of their input), which can fake a low-CCF subclone; and optionally the
same filter for the `mutect2` signature set. Methods wording: "Mutect2 PASS SNVs with ROQ < 30
were removed unless also called by another caller."

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
- `06_drivers.R` added; sarek VEP 116 annotation of all callers (`SAREK_STEP=annotate`, `../wgs/qsub_sarek_annotate.sh`).
- `05` fit quality: the FitMS cosine is labelled as the organ-signature fit (the RefSig plot is the same fit).
- `06` slide tables: `drivers_somatic_table` and `drivers_germline_table` (.csv/.md); germline includes SPRTN from the slice-BAM pileups, phasing and PURPLE CN.
- sarek annotate fixes: VEP settings via `-params-file` (Nextflow 26 passes `--x true` as a string); one samplesheet row per VCF with its own sample name (sarek requires a unique patient-sample-status-lane). Result: all 5 VCFs annotated, record counts unchanged (`--filter_common` left out).
- `01` section 5b: Mutect2 ROQ analysis - `roq_qc_plot()` (ROQ distribution, VAF low vs high ROQ with median/mean and Wilcoxon, ALT reads, ROQ vs VAF, classes, % confirmed) and `roq_cutoff_plot()` (evidence per ROQ bin, what each cut-off removes). See Method decisions and Open questions 2.
- Doublets counted once (see Method decisions): `vcf_sbs/` signature inputs, 01b runs the matrix generator twice (`matrix_generator/` for SBS96, `matrix_generator_all/` for DBS78 + ID83). Expected effect: SAGE SBS96 2,158 -> 2,124, identical to SIGS. Rerun 01-03, 01b, 01c, 05, 05b.
