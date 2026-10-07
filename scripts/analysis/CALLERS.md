# Variant callers: filters and fields

How each caller's output is read by the analysis scripts (`00_setup.R`, `01`-`04`).

## MuSE (sarek, extra run `results/wgs/sarek/RJALS_vc`)

Somatic SNV caller (no indels). It runs in two steps:

1. `MuSE call` pre-filters candidate positions and computes position-specific summary
   statistics with a Markov substitution model, from the tumour and normal BAM/CRAM.
2. `MuSE sump` computes tier-based cut-offs from a **sample-specific error model** and writes
   the VCF. It has one mode for WGS (`-G`) and one for WES (`-E`), and uses dbSNP (`-D`).

The VCF header records which mode ran (`##MuSE_sump=...`, look for `-G`).

**This project (RJALS_vc run):** `sump -G -I RJALS_Tm_vs_RJALS_N.MuSE.txt -n 12 -D dbsnp_146.hg38.vcf.gz`.
That's WGS mode with dbSNP 146 (hg38), taken from the VCF header. So the WGS recommendation
below applies, and all tiers are kept.

### FILTER: confidence tiers, not pass / fail

| FILTER | MuSE's description | meaning |
|---|---|---|
| `PASS` | "Accept as a confident somatic mutation" | highest confidence |
| `Tier1` | "Confident level 1" | |
| `Tier2` | "Confident level 2" | |
| `Tier3` | "Confident level 3" | decreasing confidence |
| `Tier4` | "Confident level 4" | |
| `Tier5` | "Confident level 5" | lowest confidence |

MuSE's authors recommend: *"For WGS data, we recommend to include the calls of all categories
from MuSE for downstream analysis. For WES data, we recommended to include the calls of all
categories except Tier 5"* (README, below).

**In this analysis** a MuSE record counts as a call if its FILTER is `PASS` or `Tier1`-`Tier5`.
That's `is_pass()` in `00_setup.R`. To keep fewer tiers, set e.g. `SPRTN_MUSE_TIERS=1,2,3,4`
before running the R scripts. `qc_vcf/raw/muse_tiers.csv` lists how many records are in each
tier and whether they were kept.

For every other caller (Mutect2, Strelka2, SAGE), only `FILTER = PASS` counts.

### FORMAT fields

The sample columns are named `TUMOR` and `NORMAL`.

| field | MuSE's description | notes |
|---|---|---|
| `GT` | Genotype | `0/1` in the tumour, `0/0` in the normal for a somatic call |
| `DP` | Read depth at this position in the sample | |
| `AD` | Depth of reads supporting alleles 0/1/2/3... | REF first, then ALT: `80,4` = 80 REF reads, 4 ALT reads |
| `BQ` | Average base quality for reads supporting alleles | same order as `AD`; `0` where there are no reads for that allele |
| `SS` | Variant status relative to non-adjacent Normal: 0 = wildtype, 1 = germline, 2 = somatic, 3 = LOH, 4 = post-transcriptional modification, 5 = unknown | `2` in the tumour column for somatic calls; `.` in the normal |

INFO is always `SOMATIC`. QUAL is `.`, since MuSE gives no per-call score.

Example (made-up numbers): `GT:DP:AD:BQ:SS  0/1:84:80,4:31,32:2` is a heterozygous call with
84 reads (80 REF, 4 ALT), average base quality 31 (REF) and 32 (ALT), status somatic.

**In this analysis** MuSE writes no allele fraction, so the tumour VAF is computed from `AD`:
ALT / (REF + ALT) (`read_vcf_table()` in `00_setup.R`). The normal VAF is computed the same way,
and `T_ALT` is the second `AD` value.

### Sources

- MuSE README (filter tiers, WGS/WES recommendation, `call` / `sump`):
  https://github.com/wwylab/MuSE/blob/master/README.md
- FILTER and FORMAT definitions, from the header of MuSE's example VCF:
  https://github.com/wwylab/MuSE/blob/master/example/COLO829_illumina_1_100000_1000000.vcf
- MuSE 1: Fan Y, et al. *MuSE: accounting for tumor heterogeneity using a sample-specific error
  model improves sensitivity and specificity in mutation calling from sequencing data.*
  Genome Biol 2016;17:178. doi:10.1186/s13059-016-1029-6
- MuSE 2: Ji S, Zhu T, Sethia A, Wang W. *Accelerated somatic mutation calling for whole-genome
  and whole-exome sequencing data from heterogenous tumor samples.* Genome Res 2024.
  doi:10.1101/gr.278456.123
