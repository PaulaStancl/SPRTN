# Human tumour → expected mouse alterations (H2M)

These scripts turn the RJALS WGS results (nf-core/oncoanalyser 3.0.0, `--mode wgts`, GRCh38_hmf, DNA only) into an **expectation map**. For every human somatic alteration they predict what the equivalent event in mouse would be, and where it would sit in GRCm39 and GRCm38. The predictions are frozen before any mouse data are looked at. Later, the same code base compares them with the three GEMM tumours (nf-core/sarek). The human reference frame is HCC and hepatoblastoma.

## Principles

- **Gene-centred.** Everything is keyed on Ensembl gene IDs and orthologs. Raw coordinates are never compared across species; liftOver is used only to place non-coding variants and breakends and to cross-check the CN map.
- **No counterpart is its own category.** A human event can have no possible mouse counterpart because there is no ortholog, the gene is absent from the build, the variant cannot be modelled, or the region does not lift over. That is recorded as such, never as a discordance.
- **`neutral` is not `not_assessable`.** `neutral` means assessed with no event; `not_assessable` means it could not be assessed. The distinction is kept everywhere.
- **Fusions are a feasibility flag, never a prediction.**
- **Predictions are frozen** (`07b_freeze.R`) before the mouse data. The mouse modules read only the frozen copy, and refuse to run if its checksums differ.

## Setup (server)

```bash
cd /common/WORK/pstancl/projects/MariaBoskovic/SPRTN/scripts/cross_species/h2m
micromamba create -y -p /common/WORK/pstancl/envs/cross-species-h2m -f envs/environment.yml
micromamba activate /common/WORK/pstancl/envs/cross-species-h2m
Rscript tests/run_tests.R            # unit tests (synthetic inputs)
Rscript scripts/00_resources.R --dry-run   # what would be downloaded, and how much
```

Everything is R (data.table) except the H2M call itself. H2M (`bioh2m`) exists only as a Python package, so `01` writes H2M's input table and runs `scripts/h2m_model.py` with the interpreter set in `config.yaml` (`python:`).

Paths, sample IDs, builds, thresholds and releases are all in `config/config.yaml`. Outputs go to `results/wgs/cross_species/h2m/ens<release>/` on the server. Downloads go to `/common/WORK/pstancl/references/cross_species/`. Neither is in git, and the scripts refuse cloud-synced output folders.

## Run

```bash
./run_all.sh                        # 00 -> 07b with Ensembl 110
./run_all.sh --release latest       # the same chain on the current Ensembl, separate ens<N>/ folder
FROM=03 ./run_all.sh                # restart at a later step
FREEZE=0 ./run_all.sh               # everything except the freeze
```

Each script also runs on its own, for example `Rscript scripts/03_cn_to_mouse.R [--config ...] [--release ...]`. Each one writes a log to `logs/`, its rows to `human_predictions/summary/<layer>.tsv`, and its inputs, resource checksums and software versions to `provenance.json`. Patient files are recorded by path and size only.

## Scripts, inputs and outputs

All outputs are under `<out_dir>/human_predictions/`. Column lists are in each script's header.

| step | reads | writes |
|---|---|---|
| `00_resources.R` | Ensembl FTP (GTFs, FASTAs, Compara homologies), UCSC, Reactome | in `resources_dir`: ortholog table (gene IDs, symbols, homology type, % identity, confidence, coordinates in both species), human GTF (release R), mouse GTFs (GRCm39 R, GRCm38 102), mouse FASTAs, chains hg38→mm39/mm10 and mm39→mm10, Reactome; `ortholog_release_diff_<R>_vs_<R2>.tsv` |
| `01_snv_coding_h2m.R` (+ `h2m_model.py`) | `purple/*.purple.somatic.vcf.gz`, `*.purple.driver.catalog.somatic.tsv` | `expected_mouse_snv.tsv`, `h2m_raw.tsv`, `h2m_failed.tsv`, `h2m_validation.tsv` |
| `02_snv_noncoding_liftover.R` | PURPLE somatic VCF | `expected_mouse_noncoding.tsv`, `expected_mouse_noncoding_multihits.tsv` |
| `03_cn_to_mouse.R` | `*.purple.cnv.gene.tsv`, `*.purple.cnv.somatic.tsv`, `*.purple.purity.tsv` | `expected_mouse_cn.<build>.bed`, `expected_mouse_cn_genes.tsv`, `cn_one2many.tsv`, `cn_segment_mapping.tsv`, `cn_liftover_crosscheck.tsv` |
| `04_sv_to_mouse.R` | `linx/*.linx.{vis_sv_data,clusters,breakend,fusion,driver.catalog}.tsv` | `expected_mouse_sv.<build>.bedpe`, `expected_mouse_sv_genes.tsv`, `expected_mouse_fusions.tsv` |
| `05_signatures_expected.R` | PURPLE somatic VCF, `sigs/*.sig.allocation.tsv` | `expected_mouse_sbs96.tsv`, `expected_signatures.tsv`, `trinucleotide_frequencies.tsv`, `sbs96_crosscheck.tsv` |
| `05b_pathways.R` | outputs of 01-04, COSMIC CGC (yours), Reactome | `gene_pathways.tsv`, `pathways_affected.tsv` |
| `06_build_matrix.R` | outputs of 01-04 | `gene_matrix.tsv` (gene × human / expected_mouse / mouse_1..3 / freq_HCC / freq_HB) |
| `06b_build_event_table.R` | outputs of 01-05b | `event_table.tsv` (one row per human alteration or mouse fragment, mouse_1..3 = `pending`) |
| `07_report.R` | everything above | `report.html` |
| `07b_freeze.R` | everything above | `FROZEN/` (copies, `MANIFEST.sha256`, `FROZEN_AT.txt`, read-only); `--verify` checks it |

### Categories used throughout

- **Counterpart (summary tables, report):**
  - `yes` means a predicted mouse event exists.
  - `uncertain` covers 1:many orthologs, multi-mapped regions, H2M class 2 (effect may differ), and LOH (inbred mice have no heterozygous SNPs, so mouse LOH is `not_assessable` unless heterozygous regions are found).
  - `no` means no possible counterpart.
- **Mapping confidence (event table):**
  - `high` and `partial` are self-explanatory.
  - `ambiguous` covers 1:many orthologs and multi-mapped regions.
  - `none` means no counterpart.
  - For CN segments, confidence comes from the fraction of the segment's genes placed in mouse and the fraction of its length that lifts over (thresholds `cn.conf_*`).
- **H2M classes:**
  - 0: modelled with the same DNA change.
  - 1: modelled by an alternative change giving the same protein change.
  - 2: modelled, but the effect may differ.
  - 3–6: not modelable.
- **CN states (ploidy-relative, `cn.*`):**
  - `DEEP_DEL` below 0.5 copies.
  - `LOSS` below 0.75 × ploidy.
  - `GAIN` at 1.25 × ploidy or more.
  - `AMP` at 3 × ploidy or more (Hartwig's AMP).
  - Losses are called from the gene's minimum CN and gains from its maximum. `partial` means the two disagree.

## Decisions and deviations from the original spec

- **Plain scripts, not Snakemake** (2026-10-28, Paula). R wherever possible, Python only for H2M.
- **Ensembl 110** is the pinned release, because H2M 1.0.3's bundled tables were built from GENCODE v44 / vM33, which is Ensembl 110. `--release latest` reruns the whole chain on the current release into its own folder, and `00` writes a per-gene ortholog diff between the two releases.
- **H2M** is used only for its modelling engine (`model_batch`, `ver=38`):
  - Its mouse side is **GRCm39 only**; GRCm38 coordinates come from an mm39→mm10 liftOver.
  - Orthologs come from the Ensembl Compara homology files of the release (not BioMart: the archive BioMart returned HTTP 500), not H2M's own name table.
  - Transcripts are PAVE's transcript, otherwise the Ensembl canonical one.
  - Its annotation DBs are built from only the transcripts in use, and its genome is read lazily from the iGenomes GRCh38 FASTA (same coordinates as GRCh38_hmf, no masked regions). That avoids H2M's Dropbox downloads and whole genomes held in memory.
  - Five known hotspots (KRAS G12D, TP53 R175H→Trp53 R172H, CTNNB1 S45F, PIK3CA H1047R, BRAF V600E→Braf V637E) are run first. Any failure stops `01` unless `--skip-validation` is given.
- **SBS96** is counted directly from the VCF and FASTA (Rsamtools / Biostrings) instead of a new SigProfilerMatrixGenerator install. It is compared with the existing SigProfiler matrix from `analysis/01b` (`sbs96_reference`). MNV records are left out, as in `analysis/01b`.
- **Synteny blocks** come from the ortholog order: runs of 1:1 orthologs that keep their human order on one human chromosome. They are not downloaded from Ensembl.
- **SV breakends** come from `linx.vis_sv_data.tsv`, the same columns `analysis/run_chromothripsis.R` uses. Chromothripsis regions (e.g. chr9p from that analysis) go in `sv.chromothripsis_regions`.
- **The COSMIC Cancer Gene Census** needs a COSMIC login, so you download it yourself and set `cosmic_cgc` / `cosmic_cgc_version`. Without it, missense and promoter directions stay `unknown`.

## Mouse modules (10-12, not written yet)

Agreed design, to be implemented once the mouse data details are known:

1. **`10_mouse_ingest.R`** reads sarek outputs per mouse: Strelka2 / Mutect2 VCFs, Manta VCF, CNVkit `.cns` / `.call.cns`, whichever exist. A missing layer is `not_assessable`, never `absent`. A configurable BED of engineered loci labels events as engineered / secondary / other.
2. **`11_compare.R`**:
   - First, the mouse analysis is run genome-wide, independent of the predictions. Only then are the frozen predictions tested.
   - Outcome per prediction and mouse is one of: `equivalent_observed`, `related_gene_or_pathway_observed` (via the Reactome R-MMU IDs from `05b`), `absent_adequately_assessed`, or `not_assessable`.
   - Directional CN concordance with block weighting, Jaccard, and a permutation test by circular shifts of mouse segments (regioneR-style).
   - Mouse-only secondary events (not engineered, not predicted) are reported separately when they recur in 2 or 3 of the 3 mice.
3. **`12_callability.R`** gives, per predicted locus, mouse coverage and mappability from the CRAMs or mosdepth output, which decides `absent_adequately_assessed` vs `not_assessable`.

All three call `verify_frozen()` (`common.R`) first and read predictions only from `FROZEN/`.

**Open, to confirm before writing them:**
- sarek version and `--tools`.
- Mouse genome build(s).
- The engineered allele(s) and their coordinates.
- Whether any mouse has a matched normal of a different strain.
