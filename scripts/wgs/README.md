# WGS on the real data: RJALS tumour/normal

These scripts run nf-core/sarek, oncoanalyser and tumourevo on the SPRTN patient pair. One-time setup (envs, pipeline pulls, references) and the HCC1395 test run live in `../wgs_test/`, and nothing here downloads anything.

## Data

| | sample | FASTQ | lane-pairs |
|---|---|---|---|
| tumour | RJALS_Tm | ~166 GB | 3 |
| normal | RJALS_N | ~164 GB | 3 |

- **Location:** `/common/RAW/pstancl/MariaBoskovic/SPRTN/wgs/X208SC25056159-Z01-F001`, the Novogene delivery. It passed `md5sum -c MD5.txt` (15/15) on arrival from the drive. The check after the move to `/common/RAW` is logged in `md5check_RAW.log` in that folder.
- **Input:** only `01.RawData/` is used. `02.Bam/` holds Novogene's own bwa BAMs (no GATK, no markdup), and both pipelines realign from FASTQ instead.
- **Lanes:** each sample has 3 lane-pairs on 2 flowcells: `22VTTHLT4` L7 and L8, and `22VTVTLT4` L6. `01_make_samplesheets.sh` finds them from the filenames and writes one row per lane-pair, with lane `<flowcell>_L<n>`. That keeps read groups separate for sarek's BQSR.
- **Depth:** roughly 90x per sample. This is an estimate from file size compared with the HCC1395 test data (94 GB ≈ 53x).

## Before the first run

1. **Set `SEX` in `00_config.sh`** to `XX` or `XY`. sarek needs it for ASCAT, and it is part of every task's inputs, so changing it after sarek has started restarts the whole run. oncoanalyser does not need it.
2. **Set `CANCER_TYPE`** (the IntOGen code) before `04`. Only tumourevo uses it.
3. **Copy the oncoanalyser reference config** once: `cp ../wgs_test/conf/oncoanalyser_refdata.config conf/`. It is gitignored because `04_download_references.sh` writes it on the server.
4. **Run on a machine that mounts `/common/RAW`.** Every script checks this first and stops if it can't read the FASTQ.

## Order

Run from `scripts/wgs/`, inside `screen` or `tmux`.

```bash
./01_make_samplesheets.sh     # sarek + oncoanalyser sheets (sarek's only once SEX is set)
./02_run_sarek.sh             # FASTQ -> Mutect2, Strelka, Manta, ASCAT
./03_run_oncoanalyser.sh      # FASTQ -> Hartwig WiGiTS + ORANGE report; independent of sarek
./04_run_tumourevo.sh         # after 02; builds its sheet from sarek output
```

- **Run 02 and 03 one after the other.** Each one takes all 8 cores. Either order works; if the patient's sex isn't confirmed yet, start with 03.
- **Resuming:** re-run the same script after a crash and it continues from the last finished task.

## Time and disk

- **Time:** expect roughly 1–2 weeks each for sarek and oncoanalyser on 8 cores. This is scaled from about 3.5× the test depth, and no full-depth run has been timed yet. Alignment and Mutect2 dominate.
- **Disk:** each work dir can reach several TB, and `/common/WORK` had 22 TB free. Delete `work/wgs/<run>` once that run's results are checked.

## Layout on the server

```
/common/RAW/pstancl/MariaBoskovic/SPRTN/wgs/X208SC25056159-Z01-F001/   input, read-only to the pipelines
/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/
├── scripts/wgs/           these scripts
├── scripts/wgs_test/      setup + HCC1395 test run
├── results/wgs/           sarek/RJALS  oncoanalyser/RJALS  tumourevo/RJALS
├── work/wgs/<run>/        .nextflow/ history + work/   (delete after each run)
└── logs/wgs/
```

`results/`, `work/` and `logs/` are shared with the test run, and the dataset name (`RJALS` vs `HCC1395`) keeps them apart.

## Differences from `../wgs_test`

- **Each run has its own launch directory** (`work/wgs/<run>/`). A plain `-resume` resumes whichever run was started *last* in the launch directory. With a shared directory, running oncoanalyser between two sarek attempts would make sarek start over.
- **`/common/RAW` is bound read-only** into the containers (`conf/lobsang.config`), so no pipeline step can write into the archive.
- **`umask 077`:** everything these runs write is readable only by you.
- The pipeline versions, references, container engine and resources are the same as in the test run. See `../wgs_test/README.md` for why `singularity` is used rather than `apptainer`, and for the tumourevo and VEP pins.
