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

1. **`SEX` is set to `XY`** (male) in `00_config.sh`. sarek needs it for ASCAT, and it is part of every task's inputs, so changing it after sarek has started restarts the whole run. oncoanalyser does not need it.
2. **`CANCER_TYPE` is set to `HCC`** (hepatocellular carcinoma). Only tumourevo uses it, for driver annotation. An IntOGen code that isn't in tumourevo's driver table doesn't fail: it silently falls back to pan-cancer drivers, so check any new code against `Compendium_Cancer_Genes.tsv`.
3. **Copy the oncoanalyser reference config** once: `cp ../wgs_test/conf/oncoanalyser_refdata.config conf/`. It is gitignored because `04_download_references.sh` writes it on the server.
4. **Run on a machine that mounts `/common/RAW`.** Every script checks this first and stops if it can't read the FASTQ.

## Order

Run from `scripts/wgs/`. Either submit the steps as PBS jobs (next section) or run them directly inside `screen` or `tmux`.

```bash
./01_make_samplesheets.sh     # sarek + oncoanalyser sheets (sarek's only once SEX is set)
./02_run_sarek.sh             # FASTQ -> Mutect2, Strelka, Manta, ASCAT
./03_run_oncoanalyser.sh      # FASTQ -> Hartwig WiGiTS + ORANGE report; independent of sarek
./04_run_tumourevo.sh         # after 02; builds its sheet from sarek output
```

- **Run directly on lobsang, 02 and 03 must go one after the other,** because each one takes all 8 cores. As PBS jobs they can run at the same time. Either order works; if the patient's sex isn't confirmed yet, start with 03.
- **Resuming:** re-run the same script after a crash and it continues from the last finished task.

## Running as PBS jobs (queue q2)

There is one job per pipeline. Each job runs Nextflow with the local executor on one node, capped at that job's allocation. PBS sets `NCPUS`, and each `qsub_*.sh` sets `JOB_MEMORY_GB`, which must equal its `mem=`. `nf_run` writes these caps to `work/wgs/<run>/resources.config`, leaving 12 GB for Nextflow itself.

| job | select | walltime | notes |
|---|---|---|---|
| `qsub_sarek.sh` | `ncpus=24:mem=250gb` | 240 h | bwa-mem2 (which asks for exactly 24 cpus) and Mutect2 dominate. |
| `qsub_oncoanalyser.sh` | `ncpus=24:mem=250gb` | 240 h | bwa-mem2, REDUX, SAGE, ESVEE (asks for 32 cpus, capped to 24) |
| `qsub_tumourevo.sh` | `ncpus=8:mem=48gb` | 48 h | VEP plus clonality; mostly single-threaded, so more cpus don't help |
| `qsub_check.sh` | `ncpus=1:mem=2gb` | 15 min | pre-flight: node, singularity, mounts, internet |

Memory needs at least about 8 GB per cpu, because more cpus means more tasks running at once; 250 GB for 24 cpus leaves headroom. The whole job must fit on **one** node, so check the node sizes first (`pbsnodes -a | grep -E 'resources_available.(ncpus|mem) ='`). A request bigger than any node either gets rejected or sits in the queue forever. To change the size, edit `ncpus=`/`mem=` and set `JOB_MEMORY_GB` to the same value as `mem=`. `qsub` rejects a request above the queue's limits straight away, so a walltime that is too long fails at submission rather than days later.

```bash
cd scripts/wgs
./01_make_samplesheets.sh                          # on the login node, once
qsub qsub_check.sh                                 # read check_q2.o<id>: every line OK?
qstat -Qf q2 | grep -E 'resources_(max|default)'   # queue limits: walltime, ncpus, mem
SAREK=$(qsub qsub_sarek.sh)
qsub qsub_oncoanalyser.sh
qsub -W depend=afterok:$SAREK qsub_tumourevo.sh    # starts when sarek finishes OK
```

- **Out of walltime:** `qsub` the same script again. The run resumes from the last finished task, because every run has its own launch directory.
- **Job output:** PBS writes `<name>.o<jobid>` into `scripts/wgs/` (gitignored). The full Nextflow log is in `logs/wgs/`.

## Time and disk

- **Time:** roughly 2,000–3,000 CPU-hours each for sarek and oncoanalyser at this depth. That is about 1–2 weeks on 8 cores and roughly 4–6 days on 24. Steps that don't parallelise (markdup, ASCAT, some Manta stages) keep it from scaling perfectly. These are estimates scaled from the test depth; no full-depth run has been timed yet. Alignment and Mutect2 dominate.
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
