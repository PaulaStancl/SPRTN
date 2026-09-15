# WGS processing on lobsang: sarek, oncoanalyser, tumourevo

These scripts set up and run three nf-core pipelines on one tumour/normal whole-genome pair:
- **sarek** calls somatic variants.
- **oncoanalyser** runs the Hartwig WiGiTS toolchain.
- **tumourevo** does clonal evolution analysis on sarek's output.

## Pipelines

| pipeline | version | input | output used downstream |
|---|---|---|---|
| nf-core/sarek | 3.10.0 | FASTQ | Mutect2 VCF + ASCAT copy number go to tumourevo |
| nf-core/oncoanalyser | 3.0.0 | FASTQ | stand-alone (PURPLE, LINX, ORANGE report) |
| nf-core/tumourevo | dev `738cb05` (no release yet) | sarek VCF + ASCAT | subclonal deconvolution and signatures |

All three use the existing Nextflow env `envs/nextflow-26.04.6`, the same one CHLOCK uses, and run their tools in apptainer containers.

## Test data: SEQC2 HCC1395 (breast cancer cell line)

| | sample | SRA run | depth | size |
|---|---|---|---|---|
| tumour | HCC1395T | SRR7890829 | ~53x | 94 GB |
| normal | HCC1395BL | SRR7890826 | ~55x | 97 GB |

- **Source:** public, from BioProject PRJNA489865, sequenced at Fudan on a HiSeq X, 2×150, one run per sample.
- **Subsampling:** script 06 reduces the reads to about 30x tumour and 20x normal.
- **Truth set:** SEQC2 v1.2.1 somatic SNVs and indels on GRCh38, downloaded to `data/wgs_test/HCC1395/truth_set/`.

## Layout on the server

```
/common/WORK/pstancl/
├── envs/nextflow-26.04.6/      reused
├── envs/wgs-tools/             seqtk, pigz, awscli, samtools, bcftools
├── singularity_cache/          shared container cache
├── PROGRAMI/nextflow/          NXF_HOME (pulled pipelines)
├── references/
│   ├── igenomes/Homo_sapiens/GATK/GRCh38/   sarek + tumourevo fasta
│   ├── hmf/oncoanalyser/                     GRCh38_hmf + WiGiTS resources
│   └── vep_cache/homo_sapiens/115_GRCh38/    tumourevo
└── projects/.../SPRTN/         this folder, copied to the server
    ├── scripts/wgs/            these scripts
    ├── data/wgs_test/HCC1395/  fastq_raw/  fastq_subsampled/  truth_set/
    ├── results/wgs/            sarek/  oncoanalyser/  tumourevo/  _install_tests/
    ├── work/wgs/               Nextflow work dirs (delete after each run)
    └── logs/wgs/
```

To change any path or version, edit `00_config.sh`.

## Order

Run from `scripts/wgs/`. Start steps 04 onwards inside `tmux`.

```bash
./01_check_system.sh                  # read-only: engine, CPU/RAM vs config, disk, network
./02_create_envs.sh                   # reuses nextflow-26.04.6, creates wgs-tools
./03_pull_pipelines.sh                # pulls the 3 pipelines at the pinned versions
./04_download_references.sh all       # ~106 GB: iGenomes, VEP 115, oncoanalyser refs
./05_download_test_data.sh            # 191 GB FASTQ + truth set, md5-checked
./06_subsample_fastq.sh               # -> ~30x T / ~20x N
./07_make_samplesheets.sh             # sarek + oncoanalyser sheets (--full = no subsampling)
./08_test_pipelines.sh all            # nf-core's own small tests: checks install, fills cache
./10_run_sarek.sh                     # 1-3 days
./11_run_oncoanalyser.sh              # 1-2 days, independent of sarek
./12_run_tumourevo.sh                 # after 10; builds its sheet from sarek output
```

- **Parallel work:** 04 and 05 are independent and can run at the same time in two tmux windows.
- **Run 10 and 11 one after the other.** Each one takes all 8 cores.
- **Resuming:** every run uses `-resume`, so after a crash re-run the same script.

## Before the first run

- **Resources:** check `conf/lobsang.config`. It is set to 8 CPUs and 400 GB RAM (lobsang: 8 CPUs, 7 TB RAM). `01_check_system.sh` compares that with the machine and fails if the config asks for more than the server has.
- **Time:** at 30x/20x on 8 cores, alignment alone takes most of a day. Mutect2 is the slowest step in sarek. For a faster first pass run `SAREK_TOOLS=strelka,manta,ascat ./10_run_sarek.sh`, but tumourevo needs the Mutect2 VCF.

## Things that are easy to trip over

- **tumourevo has no release.** It is pinned to a dev commit. Nextflow 26.04 parses strict syntax by default and tumourevo dev does not yet support it, so scripts 08 and 12 set `NXF_SYNTAX_PARSER=v1` for tumourevo only.
- **Two VEP versions.** sarek 3.10 ships VEP 116 and tumourevo ships VEP 115. Only tumourevo annotates here, so only the 115 cache is downloaded, and sarek runs without `vep`.
- **VCF sample names.** sarek writes them as `<patient>_<sample>`, e.g. `HCC1395_HCC1395T`. tumourevo's `tumour_sample` / `normal_sample` must match exactly, and script 12 checks this with `bcftools query -l`.
- **tumourevo can't use oncoanalyser output.** It accepts ASCAT, sequenza, Battenberg or facets for copy number, not PURPLE.
- **Two different GRCh38 builds.** oncoanalyser uses Hartwig's `GRCh38_masked_exclusions_alts_hlas`, not the GATK `Homo_sapiens_assembly38`. Both use `chr` names, but BAMs are not interchangeable between the two pipelines.
- **Single-sample signatures.** With one tumour, the signature tools in tumourevo only test that the pipeline runs.
- **Disk.** Delete `work/wgs/<run>` once a run's results are checked. Raw FASTQ (191 GB) can go after script 06.
