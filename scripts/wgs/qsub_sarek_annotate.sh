#!/bin/bash
# ---------------------------------------------------------------------------
# PBS job: VEP annotation of every sarek caller's somatic VCF (Mutect2, Strelka2,
# MuSE - bcftools-normalised by 06 - and Manta SVs) - ./02_run_sarek.sh with
# SAREK_STEP=annotate. No realignment or calling; a few hours.
#
#   VEP_CACHE_VERSION=116 ../wgs_test/04_download_references.sh vep   # once, ~25 GB
#   cd scripts/wgs && qsub qsub_sarek_annotate.sh
#
# Writes to results/wgs/sarek/RJALS_annotate/annotation/<caller>/RJALS_Tm_vs_RJALS_N/.
# ---------------------------------------------------------------------------
#PBS -q q2
#PBS -l select=1:ncpus=16:mem=128gb
#PBS -l walltime=48:00:00
#PBS -N sarek_RJALS_annotate
#PBS -M volimpapat22@gmail.com
#PBS -m ae
#PBS -j oe

export JOB_MEMORY_GB=128        # keep equal to mem= above
export SAREK_STEP=annotate

cd "$PBS_O_WORKDIR" && [[ -f 00_config.sh ]] \
    || { echo "ERROR: submit from scripts/wgs/  (cd scripts/wgs && qsub qsub_sarek_annotate.sh)" >&2; exit 1; }
./02_run_sarek.sh
