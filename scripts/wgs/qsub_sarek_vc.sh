#!/bin/bash
# ---------------------------------------------------------------------------
# PBS job: extra sarek callers on the finished RJALS run, from its recalibrated
# CRAMs (no realignment) - ./02_run_sarek.sh with SAREK_STEP=variant_calling.
#
#   cd scripts/wgs && qsub qsub_sarek_vc.sh
#
# Writes to results/wgs/sarek/RJALS_vc; the main run's outdir is left alone.
# Change the callers in SAREK_TOOLS below - not with qsub -v, which this PBS
# refuses ("qsub: cannot send environment with the job").
# ---------------------------------------------------------------------------
#PBS -q q2
#PBS -l select=1:ncpus=40:mem=400gb
#PBS -l walltime=240:00:00
#PBS -N sarek_RJALS_vc
#PBS -M volimpapat22@gmail.com
#PBS -m ae
#PBS -j oe

export JOB_MEMORY_GB=400        # keep equal to mem= above
export SAREK_STEP=variant_calling
export SAREK_TOOLS=muse,msisensorpro

cd "$PBS_O_WORKDIR" && [[ -f 00_config.sh ]] \
    || { echo "ERROR: submit from scripts/wgs/  (cd scripts/wgs && qsub qsub_sarek_vc.sh)" >&2; exit 1; }
./02_run_sarek.sh
