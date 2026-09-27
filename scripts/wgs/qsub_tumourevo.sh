#!/bin/bash
# ---------------------------------------------------------------------------
# PBS job: nf-core/tumourevo on RJALS - runs ./04_run_tumourevo.sh inside one
# job. Needs sarek's Mutect2 + ASCAT output, so start it after sarek - or let
# PBS wait for it:
#
#   cd scripts/wgs
#   SAREK=$(qsub qsub_sarek.sh)
#   qsub -W depend=afterok:$SAREK qsub_tumourevo.sh
#
# Light compared with sarek/oncoanalyser (VEP + clonal deconvolution), but sized
# for what its tasks ask: process_high wants 84 GB and process_high_memory 200 GB,
# so 256gb leaves 244 GB for tasks; 12 cpus matches process_high's request.
# tumourevo's per-task time limits are tighter than sarek's (2 h by default,
# 6/8/10 h per label) and Nextflow kills a local task that passes them - see the
# markdup failures in qsub_sarek.sh - so 04 lifts every task to 240 h.
# JOB_MEMORY_GB below must equal mem= above. Needs CANCER_TYPE in 00_config.sh.
# ---------------------------------------------------------------------------
#PBS -q q2
#PBS -l select=1:ncpus=12:mem=256gb
#PBS -l walltime=240:00:00
#PBS -N tevo_RJALS
#PBS -M volimpapat22@gmail.com
#PBS -m ae
#PBS -j oe

export JOB_MEMORY_GB=256        # keep equal to mem= above

cd "$PBS_O_WORKDIR" && [[ -f 00_config.sh ]] \
    || { echo "ERROR: submit from scripts/wgs/  (cd scripts/wgs && qsub qsub_tumourevo.sh)" >&2; exit 1; }
./04_run_tumourevo.sh
