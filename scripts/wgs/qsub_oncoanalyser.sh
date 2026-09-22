#!/bin/bash
# ---------------------------------------------------------------------------
# PBS job: nf-core/oncoanalyser on RJALS - runs ./03_run_oncoanalyser.sh
# inside one job. Independent of sarek, so both jobs can run at the same time.
#
#   cd scripts/wgs && qsub qsub_oncoanalyser.sh
#
# Nextflow uses the local executor inside this job, capped at the job's
# allocation: NCPUS comes from PBS, JOB_MEMORY_GB below must equal mem= above.
# Out of walltime? qsub again - the run resumes from the last finished task.
# Needs ./01_make_samplesheets.sh run beforehand (SEX not required).
# ---------------------------------------------------------------------------
#PBS -q q2
#PBS -l select=1:ncpus=40:mem=400gb
#PBS -l walltime=240:00:00
#PBS -N onco_RJALS
#PBS -M volimpapat22@gmail.com
#PBS -m ae
#PBS -j oe

export JOB_MEMORY_GB=400        # keep equal to mem= above

cd "$PBS_O_WORKDIR" && [[ -f 00_config.sh ]] \
    || { echo "ERROR: submit from scripts/wgs/  (cd scripts/wgs && qsub qsub_oncoanalyser.sh)" >&2; exit 1; }
./03_run_oncoanalyser.sh
