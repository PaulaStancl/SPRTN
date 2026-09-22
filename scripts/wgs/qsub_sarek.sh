#!/bin/bash
# ---------------------------------------------------------------------------
# PBS job: nf-core/sarek on RJALS - runs ./02_run_sarek.sh inside one job.
#
#   cd scripts/wgs && qsub qsub_sarek.sh
#
# Nextflow uses the local executor inside this job, capped at the job's
# allocation: NCPUS comes from PBS, JOB_MEMORY_GB below must equal mem= above.
# Out of walltime? qsub again - the run resumes from the last finished task.
# Resubmit with the SAME size: a task's cpu count is part of its command line,
# so a different size re-runs every task whose cpus change (e.g. bwa-mem2).
# Needs SEX set in 00_config.sh and ./01_make_samplesheets.sh run beforehand.
# ---------------------------------------------------------------------------
#PBS -q q2
#PBS -l select=1:ncpus=40:mem=400gb
#PBS -l walltime=240:00:00
#PBS -N sarek_RJALS
#PBS -M volimpapat22@gmail.com
#PBS -m ae
#PBS -j oe

export JOB_MEMORY_GB=400        # keep equal to mem= above

cd "$PBS_O_WORKDIR" && [[ -f 00_config.sh ]] \
    || { echo "ERROR: submit from scripts/wgs/  (cd scripts/wgs && qsub qsub_sarek.sh)" >&2; exit 1; }
./02_run_sarek.sh
