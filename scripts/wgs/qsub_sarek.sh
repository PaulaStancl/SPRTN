#!/bin/bash
# ---------------------------------------------------------------------------
# PBS job: nf-core/sarek on RJALS - runs ./02_run_sarek.sh inside one job.
#
#   cd scripts/wgs && qsub qsub_sarek.sh
#
# Nextflow uses the local executor inside this job, capped at the job's
# allocation: NCPUS comes from PBS, JOB_MEMORY_GB below must equal mem= above.
# Out of walltime? qsub again - the run resumes from the last finished task.
# Resubmit with the SAME cpu count: a task's cpu count is part of its command
# line, so changing it re-runs every task whose cpus change (e.g. bwa-mem2).
# Memory can be changed freely - no task here asks for more than 30 GB.
#
# Why 1000gb for a run whose biggest task asks 30 GB: the 2026-09-23 run was
# killed (SIGTERM, exit 143) after 23 h with PBS reporting 400 GB used against a
# 400 GB request, while MarkDuplicates was writing ~200 GB BAMs for both samples
# at once. Page cache from those writes counts towards the job's memory cgroup.
# Needs SEX set in 00_config.sh and ./01_make_samplesheets.sh run beforehand.
# ---------------------------------------------------------------------------
#PBS -q q2
#PBS -l select=1:ncpus=40:mem=1000gb
#PBS -l walltime=240:00:00
#PBS -N sarek_RJALS
#PBS -M volimpapat22@gmail.com
#PBS -m ae
#PBS -j oe

export JOB_MEMORY_GB=1000       # keep equal to mem= above

cd "$PBS_O_WORKDIR" && [[ -f 00_config.sh ]] \
    || { echo "ERROR: submit from scripts/wgs/  (cd scripts/wgs && qsub qsub_sarek.sh)" >&2; exit 1; }
./02_run_sarek.sh
