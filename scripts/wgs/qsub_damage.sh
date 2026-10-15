#!/bin/bash
# ---------------------------------------------------------------------------
# PBS job: 13_damage_metrics.sh - Picard artefact metrics + GIV on tumour and normal (a chr20
# slice by default), then read positions at 07's sites. Single-threaded Picard and Python;
# the cpus only speed up the CRAM slicing.
#
#   cd scripts/wgs && qsub qsub_damage.sh
#   qsub -v DAMAGE_REGIONS=chr2,chr20 qsub_damage.sh
#   qsub -v STEPS=reads qsub_damage.sh         # only read positions (after 07 wrote the sites)
# ---------------------------------------------------------------------------
#PBS -q q2
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=48:00:00
#PBS -N damage_RJALS
#PBS -M volimpapat22@gmail.com
#PBS -m ae
#PBS -j oe

cd "$PBS_O_WORKDIR" && [[ -f 00_config.sh ]] \
    || { echo "ERROR: submit from scripts/wgs/  (cd scripts/wgs && qsub qsub_damage.sh)" >&2; exit 1; }
./13_damage_metrics.sh
