#!/bin/bash
# Chain tasks for the cycle_length simulation study (_scripts/simstudy.R):
# 240 fits x 8 chains = 1920 tasks, one chain each, pinned to one logical CPU
# of one physical core (2 CPUs = one core and its hyperthread sibling; the
# sibling is left idle, since siblings slow a chain).
#
# Node budget for THIS STUDY ONLY (Lucas, 2026-10-15): at most 75% of
# cbsugreischar, 192 CPUs and 773,505 MB.
#   M = floor(min(192 / 2, 773505 / 4096)) = min(96, 188) = 96
#
# Runs a frozen copy of the script, so _scripts/simstudy.R can be edited while
# tasks are still starting. Make the copy once, before submitting:
#   cp _scripts/simstudy.R _data/simstudy/simstudy-run.R
#   sbatch _scripts/simstudy.sh
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-1920%96
#SBATCH --job-name=simstudy
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test
#SBATCH --output=_data/simstudy/logs/%a.out
#SBATCH --error=_data/simstudy/logs/%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

set -u
CPU=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}' | cut -d, -f1 | cut -d- -f1)
echo "task ${SLURM_ARRAY_TASK_ID}: pinned to cpu ${CPU}"
SIMSTUDY_MODE=chain taskset -c "$CPU" \
    /programs/R-4.6.1/bin/Rscript --vanilla _data/simstudy/simstudy-run.R
