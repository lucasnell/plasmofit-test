#!/bin/bash -l

# Schedule-bias simulation. Submit from the repo root, on the cluster, with
# _data/ populated:
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim.sh
#
# Each task is one (arm, replicate) pair and takes ~2 h, the same cost as a
# real fit. See _scripts/wockner-schedule-sim.R for what it does and why.
#
# TASK INDEX = row of expand_grid(arm = names(ARMS), rep = seq_along(REP_SEEDS)),
# so REP VARIES FASTEST and adding a seed SHIFTS every index after `default`.
# With 7 arms x 5 replicates:
#    1-5  default        6-10  wide          11-15 tight_sigma (cancelled)
#   16-20 no_hier       21-25  wide_bshape   26-30 wide_total0
#   31-35 wide_nuis
# Do not read these off by hand; regenerate the table and check it against
# what is on disk, as claude/handoff-2026-09-25.md shows.
#
# The array below is the 13 budget-arm tasks that were NOT on disk on
# 2026-09-25: replicates 4-5 of every arm, plus no_hier-rep1, which was never
# run (findings.md's hierarchy table only ever had rep2 and rep3).
# `tight_sigma` is excluded -- cancelled, see findings.md.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=4,5,9,10,16,19,20,24,25,29,30,34,35
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-schedsim
#SBATCH --output=_data/wock-schedsim-%a.out
#SBATCH --error=_data/wock-schedsim-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/wockner-schedule-sim.R
