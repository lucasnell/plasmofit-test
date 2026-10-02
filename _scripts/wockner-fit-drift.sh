#!/bin/bash -l

# Thread 6: how much do two BUILDS of the same source disagree?
#
# The b_shape regression turned up an excess that tracks the package build
# rather than the change: mu_logit_R at 3.5x the seed-only null, the
# (log10_total0, R) ridge at 2.8x, and b_shape itself at 1.0. That was never
# separated from the change, because the two builds compared differed by
# both. This separates it.
#
#   np_wide_total0          source e61fb47, built 2026-09-24, seed 538065874
#   np_wide_total0-drift    source e61fb47, built 2026-10-02, seed 538065874
#
# IDENTICAL SOURCE, identical seed, identical data. e61fb47 is dated
# 2026-09-24 08:52 and the reference fit was written at 14:50 the same day,
# so it is the source that produced it. Every difference this shows is the
# rebuild alone: a recompile reorders floating-point operations and HMC turns
# a last-bit difference into a different trajectory within a few leapfrog
# steps.
#
# Read it against the seed-only null (mean |z| 0.85, well calibrated) from
# _data/bshape-regression.log. If this reaches the ~1.45x the b_shape test
# showed, the excess is rebuild drift and the change is fully exonerated. If
# it sits near 1.0, the excess belongs to the change after all and the merge
# needs revisiting.
#
# WOCKFIT_LIB points at a build of e61fb47 installed to its own library, so
# nothing about the current install is disturbed; dependencies still resolve
# from the main library. Built by:
#   git worktree add --detach /home2/lan68/plasmofit/.prechange/src e61fb47
#   R CMD INSTALL --preclean -l /home2/lan68/plasmofit/.prechange/lib <that>
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-drift.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=9
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-fit-drift
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-fit-drift-%a.out
#SBATCH --error=wock-fit-drift-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_LIB=/home2/lan68/plasmofit/.prechange/lib
export WOCKFIT_SUFFIX=-drift

Rscript --vanilla ../_scripts/wockner-fit.R
