#!/bin/bash -l

# Reseed-only retry of entry 40 (np_bs400_nc384), the second and last rung of
# the production n_c ladder with b_shape pinned at 400.
#
# WHY THIS IS THE ONLY THING BLOCKING A REPORTABLE CYCLE LENGTH. Entry 39
# (n_c = 192) converged on a reseed -- SLURM 30525, max R-hat 1.0432, mean
# cycle_length 41.203 h over 13 trials -- so ONE rung of the ladder exists and
# the other has never converged. 29635_40 returned max R-hat 8.34 with one
# chain stuck at lp__ -996.0 against -925.6, -955.3 and -957.8, and its
# fitted values and loo_compare are unusable. The question the ladder exists
# to answer -- does the n_c effect on cycle_length survive a pinned b_shape?
# -- needs both rungs and cannot be answered from one.
#
# WHY RESEED RATHER THAN TUNE, AGAIN. It worked on entry 39, for the same
# failure mode (one stuck chain), and it changes one thing. The tuned retry
# (_scripts/wockner-fit-nc-bs400-retry.sh, entries 41-42) is now LESS
# justified than when it was written: it raises max_treedepth on the strength
# of 46-61% saturation recorded from 29635, but 30525 saturated treedepth in
# only 2% of transitions with a maximum of 10. Keep it unsubmitted.
#
# WOCKFIT_SUFFIX IS NOT OPTIONAL. Output names are built from the config name
# alone and the seed never appears in them, so without a suffix this would
# OVERWRITE _data/wock-fit-np_bs400_nc384.rds -- the non-converged fit that is
# the evidence for why this rerun exists. See claude/gotchas.md.
#
# WHAT TO CHECK WHEN IT LANDS. Max R-hat < 1.05 is the gate, but check lp__
# per chain first: one stuck chain and a general failure want different fixes.
# Then compare cycle_length against entry 39's 41.203 h. The deterministic
# tests predict the n_c effect PERSISTS AND GROWS with b_shape pinned -- 0.98 h
# at b_shape 400 against 0.72 h at 15, in all twelve cells of the least-squares
# table. If it converges and the effect is gone, the mechanistic account is
# wrong and should be revisited rather than patched. Note also that the ML
# profile, with b_shape free and no priors, puts cycle_length at 40.93-41.15 h
# once n_c >= 384.
#
# IF A CHAIN STICKS AGAIN at a similar lp gap, the mode is real: report it as
# multimodality, show both modes, and do not reseed a third time.
#
# COST. 29635_40 ran 24:34:15 and completed; it failed the convergence gate,
# not the walltime. Walltime here is 3 days, the measurement times three,
# because every cost extrapolation in this project has come in low.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-nc-bs400-reseed40.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=40
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=wock-ncbs-rs40
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-ncbs-rs40-%a.out
#SBATCH --error=wock-ncbs-rs40-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_SEED=20261009
export WOCKFIT_SUFFIX=-seed2

Rscript --vanilla ../_scripts/wockner-fit.R
