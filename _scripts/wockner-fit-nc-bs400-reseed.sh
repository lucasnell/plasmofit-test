#!/bin/bash -l

# Reseed-only retry of entry 39 (np_bs400_nc192), which did not converge.
#
# WHY ENTRY 39 AND NOT 41. Entry 41 (np_bs400_nc192_t) is the *same* config
# with adapt_delta 0.95 and max_treedepth 12 already baked in, so "a
# reseed-only run of entry 41" is a contradiction -- earlier notes said that
# and were wrong. Reseed-only means entry 39 with WOCKFIT_SEED changed and
# nothing else. _scripts/wockner-fit-nc-bs400-retry.sh runs 41-42 and stays
# unsubmitted.
#
# WHY RESEED FIRST. SLURM 29635 returned max R-hat 6.13 (n_c = 192) with ONE
# chain stuck far below the others: lp__ -929.0 against -832.3, -833.2 and
# -836.7. The majority chains agreed (cycle_length[1] 42.18, 42.28, 42.10,
# outlier 43.09). A single stuck chain is exactly what a different seed fixes,
# and reseeding changes one thing rather than three. The sampler already
# saturated max_treedepth in 46-61% of transitions, so raising adapt_delta to
# 0.95 and treedepth to 12 costs 2.5-5x: entry 39 took 8:20:35 here, so the
# tuned version would run 21-42 h and entry 42 would likely exceed its own
# 3-day walltime.
#
# WHAT IT DECIDES. Whether the n_c effect on cycle_length survives a pinned
# b_shape. The deterministic tests predict it should PERSIST AND GROW: 0.98 h
# at b_shape 400 against 0.72 h at 15, in all twelve cells of the
# least-squares table. If this converges and the effect is gone, the whole
# mechanistic account is wrong and should be revisited, not patched.
#
# IF A CHAIN STICKS AGAIN at a similar lp gap, the mode is real: report it as
# multimodality, show both modes, and move to entry 41 rather than reseeding
# a third time.
#
# WOCKFIT_SUFFIX IS NOT OPTIONAL. The output name is built from the config
# name alone -- the seed does not appear in it -- so without a suffix this
# would OVERWRITE _data/wock-fit-np_bs400_nc192.rds, the non-converged fit
# that is the evidence for why this rerun exists.
#
# Walltime 1 day against 8:20:35 observed, for headroom on a different seed.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-nc-bs400-reseed.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=39
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-ncbs-rs
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-ncbs-rs-%a.out
#SBATCH --error=wock-ncbs-rs-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_SEED=20261008
export WOCKFIT_SUFFIX=-seed2

Rscript --vanilla ../_scripts/wockner-fit.R
