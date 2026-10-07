#!/bin/bash -l

# Retry of entries 39-40, which did not converge.
#
# SLURM 29635 returned max R-hat 6.13 (n_c = 192) and 8.34 (384), with ONE
# chain stuck far below the others in each: lp__ -929.0 against -832.3,
# -833.2, -836.7 at 192, and -996.0 against -925.6, -955.3, -957.8 at 384.
# The majority chains agreed (cycle_length[1] = 42.18, 42.28, 42.10, outlier
# 43.09), so the information is probably there -- but a run with a stuck chain
# is not a posterior, and reporting the 3-chain subset would be choosing the
# chains that give a tidy answer.
#
# Pinning b_shape is already known to harden the geometry: fix_bshape found
# the same thing in simulation, where two seeds both stalled just above the
# convergence gate. Pinning at 400, very tight synchrony, AND raising n_c
# appears to compound it. Plausibly phase multimodality -- sharply
# synchronised parasites with a slightly mismatched period admit more than one
# local phase alignment, and b_offset is a unit_vector, so those are separate
# modes rather than a ridge.
#
# Changes two things at once, deliberately, because the aim is to get a
# converged answer rather than to attribute the failure: adapt_delta 0.95 with
# max_treedepth 12 (entries 41-42), and a different seed.
#
# IF A CHAIN STILL STICKS at a similar lp gap, the mode is real. Report it as
# multimodality, show both modes, and do not chase it with a third
# configuration.
#
# Walltime is 3 days: entry 40 took 24.5 h at adapt_delta 0.8, and raising it
# to 0.95 makes each gradient path longer.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-nc-bs400-retry.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=41-42
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=wock-ncbs2
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-ncbs2-%a.out
#SBATCH --error=wock-ncbs2-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_SEED=20261007

Rscript --vanilla ../_scripts/wockner-fit.R
