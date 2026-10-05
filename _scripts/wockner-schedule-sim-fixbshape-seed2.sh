#!/bin/bash -l

# fix_bshape, reps 2 and 3, reseeded.
#
# The first pass (SLURM 29628) returned a clean and consistent paired result --
# fix_bshape minus default of -0.92, -1.15, -0.97 h across reps 1-3 -- but
# fixing b_shape degraded sampling in every replicate relative to `default`
# (divergences 53->67, 25->586, 35->222). Reps 2 and 3 are not usable as
# posteriors on their own terms: rep2 is max R-hat 1.343 with min ESS 10.2 and
# 14.7% divergences, rep3 is R-hat 1.066 with 5.5% divergences.
#
# This reseeds the sampler on those two with SCHEDSIM_FIT_SEED=2, leaving the
# simulated data identical (seed_noise is unchanged), so they still pair against
# the same `default` fits. The reseed is triggered by the diagnostics alone and
# would have been run whichever way the biases pointed; precedent is the
# daH2_no_pool reseed, where R-hat 1.154 -> 1.017 moved a rung from z -3.08 to
# z -1.78 and flattened the horizon ladder.
#
# Configs land as fix_bshape-rep2-seed2 and fix_bshape-rep3-seed2, so the
# originals are not overwritten and both can be compared.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-fixbshape-seed2.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=86-87
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-fixbs-s2
#SBATCH --output=_data/wock-fixbshape-seed2-%a.out
#SBATCH --error=_data/wock-fixbshape-seed2-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export SCHEDSIM_FIT_SEED=2

Rscript --vanilla _scripts/wockner-schedule-sim.R
