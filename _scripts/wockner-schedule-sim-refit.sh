#!/bin/bash -l

# Refit the three replicates from SLURM 28940 that failed the convergence
# gate, on a different sampler seed.
#
#   default-rep4      max R-hat 1.22, 478 divergences, min ESS 13
#   default-rep5      max R-hat 1.07, min ESS 57
#   wide_total0-rep4  max R-hat 1.06, min ESS 63
#
# Thread 9's finding is that both previous failures of this kind were bad FIT
# SEEDS rather than hard datasets, so this changes the seed and nothing else:
# the data and the truth are rebuilt deterministically from the same noise
# seeds. Output names carry the seed, so the failed runs stay on disk beside
# the new ones rather than being overwritten.
#
# This matters more than two lost arms. `default` is the baseline every
# paired contrast in the cycle-length budget subtracts, so with reps 4 and 5
# missing the batch added no PAIRED replicates at all and the budget is still
# n=3 (see findings.md, "Replicates 4-5").
#
# Tasks are rows of expand_grid(arm, rep) with 7 arms x 5 replicates:
#   4 = default-rep4, 5 = default-rep5, 29 = wide_total0-rep4.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-refit.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=4,5,29
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-schedsim-refit
#SBATCH --output=_data/wock-schedsim-refit-%a.out
#SBATCH --error=_data/wock-schedsim-refit-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export SCHEDSIM_FIT_SEED=1618033989

Rscript --vanilla _scripts/wockner-schedule-sim.R
