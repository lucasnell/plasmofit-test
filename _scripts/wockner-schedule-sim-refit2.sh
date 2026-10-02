#!/bin/bash -l

# New-seed refits of the four replicates still failing the convergence gate.
#
#   task  4  default-rep4       R-hat 1.06 on seed 1618033989, 1.22 before it
#   task 51  cl_move-rep2       R-hat 1.32, 513 divergences, min ESS 10
#   task 58  cl_wide_move-rep2  R-hat 1.07
#   task 59  cl_wide_move-rep3  R-hat 1.06
#
# default-rep4 is the one that matters most: `default` is the baseline every
# paired contrast subtracts, so without it replicate 4 contributes nothing to
# the budget however many other arms converged. It has now failed on two
# seeds, though the second was a near miss -- divergences fell 478 -> 39 and
# ESS rose 13 -> 117 -- so a third seed is reasonable rather than stubborn.
#
# IF IT FAILS A THIRD TIME, stop reseeding and record it: that would be real
# evidence the dataset is hard, which would contradict thread 9's finding on
# the two earlier cases and is worth knowing in its own right.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-refit2.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=4,51,58,59
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-ss-refit2
#SBATCH --output=_data/wock-ss-refit2-%a.out
#SBATCH --error=_data/wock-ss-refit2-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export SCHEDSIM_FIT_SEED=1123581321

Rscript --vanilla _scripts/wockner-schedule-sim.R
