#!/bin/bash -l

# The horizon ladder: separate "how far ahead" from "how much was withheld".
#
# Design A gave a mask-dependent answer -- pooled_cl 3.27 log units behind
# no_pool at the last-third mask, 0.04 behind at the last-quarter. Decomposing
# the third-mask fits showed 89% of the deficit sits on the points the quarter
# mask ALSO holds out, so it is not a horizon effect: the two differ in what
# the models were TRAINED on, not what they were scored on.
#
# This ladder isolates that. k = 1, 2, 3 points from the end of each series,
# capped at n - 3 so no series keeps fewer than 3, holding out 177, 350 and
# 479 of 1130. The masks NEST, so the k=1 points are held out by all three
# rungs -- score every rung on THAT common window and the scored observations
# are fixed while only the training set varies.
#
#   31,32 daH1 (k=1)   33,34 daH2 (k=2)   35,36 daH3 (k=3)
#
# Two models only, no_pool and pooled_cl. pooled_R and pooled_both would
# double the cost to answer a different question.
#
# Everything else matches Design A: corrected log10_total0 prior, b_shape
# FIXED at 400, calc_log_lik on.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-horizon.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=31-36
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-horizon
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-horizon-%a.out
#SBATCH --error=wock-horizon-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla ../_scripts/wockner-fit.R
