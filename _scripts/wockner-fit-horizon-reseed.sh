#!/bin/bash -l

# Reseed daH2_no_pool, the one horizon rung that failed the convergence gate.
#
# The ladder scored k = 1, 2, 3 on a common 177-point window and read null at
# k = 1 (z -0.46) and k = 3 (z +0.70). Only k = 2 showed a signal (z -3.08)
# and only k = 2's no_pool fit failed: max R-hat 1.15, 312 divergences. By the
# project's gate that row is void, so the rung that disagrees with the other
# two is the rung that cannot be read.
#
# Output lands as daH2_no_pool-seed2, beside the failed fit rather than over
# it. Then put it into _scripts/wockner-horizon-score.R's rung 2 in place of
# the original and re-run.
#
# If it converges and still reads -1.7, the ladder is non-monotone and that is
# a finding. If it reads null like the other two, the ladder is flat and
# Design A's third-mask result has no support anywhere.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-horizon-reseed.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=33
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-h2reseed
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-h2reseed-%a.out
#SBATCH --error=wock-h2reseed-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_SEED=1123581321
export WOCKFIT_SUFFIX=-seed2

Rscript --vanilla ../_scripts/wockner-fit.R
