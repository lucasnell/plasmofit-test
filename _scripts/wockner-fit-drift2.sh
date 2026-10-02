#!/bin/bash -l

# Thread 6's second drift measurement, so it is n=2 rather than n=1.
#
# np_bs400 (task 19) rebuilt from the PRE-CHANGE source and refitted with the
# same seed, exactly as wockner-fit-drift.sh does for np_wide_total0. Task 19
# passes only arguments the pre-change archer_stan_data() already had --
# mean_log_b_shape, sd_log_b_shape, max_shape, sd_log10_total0 -- so it runs
# unmodified against that build.
#
# It is the useful second case because b_shape is PINNED there. If drift is a
# property of the weakly identified directions, it should appear whether
# b_shape is free or not, as the original excess did.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-drift2.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=19
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-drift2
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-drift2-%a.out
#SBATCH --error=wock-drift2-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_LIB=/home2/lan68/plasmofit/.prechange/lib
export WOCKFIT_SUFFIX=-drift

Rscript --vanilla ../_scripts/wockner-fit.R
