#!/bin/bash -l

# Task 22, np_bs400_data: the first fit that actually exercises
# archer_stan_data(b_shape = ), the new argument on the package's
# fixed-b-shape branch.
#
# Regression target is np_bs400 (task 19), which pins b_shape at 400 with a
# tight prior instead. sd_log_b_shape = 0.05 is +-10% at 95%, worth under
# 0.2 h of starting-stage spread at this centre, so the two are nearly the
# same model and should agree to within Monte Carlo error. They will NOT
# agree exactly, and not only because of the recompile: a tight prior is not
# a point mass.
#
# WHY THIS RUNS FROM A COPY OF wockner-fit.R.
#
# SLURM 29526 was reading _scripts/wockner-fit.R when this config was added,
# and Rscript reads a source file incrementally -- editing it under a running
# job shifts byte offsets and the job parses garbage, which has already eaten
# two 1.5 h fits here. claude/gotchas.md's sanctioned workaround is to copy
# the script and submit the copy, which is what this is.
#
# **_scripts/wockner-fit-bshapedata.R IS TRANSIENT.** Once the queue drains,
# fold `np_bs400_data` into wockner-fit.R's CONFIGS as entry 22 -- the copy
# is identical to it apart from that entry -- and delete both this script and
# the copy. Until that is done, wockner-fit.R cannot reproduce
# wock-fit-np_bs400_data.rds, which is exactly the provenance gap the project
# conventions exist to prevent.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-bshapedata.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=22
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-fit-bsdata
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-fit-bsdata-%a.out
#SBATCH --error=wock-fit-bsdata-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla ../_scripts/wockner-fit-bshapedata.R
