#!/bin/bash -l

# The NULL RUN the b_shape regression needs.
#
# A regression test compares a posterior shift against something. Comparing it
# against Monte Carlo standard error is anti-conservative: two runs also
# differ in step-size and mass-matrix adaptation, which MCSE does not capture.
# The honest reference is two fits of IDENTICAL code and data differing only
# in the sampler seed, which is this.
#
#   np_wide_total0          old code, seed 538065874   (the reference)
#   np_wide_total0-rebuild  new code, seed 538065874   (the test, SLURM 29526)
#   np_wide_total0-null     new code, seed 1618033989  (this job, the NULL)
#
# reference-vs-rebuilt is the test: it differs by the code change AND by the
# sampler trajectory, which a recompile alone changes. rebuilt-vs-null differs
# by the trajectory alone. So the null's max |z| is what the test's max |z|
# has to be read against.
#
# WHAT IT DECIDES. The test put max |z| at 3.52, on b_shape[2], with 17.5% of
# entries beyond 2. b_shape under a free wide prior is the most heavy-tailed
# quantity in this model, so its posterior MEAN is a high-variance statistic
# and se_mean understates run-to-run spread worst exactly there -- but it is
# also the parameter the change touches. If the null reaches a comparable
# figure, the test is Monte Carlo variation and the change is inert. If the
# null sits near 2, the test found a real shift in b_shape.
#
# Task 9 is np_wide_total0, the config the test used. WOCKFIT_SUFFIX keeps
# this beside the others instead of overwriting either.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-null.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=9
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-fit-null
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-fit-null-%a.out
#SBATCH --error=wock-fit-null-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_SEED=1618033989
export WOCKFIT_SUFFIX=-null

Rscript --vanilla ../_scripts/wockner-fit.R
