#!/bin/bash -l

# End-to-end regression for the package's b_shape change (branch
# fixed-b-shape, package commits 0bf12b0 and e637d8d).
#
# The package test suite passes and a fit with b_shape = 400 returns
# b_shape = 400.0 exactly, so the NEW path works. What nothing has checked is
# that the change left the OLD path alone. Both configs here use b_shape = NULL
# and so should be unaffected by it:
#
#   9  np_wide_total0  b_shape free, corrected log10_total0 prior
#   19 np_bs400        b_shape pinned by a tight prior, not by the new argument
#
# WOCKFIT_SUFFIX makes these land as `<config>-rebuild`, BESIDE the existing
# fits rather than on top of them. Re-running a config under its own name
# would overwrite the baseline being compared against, which is the whole
# point of the exercise.
#
# HOW TO READ IT. The Stan models were recompiled, and per claude/gotchas.md
# two fits of the same model are never bit-identical across a recompile:
# rebuilding reorders floating-point operations and HMC turns a last-bit
# difference into a different trajectory within a few leapfrog steps. So the
# rebuilt fit is effectively an INDEPENDENT run of the same specification,
# and the comparison is a posterior one, scaled by Monte Carlo error --
# the method in _scripts/wockner-anchor-regression.R. Do not expect, or
# demand, identical draws.
#
# A shift of a few MCSE in any parameter is the expected result. A shift far
# outside that, or a change in sampler health, is the thing this is looking
# for.
#
# NOTE neither config exercises the new argument. Testing that path means a
# new CONFIGS entry passing b_shape = 400 through archer_stan_data(), which
# should then reproduce np_bs400 -- a tight prior at sd 0.05 is nearly the
# same model. That is a separate run and is not in this job.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-regression.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=9,19
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-fit-regr
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-fit-regr-%a.out
#SBATCH --error=wock-fit-regr-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_SUFFIX=-rebuild

Rscript --vanilla ../_scripts/wockner-fit.R
