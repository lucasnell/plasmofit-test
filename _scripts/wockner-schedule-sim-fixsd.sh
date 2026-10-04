#!/bin/bash -l

# The `fix_sd` arm: what does estimating the error scale cost?
#
# findings.md's ladder on the same simulated data reads
#
#   per-trial MLE            -0.195 h
#   pooled MLE               +0.544 h
#   fitted posterior mean    +1.56 to +1.71 h
#
# so the likelihood is not biased and the posterior is. Widening the nuisance
# priors accounts for ~0.5 h of the gap. Two named differences between the
# MLE and the fit remain, and the MAP check could not settle the first
# (marginalisation): the joint surface has no usable mode. This arm settles
# the second -- the MLE fixes sd_iRBC at the truth where the fit estimates it.
#
# `fix_sd` passes the simulation's own sd_iRBC (27 values, 0.396-0.786) into
# archer_stan_data(), so it differs from `default` in that alone and pairs
# against it on the same simulated datasets.
#
# Tasks 64-66 are replicates 1-3, the ones with a converged `default` to pair
# against. Diagnostic only: fixing an error scale at a value taken from the
# truth is available in simulation and nowhere else.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-fixsd.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=64-66
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-fixsd
#SBATCH --output=_data/wock-fixsd-%a.out
#SBATCH --error=_data/wock-fixsd-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/wockner-schedule-sim.R
