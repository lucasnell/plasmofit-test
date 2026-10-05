#!/bin/bash -l

# fix_bshape: is the cycle-length bias the cost of not knowing synchrony?
#
# The profile scans put only about a fifth of the bias in the likelihood with
# every nuisance held at the truth (+0.41 h of +1.97 h, varying by replicate).
# The other four fifths appear only when the nuisances are ESTIMATED rather
# than known, and the prior-location arms account for 0.44 h of that. So the
# question is which nuisance's estimation carries the rest.
#
# This fixes b_shape at the value the data were simulated from -- 14 values,
# range 10.41-21.78 -- so the arm differs from `default` in that alone and
# pairs against it on the same simulated datasets. b_shape is the parameter
# this project has repeatedly found unidentified: the design cannot separate
# 15 from 65, and its posterior mean runs to 65 with a wide prior.
#
# `fix_sd` is the same test for the error scale and carries none of the bias
# (+0.104 h, wrong sign). If this arm carries a large share, the bias is the
# cost of not knowing synchrony. If it does not, then no single nuisance
# carries it and the cost is joint -- which would be a different and harder
# result, and would make SBC the next instrument rather than more arms.
#
# Tasks 85-87 are replicates 1-3, the ones with a converged `default` to pair
# against. Diagnostic only: fixing a parameter at the truth is available in
# simulation and nowhere else.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-fixbshape.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=85-87
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-fixbshape
#SBATCH --output=_data/wock-fixbshape-%a.out
#SBATCH --error=_data/wock-fixbshape-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/wockner-schedule-sim.R
