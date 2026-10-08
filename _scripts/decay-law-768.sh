#!/bin/bash -l

# Thread 15, the deciding rung: extend model A's n_c profile to 768.
#
# WHAT THIS DECIDES. The IPM-versus-new-parameter choice turns on whether the
# data's preferred dispersion sits AT the Erlang floor or BELOW it. In the
# chain, transit time has mean cycle_length and coefficient of variation
# 1/sqrt(n_c), and for a fixed stage count Erlang is the MINIMUM-variance
# case, so 1/sqrt(n_c) is a floor the family cannot go under. If the profile
# has plateaued, the floor is reachable and no structural rewrite is
# justified. If it is still climbing, the family is fighting the data -- they
# want near-deterministic transit, which Erlang only approaches as n_c -> inf
# at cubic cost -- and an IPM's free dispersion width is the point.
#
# PRE-REGISTERED READING RULE, fixed before submission. The gains so far are
# +72.4 summed log-likelihood units for 96 -> 192 and +10.4 for 192 -> 384, a
# ratio of 0.14, so geometric decay predicts +1.5 for 384 -> 768.
#   - summed gain < +3 AND fewer than 9/14 units positive -> plateau is real,
#     no rewrite, pick n_c by elpd and optionally free the stage rates.
#   - summed gain >= +8 AND >= 10/14 units positive -> gains are not decaying,
#     the Erlang family is fighting the data, an IPM is justified on evidence.
#   - anything between -> ambiguous, and the decision goes back to biology.
#
# Model A only. Model B is irrelevant to this question and costs 7x.
#
# COST. One trajectory takes 1.925 s at n_c = 384 and 14.894 s at 768
# (measured 2026-10-08, ratio 7.74, consistent with cubic). Model A at 384
# took 1777 s per unit in SLURM 29684, so expect ~3.8 h per unit and 2-7 h
# across units, which differ about twofold. The 2-day walltime is ample.
# Do NOT extrapolate this to a production Stan fit: at 768 that would be
# ~540 h.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/decay-law-768.sh
#   Rscript --vanilla _scripts/nc-768-read.R

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-14
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=decay768
#SBATCH --output=_data/decay-law-768-%a.out
#SBATCH --error=_data/decay-law-768-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export DECAY_MODE=a768
Rscript --vanilla _scripts/decay-law-test.R
