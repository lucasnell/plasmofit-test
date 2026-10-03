#!/bin/bash -l

# Design A at a QUARTER mask -- the pre-registered sensitivity to 23-26.
#
# Identical to wockner-fit-designA.sh except the mask: 218 of 1130
# observations held out instead of 306, every series still keeping 3-6
# points. Both fractions are run so the fraction cannot be chosen after
# seeing the answer. Read the two side by side; a verdict that holds at one
# fraction and not the other is a finding about the design's sensitivity, not
# a reason to prefer whichever reads better.
#
# Every no_pool vs pooled_cl verdict in this project rests on a comparison its
# own documentation calls weak. Observation-level loo holds out one point at a
# time, and a held-out point is pinned down by its neighbours whatever the
# model. Trial-level loo is the comparison meant to settle it and it is broken
# here: Pareto k > 0.7 for all 13 trial units, in every model, so the
# importance-sampling approximation has failed rather than the models being
# indistinguishable.
#
# So: mask the last QUARTER of every series, fit, and score the four model
# variants on the window none of them saw. 306 of 1130 observations held out,
# 912 fitted, every series keeping 3-6 points. Each series' initial
# conditions and error scale stay informed, and no_pool can adapt eta_cl per
# trial where pooled_cl cannot -- so a cycle-length error accumulates as phase
# drift exactly inside the held-out window, which is what the comparison is
# meant to be sensitive to.
#
#   27 daQ_no_pool      28 daQ_pooled_cl
#   29 daQ_pooled_R     30 daQ_pooled_both
#
# All four carry the defensible nuisance priors (sd_log10_total0 = 1) and a
# b_shape FIXED at 400, the ladder's plateau. b_shape is not identified by
# these data, and leaving it free would give each model a different amount of
# freedom in a comparison that is not about b_shape.
#
# READ IT with the held-out elpd, not the total: sum log_lik over the
# observations with hold_out == 1 only. generated quantities computes log_lik
# for every observation whether or not it was fitted, which is the whole point
# -- the total would mix fitted and held-out points and answer nothing.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-designQ.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=27-30
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-designQ
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-designQ-%a.out
#SBATCH --error=wock-designQ-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla ../_scripts/wockner-fit.R
