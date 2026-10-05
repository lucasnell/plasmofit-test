#!/bin/bash -l

# The cycle-length prior's LOCATION, which has never been varied.
#
# On the logit scale between [35, 50], a normal centred at cl_prior_center =
# 48 implies a prior MEAN of 47.44 h, which is 2.43 h above the simulated
# truth of 45.012. Widening sd_logit_cl leaves the median pinned at 48.00 and
# drags the mean down only slowly through the bounds -- so the result that
# closed this question, "quadrupling the prior variance moves the estimate
# 0.15 h", varied the prior's WIDTH and never its LOCATION. The weight on
# record (+0.271 h) is a normal-normal approximation on the logit scale,
# where the asymmetry producing that 2.43 h gap does not appear.
#
# Three points, reading these two arms together with `default`:
#
#   default           centre 48       prior mean 47.44 h    +2.43 h
#   cl_center_truth   centre = truth  prior mean 44.61 h    -0.41 h
#   cl_center_low     centre 42       prior mean 42.09 h    -2.93 h
#
# Both displacements pre-flighted against the built data list, not assumed.
#
# READ IT as a SLOPE: posterior cycle-length bias against prior-mean
# displacement, across those three arms, paired within replicate. That slope
# IS the prior's weight, measured rather than approximated. If it is near the
# 0.15 the width test implied, the prior explains ~0.36 h of the ~1.5 h and
# thread 2 keeps most of its gap. If it is much larger, the prior explains
# more of the bias than this project has believed since the schedule-bias
# simulation, and the ~0.15 h figure was an artefact of testing width.
#
# DIAGNOSTIC, like fix_sd. Centring a prior on the truth is available in
# simulation and nowhere else; this measures a weight, it does not propose a
# prior. Deciding cl_prior_center for real data is thread 7 and is a
# different question.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-center.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=71-73,78-80
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-center
#SBATCH --output=_data/wock-center-%a.out
#SBATCH --error=_data/wock-center-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/wockner-schedule-sim.R
