#!/bin/bash -l

# Thread 2, bound geometry: the last cheap untested suspect for the ~1.5 h of
# cycle-length bias that survives correcting the nuisance priors.
#
# The default window is [35, 50] against a truth of 45.012 h -- 5 h above and
# 10 h below -- and the prior on cycle_length is normal on the LOGIT scale
# between those bounds, so an asymmetric window is also an asymmetric prior in
# hours. Two arms, both MOVING the bounds rather than widening them alone,
# since max_cl = 55 is what reintroduced the boundary mode:
#
#   cl_move       [37.5, 52.5]  same 15 h width, centred on the truth
#                               -> isolates ASYMMETRY at constant width
#   cl_wide_move  [30, 60]      moved and widened, as thread 2 proposed
#                               -> asks whether distance from any bound matters
#
# Read together. Only cl_wide_move moving the bias means it is distance from
# the bounds; both moving it means it is the asymmetry; neither means bound
# geometry is not the remaining explanation and the ~1.5 h is elsewhere.
#
# Replicates 1-3, which are the ones with a converged `default` on the same
# simulated dataset, so each arm can be paired against it. Tasks are rows of
# expand_grid(arm, rep) with 9 arms x 5 replicates:
#   36-38 = cl_move reps 1-3, 41-43 = cl_wide_move reps 1-3.
#
# cl_wide_move has a lower min_cl, which widens the Erlang window and costs
# compute; allow it to run longer than the other arms.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-bounds.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=36-38,41-43
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-schedsim-bounds
#SBATCH --output=_data/wock-schedsim-bounds-%a.out
#SBATCH --error=_data/wock-schedsim-bounds-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/wockner-schedule-sim.R
