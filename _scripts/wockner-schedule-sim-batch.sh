#!/bin/bash -l

# Replicates 6-7 of the six budget arms, and replicates 4-5 of the two
# bound-geometry arms.
#
# The paired cycle-length budget is the widest error bar left in this project.
# At n=4 its own spread across replicates is ~0.3 h, against component effects
# of 0.43 h (b_shape prior) and 0.10 h (log10_total0 prior) -- so the
# decomposition is still comparing effects to a noise floor of the same size.
# This is the one open number more compute can narrow; the ~1.5 h of
# unexplained cycle-length bias is not fit-limited and these do not touch it.
#
# The bound arms are at n=2 and n=1 after three of six replicates failed the
# convergence gate, so they get replicates 4-5 here and new-seed refits of the
# failures in wockner-schedule-sim-refit2.sh.
#
# Tasks are rows of expand_grid(arm, rep) with 9 arms x 7 replicates; adding
# the two seeds SHIFTED every index above arm 1. Regenerate the table rather
# than reading indices off by hand.
#   6,7 default      13,14 wide         27,28 no_hier
#   34,35 wide_bshape 41,42 wide_total0 48,49 wide_nuis
#   53,54 cl_move     60,61 cl_wide_move
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-batch.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=6,7,13,14,27,28,34,35,41,42,48,49,53,54,60,61
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-ss-batch
#SBATCH --output=_data/wock-ss-batch-%a.out
#SBATCH --error=_data/wock-ss-batch-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/wockner-schedule-sim.R
