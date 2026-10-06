#!/bin/bash -l

# Thread 5 in simulation: what does misspecifying n_c do to cycle_length?
#
# The real-data ladder (SLURM 29633) found cycle_length moving -3.20 h from
# n_c = 96 to 384, with b_shape falling alongside, and n_c = 96 losing 71.8
# elpd. Real data cannot say which end is right: there is no truth to miss,
# and elpd scores prediction on the observed window, which an amplitude/period
# tradeoff can win while getting the period wrong.
#
# These arms put a truth in. Generating n_c and fitting n_c are decoupled, so
# with `default` (96/96) already run the four cells are:
#
#              fit 96              fit 192
#    gen  96   default             sim96_fit192
#    gen 192   sim192_fit96        sim192_fit192
#
# The diagonal is correctly specified and measures estimator bias at each n_c
# -- `default` is +1.97 h. The off-diagonal measures what misspecifying n_c
# manufactures.
#
# sim96_fit192 is the cell that bears on the real-data result most directly:
# if fitting ABOVE the true n_c drags cycle_length down, the -3.20 h on real
# data may be an artefact of over-large n_c rather than a correction.
#
# Tasks: 92-94 (sim192_fit96), 99-101 (sim192_fit192), 106-108 (sim96_fit192),
# reps 1-3, the replicates with a converged `default` to read against.
# Fits at n_c = 192 run roughly 3x the n_c = 96 cost, so ~3 h rather than ~1 h.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-nc2x2.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=92-94,99-101,106-108
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-nc2x2
#SBATCH --output=_data/wock-nc2x2-%a.out
#SBATCH --error=_data/wock-nc2x2-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/wockner-schedule-sim.R
