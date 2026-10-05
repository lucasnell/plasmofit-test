#!/bin/bash -l

# Thread 5: the n_c / b_shape confound, on real data.
#
# n_c is the number of sequential exponential compartments, so transit over one
# cycle is Erlang(n_c, n_c / cycle_length) and the stage distribution's sd after
# k cycles is sqrt(k / n_c) CYCLES. Over the Wockner window (k = 4.80) that is
# 0.224 cycles at n_c = 96, 0.158 at 192, 0.112 at 384. n_c alone therefore
# fixes how fast the model lets the parasites desynchronise, and
# check_erlang_window() chooses it for numerical accuracy, never for biology.
#
# b_shape is the only free synchrony knob and the chain's own decay dominates
# the late-window amplitude, so a wrong fixed decay rate has to land in
# b_shape. fix_bshape (SLURM 29628) showed errors landing in b_shape propagate
# into cycle_length at -0.92 h, which is what makes this worth running.
#
# Entries 37 (n_c = 192) and 38 (n_c = 384) against the n_c = 96 baseline
# np_wide_both (entry 11), which has identical priors. The three are a ladder
# in n_c and nothing else.
#
# PREDICTION, recorded before the fits are run and repeated in wockner-fit.R:
# if b_shape has been absorbing a too-fast desynchronisation rate, b_shape
# should FALL monotonically along the ladder toward the 14-19 the tight prior
# gave. If it is pinned by something real, it stays at 65+. Non-monotonic means
# neither. cycle_length is watched alongside: the confound only matters if
# moving n_c moves the headline number.
#
# Cost measured, not guessed: _scripts/nc-sizing.R puts the per-leapfrog cost
# at 1.81x going 96 -> 192, against np_wide_both's 1.58 h slowest chain, so
# ~2.9 h at 192. 384 is extrapolated and could be ~6 h; the walltime below is
# generous because the extrapolation is the uncertain part.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-nc.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=37-38
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-nc
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-nc-%a.out
#SBATCH --error=wock-nc-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla ../_scripts/wockner-fit.R
