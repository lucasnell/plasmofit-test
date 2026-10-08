#!/bin/bash -l

# Thread 15: re-run the dense dispersion profile with EIGHT optimiser starts.
#
# WHY. SLURM 30576 ran with two starts and its profiles came back rougher than
# the 2-log-likelihood currency the pre-registered reading rule is written in
# (2.39 for the chain, 7.17 for the gamma IPM), so the reader returned NO
# VERDICT. The cause was MEASURED, not assumed:
# _scripts/profile-noise-check.R refit n_eff = 4096 with eight starts and
# DSM265|1800 alone gained +11.29 log-likelihood units -- the whole of the
# 11.10 pooled drop that had looked like a turnover.
#
# WHAT IS AND IS NOT CHANGED. Eight starts instead of two, and NE_CHECK fixed
# to use values that actually appear in NE_IPM (they did not, so two of three
# mesh comparisons paired against nothing and the reader crashed on NA). The
# b_shape cap stays at 5000. Raising it is a separate choice with content --
# production pins b_shape at 400 with max_shape 1000, so a screen reaching
# 5000 is already outside the production range -- and is deliberately NOT
# bundled in here, so this run changes one thing.
#
# NOTHING IS OVERWRITTEN. NCPF_TAG=-s8 writes beside the two-start results,
# which stay as the baseline: the difference between the two IS the
# measurement of how much the optimiser was stopping short. The script also
# refuses to write over an existing output unless NCPF_FORCE is set.
#
# COST. 30576 ran 17-38 min per task with two starts, and each start is an
# independent optim() call, so eight is about 4x: 70-150 min. Walltime is 12 h,
# which is the measured worst case times three -- every cost extrapolation in
# this project has come in low (see claude/gotchas.md).
#
# RESOURCES. 14 x 1 CPU x 2G = 14 CPUs and 28,672 MB. Alongside 30524, 30525
# and 30527 that is 46 of 128 CPUs and 397,312 of 515,670 MB, inside the
# half-node rule in CLAUDE.md.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/nc-profile-fast-s8.sh
#   Rscript --vanilla _scripts/nc-profile-fast-read.R | tee _data/nc-profile-fast-$(date +%F).txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-14
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --job-name=ncprof-s8
#SBATCH --output=_data/nc-profile-fast-s8-%a.out
#SBATCH --error=_data/nc-profile-fast-s8-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export NCPF_STARTS=8
export NCPF_TAG=-s8

Rscript --vanilla _scripts/nc-profile-fast.R
