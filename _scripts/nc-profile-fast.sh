#!/bin/bash -l

# Thread 15: the dispersion profile, densely, via the convolution forward map.
#
# WHY THIS AND NOT 30527. The convolution in
# _scripts/convolution-forward-map.R is the SAME model as mat_exp_series,
# agreeing to 1e-12, and a whole unit fit takes 12 s at n_c = 384 where SLURM
# 29684 took 1777 s. 30527 is still running and is deliberately NOT cancelled:
# it computes three of these rungs on the production mat_exp_series path, and
# this script checks itself against 29684's saved log-likelihoods at 96, 192
# and 384 and STOPS if any differs by more than 0.01, so the two paths keep
# each other honest.
#
# WHAT IT ADDS.
#   A. The exact chain at nine integer rungs, 64 to 1024, instead of three.
#   B. The gamma-kernel IPM at a FIXED mesh of 192 with CONTINUOUS n_eff over
#      sixteen values from 48 to 4096. Gamma(shape = n_eff*t/cl, scale =
#      cl/n_eff) has the chain's mean, variance and right skew, so this is the
#      forward map an IPM would actually use. It is the only way to ask
#      whether the dispersion is IDENTIFIED rather than merely quantised.
#   C. A mesh-convergence check: three n_eff values repeated at M = 384. If
#      those shift much, the IPM rows are mesh artefacts and not results.
#
# THE READING RULE IS NOT NEW. It is the one pre-registered in
# _scripts/nc-dispersion-profile-read.R before any of this existed: turns over
# / saturates / still climbing / flat, judged on whether the 2-log-likelihood
# interval is bounded in the fine direction.
#
# COST. Measured: 6, 10 and 12 s for one unit at n_c = 96, 192 and 384, on the
# unit with 6 distinct observation times. Scaling by distinct times, the worst
# unit (Mefloquine, 10 times) is about 9 min for all 28 rungs. Walltime is 8 h
# because cost estimates in this project have run low twice.
#
# RESOURCES. 14 x 1 CPU x 2G = 14 CPUs and 28,672 MB. Alongside 30524, 30525
# and 30527 that totals 46 of 128 CPUs and 397,312 of 515,670 MB, so 36% and
# 77% of the half-node budget in CLAUDE.md. No --array=1-N%M throttle needed.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/nc-profile-fast.sh
#   Rscript --vanilla _scripts/nc-profile-fast-read.R | tee _data/nc-profile-fast-$(date +%F).txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-14
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
#SBATCH --job-name=ncprof-fast
#SBATCH --output=_data/nc-profile-fast-%a.out
#SBATCH --error=_data/nc-profile-fast-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/nc-profile-fast.R
