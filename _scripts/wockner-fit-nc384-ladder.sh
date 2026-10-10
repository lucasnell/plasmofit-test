#!/bin/bash -l

# n_c = 384 after entry 40 failed twice: the tuned retry (entry 42) and a
# b_shape ladder at n_c = 384 (entries 43-48), submitted together 2026-10-10.
#
# WHY. Entry 40 (np_bs400_nc384) failed in 29635 and again on its reseed,
# SLURM 30829: lp__ by chain -895.5, -973.1, -1065.3, -954.2, max R-hat 7.74,
# treedepth saturated in 34.0% of transitions. All four chains apart, a
# general mixing failure, so a third reseed is not the fix.
#
#   42     np_bs400_nc384_t: adapt_delta 0.95, max_treedepth 12. Justified by
#          the 34% saturation at THIS config (the earlier "less justified" note
#          inferred from n_c = 192's 2%). Answers: can the 400 rung converge
#          at all?
#   43-48  b_shape 15, 50, 100, 250, 600 pinned, and free, at n_c = 384 with
#          DEFAULT sampler settings. Answers: is the 400 pin, chosen at
#          n_c = 96, what makes 384 hard, and where does elpd plateau at 384?
#          Rung choice and prediction are recorded in wockner-fit.R above the
#          configs, before any of these ran.
#
# HOW TO READ. Per fit: lp__ by chain, then max R-hat < 1.05, then the rest.
# A failed rung is not a posterior and is not used. Compare converged rungs by
# paired observation-level PSIS-LOO, and report cycle_length as its spread
# across converged rungs, not the value at the best one.
#
# SEED. One new seed for all seven. Output names are distinct from every
# existing fit (checked 2026-10-10), so no WOCKFIT_SUFFIX is needed.
#
# COST. 30829 (entry 40, default settings) took 20:21:41. Ladder rungs should
# be of that order. Entry 42 at max_treedepth 12 could cost 2.5-5x, 51-102 h,
# which the old 3-day walltime brackets, so walltime is 7 days; the partition
# has no limit. Every cost extrapolation here has come in low.
#
# FOOTPRINT. 7 tasks x 4 CPUs x 24G = 28 CPUs, 172,032 MB, inside the
# 128 CPU / 515,670 MB half-node budget with no throttle. Check squeue -u lan68
# before submitting anything else beside it.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-nc384-ladder.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=42-48
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=wock-nc384
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-nc384-%a.out
#SBATCH --error=wock-nc384-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_SEED=20261010

Rscript --vanilla ../_scripts/wockner-fit.R
