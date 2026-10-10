#!/bin/bash -l

# Insurance for SLURM 30830's b_shape ladder at n_c = 384: the same five
# pinned rungs (15, 50, 100, 250, 600) at entry 42's tuned sampler settings,
# adapt_delta 0.95 and max_treedepth 12. Entries 49-53. Submitted 2026-10-10
# alongside 30830 so a failure of 43-47 at default settings does not cost
# another round of waiting. Lucas will cancel this job if 43-47 converge.
#
# RUNS A COPY, _scripts/wockner-fit-nc384t.R, not wockner-fit.R: 30830 was
# reading the original when these were added, and editing a script under a
# running job corrupts it (claude/gotchas.md). The copy is identical through
# entry 48.
#
# Sampler settings change exploration, not the posterior, so a converged rung
# here is comparable with converged rungs from 30830. Read it the same way:
# lp__ by chain, then max R-hat < 1.05, then paired PSIS-LOO across converged
# rungs. The free rung (48) has no tuned twin.
#
# COST. Like entry 42, 2.5-5x the ~20 h of a default-settings fit is
# plausible, so walltime is 7 days.
#
# FOOTPRINT. 5 tasks x 4 CPUs x 24G = 20 CPUs, 122,880 MB. With 30830's 7
# tasks the total is 48 CPUs, 294,912 MB, inside the 128 CPU / 515,670 MB
# half-node budget with no throttle.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-nc384t-ladder.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=49-53
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=wock-nc384t
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-nc384t-%a.out
#SBATCH --error=wock-nc384t-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_SEED=20261010

Rscript --vanilla ../_scripts/wockner-fit-nc384t.R
