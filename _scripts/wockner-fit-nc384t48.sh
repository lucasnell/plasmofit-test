#!/bin/bash -l

# Entry 54: the free-b_shape rung at n_c = 384 (entry 48) at entry 42's tuned
# sampler settings, adapt_delta 0.95 and max_treedepth 12. The tuned twin the
# insurance job (SLURM 30837, entries 49-53) left out. Added 2026-10-10 at
# Lucas's request; cancel it with 30837 if 43-48 converge.
#
# RUNS A SECOND COPY, _scripts/wockner-fit-nc384t48.R: wockner-fit.R and
# wockner-fit-nc384t.R were both being read by running jobs. Identical to
# wockner-fit-nc384t.R through entry 53.
#
# COST. 7-day walltime, as for 42 and 49-53.
#
# FOOTPRINT. 1 task x 4 CPUs x 24G. With 30830 and 30837 the total is 52 CPUs,
# 319,488 MB, inside the 128 CPU / 515,670 MB half-node budget.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-nc384t48.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=54
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=wock-nc384t48
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-nc384t48-%a.out
#SBATCH --error=wock-nc384t48-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export WOCKFIT_SEED=20261010

Rscript --vanilla ../_scripts/wockner-fit-nc384t48.R
