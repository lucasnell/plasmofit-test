#!/bin/bash -l

# 2-D profile of the sequestration schedule, midpoint x slope; design and
# reading rule in _scripts/start-seq2.R. 91 tasks: 13 per trial for entries
# 54 and 57 over the whole 5 x 5 grid, then one per (trial, midpoint) for
# entry 59, which carries a 3 x 3 shared-start grid as well. Single-threaded
# R, no Stan sampling; 2 CPUs per task so each has a whole physical core.
#
# FOOTPRINT. 2 CPUs and 2G per task. Throttled to 24 at once (48 CPUs,
# 49,152 MB) because two pinned fits of 40 CPUs each (SLURM 35412, entries
# 60 and 61) are running: 80 + 48 = 128 CPUs, exactly the half-node budget.
# Raise the throttle only once those have finished.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/start-seq2.sh
#   SSQ2_MODE=merge Rscript --vanilla _scripts/start-seq2.R > _data/start-seq2-<date>.txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-91%24
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --job-name=start-seq2
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test
#SBATCH --output=_data/start-seq2/slurm-%a.out
#SBATCH --error=_data/start-seq2/slurm-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SSQ2_MODE=run Rscript --vanilla _scripts/start-seq2.R
