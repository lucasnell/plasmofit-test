#!/bin/bash -l

# Profile likelihood of cycle_length, 3 units x 2 b_shape modes; design and
# reading rule in _scripts/cl-profile.R. One CPU per task: submitted
# 2026-10-14 while SLURM 30858 held 120 of our 128 CPUs, so 6 x 1 CPU and
# 6 x 8G keep the total at 126 CPUs and ~147 GB, inside the half-node budget.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/cl-profile.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-6
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=cl-profile
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test
#SBATCH --output=_data/cl-profile-%a.out
#SBATCH --error=_data/cl-profile-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/cl-profile.R
