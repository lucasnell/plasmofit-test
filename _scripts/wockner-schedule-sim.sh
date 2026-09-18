#!/bin/bash -l

# Schedule-bias simulation: 2 prior arms x 3 replicates. Submit from the repo
# root, on the cluster, with _data/ populated:
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim.sh
#
# Each task is one (arm, replicate) pair and takes ~2 h, the same cost as a
# real fit. See _scripts/wockner-schedule-sim.R for what it does and why.

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-schedsim
#SBATCH --output=_data/wock-schedsim-%a.out
#SBATCH --error=_data/wock-schedsim-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/wockner-schedule-sim.R
