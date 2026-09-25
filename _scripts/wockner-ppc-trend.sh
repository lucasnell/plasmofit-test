#!/bin/bash -l

# Option 1's posterior predictive check, as a detached job so it survives the
# session that submits it. See _scripts/wockner-ppc-trend.R for what it asks.
#
# Post-hoc on four saved fits, no refitting, but 200 draws x 177 series x 4
# fits of mat_exp_series is not fast -- allow a few hours.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-ppc-trend.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --job-name=wock-ppc-trend
#SBATCH --output=_data/wock-ppc-trend.out
#SBATCH --error=_data/wock-ppc-trend.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/wockner-ppc-trend.R > _data/ppc-trend.log 2>&1
