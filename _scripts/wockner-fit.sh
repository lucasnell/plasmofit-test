#!/bin/bash -l

# sbatch wrapper for wockner-fit.R, which until now existed only as a comment
# in that script's header. Runs with the working directory set to _data/, so
# the script's bare filenames read wockner-cleaned.csv from there and write
# wock-fit-<config>.rds beside the other saved fits -- no copying to a
# separate run directory, and outputs land where the _data/ naming table says.
#
# The array range must match the CONFIGS list in wockner-fit.R:
#   1 no_pool        2 pooled_cl     3 np_wide_prior   4 pl_wide_prior
#   5 np_center45    6 np_center42   7 np_wider_prior  8 np_anchor
#
#   sbatch --array=8 _scripts/wockner-fit.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-fit
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-fit-%a.out
#SBATCH --error=wock-fit-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla ../_scripts/wockner-fit.R
