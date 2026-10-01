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
#   9 np_wide_total0 10 np_wide_bshape 11 np_wide_both
#  12 pl_wide_total0 13 pl_wide_both
#  14 np_bs50        15 np_bs84       16 np_bs100
#  17 np_bs150       18 np_bs250      19 np_bs400
#  20 np_bs600       21 np_bs250_ms1000
#
# 14-18 are the b_shape ladder: b_shape pinned by a tight prior at five
# centres, all on the corrected log10_total0 prior, so they are read against
# task 9 (np_wide_total0) and against each other. See CONFIGS for why a
# ladder rather than one fixed value, and note they raise max_shape to 400.
#
#   sbatch --array=19-21 _scripts/wockner-fit.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=19-21
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
