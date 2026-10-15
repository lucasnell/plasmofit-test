#!/bin/bash
# Merge the simulation study's chains into fits and write per-fit summaries,
# then run the reader. Submitted with --dependency=afterany:<chain array id>,
# so it runs once every chain task has ended; safe to rerun by hand at any
# time, since merge skips fits already merged and fits with missing chains.
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=simstudy-merge
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test
#SBATCH --output=_data/simstudy/logs/merge-%j.out
#SBATCH --error=_data/simstudy/logs/merge-%j.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

set -u
SIMSTUDY_MODE=merge /programs/R-4.6.1/bin/Rscript --vanilla _data/simstudy/simstudy-run.R
/programs/R-4.6.1/bin/Rscript --vanilla _scripts/simstudy-read.R \
    > "_data/simstudy-read-$(date +%F).txt"
